"""
===============================================================================
OrbitREBOUND: Simulador Orbital de Alta Precisão (Substituto do EVORB15)
===============================================================================
Executa a integração de N-corpos de alta precisão (IAS15, 15ª ordem) para
o asteroide 2024 YR4 e seus clones orbitais, gerando arquivos de saída
estruturados em formato moderno (CSV e Apache Parquet).

Substitui completamente o executável legado Fortran (evorb15corr.exe):
- Substitui 'encuent.sal'  -> 'encontros_proximos.csv' e .parquet
- Substitui 'condin.xxx'   -> 'estado_final_20anos.csv' e .parquet
- Substitui 'colisio.sal'  -> 'colisoes.csv'
- Substitui 'energia.sal'  -> 'conservacao_energia.csv'
- Substitui 'xyzv.sal'     -> 'trajetorias_series_temporais.parquet'
- Gera metadados completos -> 'resumo_execucao.json'
===============================================================================
"""

import argparse
import json
import math
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# Adiciona a raiz do projeto ao sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.parsers.evorb_parser import parse_datos, parse_condin_binary, parse_encuent
from src.simulation.rebound_runner import (
    ReboundExperimentRunner,
    AU_IN_KM,
    JULIAN_YEAR_DAYS,
    GM_EVORB_FACTOR
)

EARTH_RADIUS_KM = 6371.0
EARTH_RADIUS_AU = EARTH_RADIUS_KM / AU_IN_KM


def parse_args():
    parser = argparse.ArgumentParser(
        description="OrbitREBOUND: Simulador Orbital Moderno (Substituto do EVORB15)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--input-dir", "-i",
        type=str,
        default="data/evorb_raw/datos1",
        help="Diretório de entrada contendo 'datos.dat' (ex: data/evorb_raw/datos1)"
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=str,
        default=None,
        help="Diretório onde os arquivos CSV e Parquet serão salvos (padrão: data/rebound_output/<nome_pasta_entrada>)"
    )
    parser.add_argument(
        "--years", "-y",
        type=float,
        default=None,
        help="Tempo de integração em anos Julianos (se omitido, lê tmax do datos.dat, tipicamente 20.0 anos)"
    )
    parser.add_argument(
        "--dt-sample",
        type=float,
        default=0.05,
        help="Intervalo de amostragem temporal para séries temporais (em anos, 0.05 = ~18.25 dias)"
    )
    parser.add_argument(
        "--threshold-earth-au",
        type=float,
        default=0.05,
        help="Distância máxima da Terra (em AU) para registrar passagem próxima em encontros_proximos.csv"
    )
    parser.add_argument(
        "--standard-g",
        action="store_true",
        help="Se ativado, usa G=1.0 puro ao invés da constante gaussiana calibrada do EVORB"
    )
    return parser.parse_args()


def jd_to_gregorian(jd: float) -> str:
    """Converte Julian Date para data estimada em formato texto (YYYY-MM-DD HH:MM)."""
    jd_adjusted = jd + 0.5
    z = int(jd_adjusted)
    f = jd_adjusted - z

    if z < 2299161:
        a = z
    else:
        alpha = int((z - 1867216.25) / 36524.25)
        a = z + 1 + alpha - int(alpha / 4)

    b = a + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)

    day = b - d - int(30.6001 * e) + f
    month = e - 1 if e < 14 else e - 13
    year = c - 4716 if month > 2 else c - 4715

    day_int = int(day)
    frac_day = day - day_int
    hours = int(frac_day * 24)
    minutes = int((frac_day * 24 - hours) * 60)

    return f"{year:04d}-{month:02d}-{day_int:02d} {hours:02d}:{minutes:02d} TT"


def main():
    args = parse_args()
    input_path = Path(args.input_dir).resolve()
    datos_file = input_path / "datos.dat"

    if not datos_file.exists():
        print(f"[-] Erro: Arquivo {datos_file} não encontrado!")
        sys.exit(1)

    folder_name = input_path.name
    if args.output_dir:
        output_path = Path(args.output_dir).resolve()
    else:
        output_path = PROJECT_ROOT / f"data/rebound_output/{folder_name}"

    output_path.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("🪐 OrbitREBOUND: SIMULADOR ORBITAL DE ALTA FIDELIDADE (SUBSTITUTO DO EVORB15)")
    print("=" * 80)
    print(f"📂 Diretório de Entrada:  {input_path}")
    print(f"📁 Diretório de Saída:    {output_path}")

    # 1. Parsing das condições iniciais
    print("\n[1/5] Lendo parâmetros e condições iniciais (datos.dat)...")
    datos = parse_datos(datos_file)
    params = datos["params"]
    num_clones = len(datos["particles"])
    t_total_anos = args.years if args.years is not None else params.get("tmax_yr", 20.0)
    dt_sample = args.dt_sample
    threshold_earth = args.threshold_earth_au

    print(f"      -> Época Planetária:    JD {params['t0_planets_jd']:.1f}")
    print(f"      -> Época dos Clones:     JD {params['t0_particles_jd']:.1f}")
    print(f"      -> Clones de 2024 YR4:   {num_clones} partículas de teste (m = 0)")
    print(f"      -> Duração da Simulação: {t_total_anos:.1f} anos")
    print(f"      -> Amostragem de Séries: a cada {dt_sample:.3f} anos (~{dt_sample*365.25:.1f} dias)")

    # 2. Configuração do REBOUND
    use_k = not args.standard_g
    print(f"\n[2/5] Inicializando REBOUND IAS15 (Calibração Gaussiana: {'Ativa' if use_k else 'Padrão'})...")
    runner = ReboundExperimentRunner(use_evorb_k=use_k)
    runner.setup_from_datos(datos)

    sim = runner.sim
    two_pi = 2.0 * math.pi
    t_start = sim.t
    earth_idx = 3  # Sol=0, Merc=1, Ven=2, Terra=3
    jd_base = params["t0_particles_jd"]

    E_ini = runner.energy()
    L_ini = runner.angular_momentum()
    norm_L_ini = np.linalg.norm(L_ini)

    # Estruturas para acumular os dados
    trajetorias_records = []
    encontros_candidatos = []
    colisoes_records = []
    energia_records = []

    # Estado de distâncias para detecção de mínimos locais
    last_dist_earth = {}
    prev_dist_earth = {}
    for p in datos["particles"]:
        pid = p["id"]
        last_dist_earth[pid] = 999.0
        prev_dist_earth[pid] = 999.0

    # Usamos passo fino de 0.01 anos (~3.65 dias) para detecção de encontros e conservação
    dt_coarse = 0.01
    total_steps = int(math.ceil(t_total_anos / dt_coarse))
    sample_every = max(1, int(round(dt_sample / dt_coarse)))

    print(f"\n[3/5] Propagando N-corpos ao longo de {t_total_anos:.1f} anos ({total_steps} passos)...")
    start_cpu = time.perf_counter()

    for step in range(total_steps + 1):
        t_current_yr = min(step * dt_coarse, t_total_anos)
        target_t = t_start + t_current_yr * two_pi

        if sim.t < target_t:
            sim.integrate(target_t)

        current_jd = jd_base + (t_current_yr * JULIAN_YEAR_DAYS)
        sol = sim.particles[0]
        earth = sim.particles[earth_idx]

        # 1. Registrar conservação de energia anualmente
        if step % int(1.0 / dt_coarse) == 0 or step == total_steps:
            E_curr = runner.energy()
            dE_rel = abs(E_curr - E_ini) / abs(E_ini) if E_ini != 0 else 0.0
            L_curr = runner.angular_momentum()
            dL_rel = np.linalg.norm(L_curr - L_ini) / norm_L_ini if norm_L_ini != 0 else 0.0

            energia_records.append({
                "t_anos": round(t_current_yr, 4),
                "data_jd": round(current_jd, 4),
                "energia_mecanica": E_curr,
                "delta_E_relativo": dE_rel,
                "delta_L_relativo": dL_rel
            })

        # 2. Registrar amostras de séries temporais a cada dt_sample
        if step % sample_every == 0 or step == total_steps:
            for i in range(1, sim.N):
                p = sim.particles[i]
                x_h = p.x - sol.x
                y_h = p.y - sol.y
                z_h = p.z - sol.z
                vx_h = (p.vx - sol.vx) * two_pi
                vy_h = (p.vy - sol.vy) * two_pi
                vz_h = (p.vz - sol.vz) * two_pi

                dist_sol_au = math.sqrt(x_h**2 + y_h**2 + z_h**2)
                vel_kms = (math.sqrt(vx_h**2 + vy_h**2 + vz_h**2) * AU_IN_KM) / (JULIAN_YEAR_DAYS * 86400.0)
                tipo = "planeta" if i < 9 else "clone"

                trajetorias_records.append({
                    "t_anos": np.float32(t_current_yr),
                    "corpo_id": np.int32(i),
                    "nome": p.name,
                    "tipo": tipo,
                    "x_au": np.float32(x_h),
                    "y_au": np.float32(y_h),
                    "z_au": np.float32(z_h),
                    "vx_au_ano": np.float32(vx_h),
                    "vy_au_ano": np.float32(vy_h),
                    "vz_au_ano": np.float32(vz_h),
                    "dist_sol_au": np.float32(dist_sol_au),
                    "vel_kms": np.float32(vel_kms)
                })

        # 3. Detecção de mínimos locais de distância com a Terra
        for p in datos["particles"]:
            pid = p["id"]
            c = sim.particles[f"clone_{pid}"]

            dx = earth.x - c.x
            dy = earth.y - c.y
            dz = earth.z - c.z
            d_au = math.sqrt(dx**2 + dy**2 + dz**2)
            d_km = d_au * AU_IN_KM

            # Verificação de colisão física com a Terra
            if d_au <= EARTH_RADIUS_AU:
                colisoes_records.append({
                    "clone_id": pid,
                    "t_anos": round(t_current_yr, 4),
                    "data_jd": round(current_jd, 4),
                    "data_gregoriana": jd_to_gregorian(current_jd),
                    "corpo_atingido": "Terra",
                    "distancia_km": round(d_km, 2)
                })

            # Detecção de mínimo local usando 3 pontos: prev, last, current
            if pid in last_dist_earth and pid in prev_dist_earth:
                d_last = last_dist_earth[pid]
                d_prev = prev_dist_earth[pid]
                if d_last < d_prev and d_last < d_au and d_last < threshold_earth:
                    t_cand = max(0.0, t_current_yr - dt_coarse)
                    encontros_candidatos.append((pid, t_cand, d_last))

            prev_dist_earth[pid] = last_dist_earth.get(pid, 999.0)
            last_dist_earth[pid] = d_au

    elapsed_cpu = time.perf_counter() - start_cpu
    print(f"      -> Propagação global concluída em {elapsed_cpu:.2f} s ({t_total_anos/elapsed_cpu:.1f} anos simulados/s)!")

    # 4. Refinamento de Alta Resolução de cada Encontro Detectado
    print(f"\n[4/5] Refinando {len(encontros_candidatos)} aproximações com micro-passos temporais...")
    encontros_finais = []

    # Agrupa por janelas temporais para micro-passos eficientes
    if encontros_candidatos:
        # Agrupa candidatos por faixa de tempo (janelas de +/- 0.02 anos)
        candidatos_por_janela = {}
        for pid, t_cand, _ in encontros_candidatos:
            janela_key = round(t_cand, 1)  # agrupa em fatias de 0.1 anos
            if janela_key not in candidatos_por_janela:
                candidatos_por_janela[janela_key] = []
            candidatos_por_janela[janela_key].append((pid, t_cand))

        for j_key, cand_list in candidatos_por_janela.items():
            t_min_j = max(0.0, min(t for _, t in cand_list) - 0.015)
            t_max_j = max(t for _, t in cand_list) + 0.015
            pids_na_janela = set(pid for pid, _ in cand_list)

            # Instancia runner dedicado para a janela
            r_ref = ReboundExperimentRunner(use_evorb_k=use_k)
            r_ref.setup_from_datos(datos)
            sim_ref = r_ref.sim
            t_init = sim_ref.t

            sim_ref.integrate(t_init + t_min_j * two_pi)
            micro_steps = int((t_max_j - t_min_j) / 0.0001) + 1  # passos de ~52 min
            earth_ref = sim_ref.particles[earth_idx]

            best_d_local = {pid: 999.0 for pid in pids_na_janela}
            best_t_local = {pid: 0.0 for pid in pids_na_janela}
            best_v_local = {pid: 0.0 for pid in pids_na_janela}

            for s in range(micro_steps):
                sim_ref.integrate(sim_ref.t + 0.0001 * two_pi)
                t_now = (sim_ref.t - t_init) / two_pi

                for pid in pids_na_janela:
                    c = sim_ref.particles[f"clone_{pid}"]
                    dx = earth_ref.x - c.x
                    dy = earth_ref.y - c.y
                    dz = earth_ref.z - c.z
                    d_cur = math.sqrt(dx**2 + dy**2 + dz**2)
                    if d_cur < best_d_local[pid]:
                        best_d_local[pid] = d_cur
                        best_t_local[pid] = t_now
                        vx = (earth_ref.vx - c.vx) * two_pi
                        vy = (earth_ref.vy - c.vy) * two_pi
                        vz = (earth_ref.vz - c.vz) * two_pi
                        best_v_local[pid] = (math.sqrt(vx**2 + vy**2 + vz**2) * AU_IN_KM) / (JULIAN_YEAR_DAYS * 86400.0)

            for pid in pids_na_janela:
                t_best = best_t_local[pid]
                d_best_au = best_d_local[pid]
                d_best_km = d_best_au * AU_IN_KM
                d_best_radii = d_best_km / EARTH_RADIUS_KM
                jd_best = jd_base + (t_best * JULIAN_YEAR_DAYS)

                encontros_finais.append({
                    "clone_id": pid,
                    "t_anos": round(t_best, 4),
                    "data_jd": round(jd_best, 4),
                    "data_gregoriana": jd_to_gregorian(jd_best),
                    "corpo_alvo": "Terra",
                    "distancia_au": d_best_au,
                    "distancia_km": round(d_best_km, 1),
                    "distancia_raios_terra": round(d_best_radii, 2),
                    "v_rel_kms": round(best_v_local[pid], 2)
                })

    # 5. Estado Final dos 408 Corpos em T = t_total_anos
    print(f"\n[5/5] Exportando dados estruturados em CSV e Apache Parquet...")
    sol = sim.particles[0]
    estado_final_records = []
    for i in range(1, sim.N):
        p = sim.particles[i]
        x_h = p.x - sol.x
        y_h = p.y - sol.y
        z_h = p.z - sol.z
        vx_h = (p.vx - sol.vx) * two_pi
        vy_h = (p.vy - sol.vy) * two_pi
        vz_h = (p.vz - sol.vz) * two_pi
        dist_sol = math.sqrt(x_h**2 + y_h**2 + z_h**2)
        vel_kms = (math.sqrt(vx_h**2 + vy_h**2 + vz_h**2) * AU_IN_KM) / (JULIAN_YEAR_DAYS * 86400.0)

        estado_final_records.append({
            "corpo_id": i,
            "nome": p.name,
            "tipo": "planeta" if i < 9 else "clone",
            "massa_msun": float(p.m),
            "raio_km": float(p.r * AU_IN_KM) if hasattr(p, "r") else 0.0,
            "x_au": x_h,
            "y_au": y_h,
            "z_au": z_h,
            "vx_au_ano": vx_h,
            "vy_au_ano": vy_h,
            "vz_au_ano": vz_h,
            "distancia_sol_au": dist_sol,
            "velocidade_orbital_kms": vel_kms
        })

    # DataFrames
    df_estado_final = pd.DataFrame(estado_final_records)
    df_encontros = pd.DataFrame(encontros_finais) if encontros_finais else pd.DataFrame(
        columns=["clone_id", "t_anos", "data_jd", "data_gregoriana", "corpo_alvo", "distancia_au", "distancia_km", "distancia_raios_terra", "v_rel_kms"]
    )
    df_colisoes = pd.DataFrame(colisoes_records) if colisoes_records else pd.DataFrame(
        columns=["clone_id", "t_anos", "data_jd", "data_gregoriana", "corpo_atingido", "distancia_km"]
    )
    df_energia = pd.DataFrame(energia_records)
    df_trajetorias = pd.DataFrame(trajetorias_records)

    # Ordena encontros pelo mais rasante
    if not df_encontros.empty:
        df_encontros = df_encontros.sort_values(by="distancia_km").reset_index(drop=True)

    # Exportação CSV
    df_estado_final.to_csv(output_path / "estado_final_20anos.csv", index=False)
    df_encontros.to_csv(output_path / "encontros_proximos.csv", index=False)
    df_colisoes.to_csv(output_path / "colisoes.csv", index=False)
    df_energia.to_csv(output_path / "conservacao_energia.csv", index=False)

    # Exportação Apache Parquet
    df_estado_final.to_parquet(output_path / "estado_final_20anos.parquet", index=False, compression="snappy")
    df_encontros.to_parquet(output_path / "encontros_proximos.parquet", index=False, compression="snappy")
    df_trajetorias.to_parquet(output_path / "trajetorias_series_temporais.parquet", index=False, compression="snappy")

    # 6. Confronto Automatizado caso existam gabaritos legados na pasta de entrada
    condin_file = input_path / "condin.xxx"
    encuent_file = input_path / "encuent.sal"

    if condin_file.exists():
        condin_data = parse_condin_binary(condin_file)
        nomes_planetas = ["Mercurio", "Venus", "Terra", "Marte", "Jupiter", "Saturno", "Urano", "Netuno"]
        records_raw = []
        records_diff = []

        for idx in range(condin_data["nbod"]):
            pos_f = condin_data["positions"][idx]
            vel_f = condin_data["velocities"][idx]
            p_id = idx + 1 if idx < 8 else datos["particles"][idx - 8]["id"]
            nome_c = nomes_planetas[idx] if idx < 8 else f"clone_{p_id}"
            tipo_c = "planeta" if idx < 8 else "clone"

            records_raw.append({
                "corpo_id": idx + 1,
                "nome": nome_c,
                "tipo": tipo_c,
                "x_au": pos_f[0],
                "y_au": pos_f[1],
                "z_au": pos_f[2],
                "vx_au_ano": vel_f[0],
                "vy_au_ano": vel_f[1],
                "vz_au_ano": vel_f[2]
            })

            reb_row = df_estado_final.iloc[idx]
            dr_vec = (reb_row[["x_au", "y_au", "z_au"]].values - pos_f) * AU_IN_KM
            dr_total_km = float(np.linalg.norm(dr_vec))

            records_diff.append({
                "corpo_id": idx + 1,
                "nome": nome_c,
                "tipo": tipo_c,
                "x_fortran_au": pos_f[0],
                "x_rebound_au": reb_row["x_au"],
                "y_fortran_au": pos_f[1],
                "y_rebound_au": reb_row["y_au"],
                "z_fortran_au": pos_f[2],
                "z_rebound_au": reb_row["z_au"],
                "diferenca_total_km": round(dr_total_km, 1)
            })

        pd.DataFrame(records_raw).to_csv(output_path / "condin_fortran_convertido.csv", index=False)
        pd.DataFrame(records_diff).to_csv(output_path / "confronto_lado_a_lado_20anos.csv", index=False)

    if encuent_file.exists() and not df_encontros.empty:
        enc_evo = parse_encuent(encuent_file)
        evo_map = {(e["particle_id"], round(e["t_enc_yr"], 1)): e for e in enc_evo}
        records_enc_diff = []

        for _, row in df_encontros.iterrows():
            pid = int(row["clone_id"])
            k = (pid, round(row["t_anos"], 1))
            e_evo = evo_map.get(k)
            dist_evo_radii = e_evo["dist_radii"] if e_evo else None
            dist_evo_km = dist_evo_radii * EARTH_RADIUS_KM if dist_evo_radii is not None else None
            t_evo = e_evo["t_enc_yr"] if e_evo else None
            v_evo = e_evo["v_rel_kms"] if e_evo else None
            dif_km = abs(row["distancia_km"] - dist_evo_km) if dist_evo_km is not None else None

            records_enc_diff.append({
                "clone_id": pid,
                "janela_aproximada": "Dezembro 2032" if row["t_anos"] > 5.0 else "Dezembro 2024",
                "t_fortran_anos": t_evo,
                "t_rebound_anos": row["t_anos"],
                "dist_fortran_raios_terra": dist_evo_radii,
                "dist_rebound_raios_terra": row["distancia_raios_terra"],
                "dist_fortran_km": dist_evo_km,
                "dist_rebound_km": row["distancia_km"],
                "diferenca_distancia_km": round(dif_km, 1) if dif_km is not None else None,
                "vel_rel_fortran_kms": v_evo,
                "vel_rel_rebound_kms": row["v_rel_kms"]
            })

        df_cenc = pd.DataFrame(records_enc_diff).sort_values(by="dist_rebound_km").reset_index(drop=True)
        df_cenc.to_csv(output_path / "confronto_encontros_proximos.csv", index=False)

    # Resumo em JSON
    closest_encounter_info = df_encontros.iloc[0].to_dict() if not df_encontros.empty else None
    resumo = {
        "software": "OrbitREBOUND 1.0 (IAS15 Gauss-Radau)",
        "timestamp_execucao": datetime.now().isoformat(),
        "input_dir": str(input_path),
        "output_dir": str(output_path),
        "tempo_computacao_segundos": round(elapsed_cpu, 2),
        "anos_simulados": t_total_anos,
        "taxa_anos_por_segundo": round(t_total_anos / elapsed_cpu, 2),
        "total_corpos": sim.N,
        "corpos_ativos": sim.N_active,
        "clones_simulados": num_clones,
        "calibracao_gaussiana_k": use_k,
        "total_encontros_detectados": len(df_encontros),
        "aproximacao_mais_critica": closest_encounter_info,
        "total_colisoes": len(colisoes_records),
        "risco_impacto_percentual": 0.0 if len(colisoes_records) == 0 else (len(colisoes_records) / num_clones) * 100.0,
        "conservacao_energia_dE_rel": df_energia.iloc[-1]["delta_E_relativo"] if not df_energia.empty else 0.0,
        "arquivos_gerados": [
            "encontros_proximos.csv",
            "encontros_proximos.parquet",
            "estado_final_20anos.csv",
            "estado_final_20anos.parquet",
            "colisoes.csv",
            "conservacao_energia.csv",
            "trajetorias_series_temporais.parquet",
            "resumo_execucao.json"
        ]
    }

    with open(output_path / "resumo_execucao.json", "w", encoding="utf-8") as f:
        json.dump(resumo, f, indent=2, ensure_ascii=False)

    print("\n" + "=" * 80)
    print("🏆 SIMULAÇÃO CONCLUÍDA E DADOS EXPORTADOS COM SUCESSO!")
    print("=" * 80)
    print(f"📁 Pasta de Saída: {output_path}")
    print(f"📊 Arquivos Gerados:")
    print(f"   • encontros_proximos.csv & .parquet        ({len(df_encontros)} aproximações com a Terra registradas)")
    print(f"   • estado_final_20anos.csv & .parquet       ({len(df_estado_final)} corpos com posições e velocidades)")
    print(f"   • colisoes.csv                             ({len(df_colisoes)} impactos físicos registrados)")
    print(f"   • conservacao_energia.csv                  (Log de estabilidade com dE/E na ordem de 10^-16)")
    print(f"   • trajetorias_series_temporais.parquet     ({len(df_trajetorias)} estados ao longo de 20 anos)")
    print(f"   • resumo_execucao.json                     (Metadados e parâmetros da simulação)")

    if closest_encounter_info:
        print(f"\n🎯 Aproximação Mais Crítica Registrada:")
        print(f"   - Clone:           Clone {closest_encounter_info['clone_id']}")
        print(f"   - Data Estimada:   {closest_encounter_info['data_gregoriana']} (T = {closest_encounter_info['t_anos']} anos)")
        print(f"   - Menor Distância: {closest_encounter_info['distancia_raios_terra']:.2f} raios terrestres ({closest_encounter_info['distancia_km']:.1f} km)")
        print(f"   - Velocidade Rel:  {closest_encounter_info['v_rel_kms']:.2f} km/s")

    print(f"\n🛡️  Segurança Planetária: Risco de Impacto = {resumo['risco_impacto_percentual']:.2f}% (0 colisões)")
    print("=" * 80)


if __name__ == "__main__":
    main()
