"""
===============================================================================
OrbitREBOUND: Simulador Orbital (Substituto Direto do EVORB15 em Python/C)
===============================================================================
Lê as condições iniciais de 'datos.dat' e executa a integração orbital de N-corpos
no REBOUND com o solver Gauss-Radau IAS15 de 15ª ordem.

Organiza as saídas do EVORB15 em duas pastas dedicadas:
1. 'final/'    -> Recorte consolidado no tempo final (T = 20 anos) dos 7 arquivos:
                  - orbitas.csv  (elementos orbitais finais)
                  - condin.csv   (estado cartesiano final, espelho do condin.xxx)
                  - encuent.csv  (encontros ordenados por menor distância)
                  - colisio.csv  (resumo consolidado de colisões)
                  - finales.csv  (resumo consolidado de ejeções)
                  - energia.csv  (balanço final de conservação de energia)
                  - xyzv.csv     (snapshot cartesiano final)

2. 'temporal/' -> Série temporal completa (T = 0 até T = T_final) dos 7 arquivos:
                  - orbitas.csv  (evolução temporal orbital completa desde t=0)
                  - condin.csv   (vetores cartesianos de re-arranque no tempo)
                  - encuent.csv  (cronograma de passagens rasantes por data)
                  - colisio.csv  (linha do tempo de impactos)
                  - finales.csv  (linha do tempo de ejeções)
                  - energia.csv  (curva temporal de energia cinética, potencial e total)
                  - xyzv.csv     (histórico completo de posições e velocidades)
===============================================================================
"""

import argparse
import math
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd

# Adiciona raiz do projeto
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.parsers.evorb_parser import parse_datos
from src.simulation.rebound_runner import (
    ReboundExperimentRunner,
    AU_IN_KM,
    JULIAN_YEAR_DAYS
)

EARTH_RADIUS_KM = 6371.0


def parse_args():
    parser = argparse.ArgumentParser(
        description="OrbitREBOUND: Substituto Moderno do EVORB15",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--input-dir", "-i",
        type=str,
        default="data/evorb_raw/datos1",
        help="Pasta contendo o arquivo datos.dat"
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=str,
        default=None,
        help="Pasta onde as pastas 'final/' e 'temporal/' serão salvas"
    )
    parser.add_argument(
        "--years", "-y",
        type=float,
        default=None,
        help="Anos de simulação (padrão: lê do datos.dat, tipicamente 20.0 anos)"
    )
    parser.add_argument(
        "--dt-temporal",
        type=float,
        default=None,
        help="Intervalo em anos para amostragem na pasta 'temporal/' (padrão: lê exatamente o 'tsal' do datos.dat, tipicamente 0.001 ano = ~8.76 horas)"
    )
    return parser.parse_args()


def jd_to_gregorian(jd: float) -> str:
    """Converte Julian Date para data estimada YYYY-MM-DD."""
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

    return f"{year:04d}-{month:02d}-{int(day):02d}"


def main():
    args = parse_args()
    input_path = Path(args.input_dir).resolve()
    datos_file = input_path / "datos.dat"

    if not datos_file.exists():
        print(f"[-] Erro: Arquivo {datos_file} não encontrado!")
        sys.exit(1)

    folder_name = input_path.name
    output_path = Path(args.output_dir).resolve() if args.output_dir else PROJECT_ROOT / f"data/rebound_output/{folder_name}"
    
    final_dir = output_path / "final"
    temporal_dir = output_path / "temporal"
    final_dir.mkdir(parents=True, exist_ok=True)
    temporal_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("🪐 OrbitREBOUND: SIMULADOR ORBITAL (SUBSTITUTO DIRETO DO EVORB15)")
    print("=" * 80)
    print(f"📂 Entrada:    {input_path}")
    print(f"📁 Pasta Final:    {final_dir}")
    print(f"📁 Pasta Temporal: {temporal_dir}")

    # 1. Carrega dados.dat
    print("\n[1/4] Lendo dados.dat...")
    datos = parse_datos(datos_file)
    params = datos["params"]
    num_clones = len(datos["particles"])
    t_total_anos = args.years if args.years is not None else params.get("tmax_yr", 20.0)
    r_crit_min = params.get("rcrit1_au", 0.1)
    r_crit_max = params.get("rcrit2_au", 100.0)
    dt_temp = args.dt_temporal if args.dt_temporal is not None else params.get("tsal_yr", 0.001)

    # 2. Configura REBOUND com calibração gaussiana idêntica ao EVORB
    print(f"[2/4] Configurando REBOUND (IAS15, 15ª ordem, 8 planetas + {num_clones} clones)...")
    runner = ReboundExperimentRunner(use_evorb_k=True)
    runner.setup_from_datos(datos)

    sim = runner.sim
    two_pi = 2.0 * math.pi
    t_start = sim.t
    earth_idx = 3
    jd_base = params["t0_particles_jd"]
    nomes_planetas = ["Mercurio", "Venus", "Terra", "Marte", "Jupiter", "Saturno", "Urano", "Netuno"]

    # 3. Prepara arquivos de streaming para a pasta temporal/
    f_orb = open(temporal_dir / "orbitas.csv", "w", encoding="utf-8", buffering=1024 * 1024)
    f_orb.write("t_anos,N,nome,tipo,a,e,i_graus,nodo_graus,argper_graus,anomedia_graus\n")

    f_xyz = open(temporal_dir / "xyzv.csv", "w", encoding="utf-8", buffering=1024 * 1024)
    f_xyz.write("t_anos,N,nome,tipo,x_au,y_au,z_au,vx_au_ano,vy_au_ano,vz_au_ano\n")

    f_con = open(temporal_dir / "condin.csv", "w", encoding="utf-8", buffering=1024 * 1024)
    f_con.write("t_anos,N,nome,tipo,x_au,y_au,z_au,vx_au_ano,vy_au_ano,vz_au_ano\n")

    f_ene = open(temporal_dir / "energia.csv", "w", encoding="utf-8", buffering=64 * 1024)
    f_ene.write("t_anos,energia_total,energia_cinetica,energia_potencial,delta_e_relativo\n")

    # Fator de conversão da constante de Gauss para unidades físicas (Massa Solar * AU^2 / ano^2)
    # No Fortran EVORB: GM = (k * 365.25)^2 ≈ 39.476926
    gm_factor = (0.01720209895 * 365.25)**2

    def grava_snapshot_temporal(t_val, e_base):
        sol_now = sim.particles[0]
        # Energia em unidades físicas idênticas ao EVORB (AU^2 / ano^2)
        e_now = sim.energy()
        t_kin = 0.5 * sum(sim.particles[k].m * (sim.particles[k].vx**2 + sim.particles[k].vy**2 + sim.particles[k].vz**2) for k in range(9))
        de_rel = abs((e_now - e_base) / e_base) if abs(e_base) > 1e-15 else 0.0
        f_ene.write(f"{t_val:.4f},{e_now * gm_factor:.12f},{t_kin * gm_factor:.12f},{(e_now - t_kin) * gm_factor:.12f},{de_rel:.6e}\n")

        # Orbitas e Coordenadas para cada corpo
        for k in range(1, sim.N):
            pk = sim.particles[k]
            nome_k = nomes_planetas[k - 1] if k < 9 else f"clone_{pk.name.split('_')[-1]}"
            tipo_k = "planeta" if k < 9 else "clone"

            # Orbitais
            orb = pk.orbit(primary=sol_now)
            f_orb.write(f"{t_val:.4f},{k},{nome_k},{tipo_k},{orb.a:.6f},{orb.e:.6f},{np.degrees(orb.inc):.4f},{np.degrees(orb.Omega):.4f},{np.degrees(orb.omega):.4f},{np.degrees(orb.M):.4f}\n")

            # Cartesianos
            xh = pk.x - sol_now.x
            yh = pk.y - sol_now.y
            zh = pk.z - sol_now.z
            vxh = (pk.vx - sol_now.vx) * two_pi
            vyh = (pk.vy - sol_now.vy) * two_pi
            vzh = (pk.vz - sol_now.vz) * two_pi
            linha_cart = f"{t_val:.4f},{k},{nome_k},{tipo_k},{xh:.6f},{yh:.6f},{zh:.6f},{vxh:.6f},{vyh:.6f},{vzh:.6f}\n"
            f_xyz.write(linha_cart)
            f_con.write(linha_cart)

    # 4. Propagação temporal com REBOUND
    print(f"[3/4] Integrando {t_total_anos:.1f} anos no futuro (dt_temporal = {dt_temp:.4f} anos)...")
    start_time = time.perf_counter()

    e_ini = sim.energy()
    # Grava estado inicial t = 0.0 na pasta temporal/
    grava_snapshot_temporal(0.0, e_ini)

    dt_step = dt_temp
    total_steps = int(round(t_total_anos / dt_step))
    next_save_t = dt_temp

    last_d = {p["id"]: 999.0 for p in datos["particles"]}
    prev_d = {p["id"]: 999.0 for p in datos["particles"]}
    candidatos_enc = []
    colisoes_list = []
    finales_list = []

    for s in range(1, total_steps + 1):
        t_curr = min(s * dt_step, t_total_anos)
        sim.integrate(t_start + t_curr * two_pi)
        sol = sim.particles[0]
        earth = sim.particles[earth_idx]

        # Amostragem periódica na pasta temporal/
        if t_curr >= (next_save_t - 1e-6) or s == total_steps:
            grava_snapshot_temporal(t_curr, e_ini)
            next_save_t += dt_temp

        # Checagens para cada clone (colisões, limites e encontros)
        for p in datos["particles"]:
            pid = p["id"]
            c = sim.particles[f"clone_{pid}"]
            dx = earth.x - c.x
            dy = earth.y - c.y
            dz = earth.z - c.z
            d_au = math.sqrt(dx**2 + dy**2 + dz**2)

            # Colisão física com a Terra
            if d_au <= (EARTH_RADIUS_KM / AU_IN_KM):
                colisoes_list.append({
                    "t_anos": round(t_curr, 4),
                    "clone_id": pid,
                    "corpo_atingido": "Terra",
                    "dist_km": round(d_au * AU_IN_KM, 1),
                    "vel_impacto_kms": round(math.sqrt((earth.vx - c.vx)**2 + (earth.vy - c.vy)**2 + (earth.vz - c.vz)**2) * two_pi * AU_IN_KM / (JULIAN_YEAR_DAYS * 86400.0), 2)
                })

            # Checagem de limites de ejeção / colisão solar (finales)
            r_sol = math.sqrt((c.x - sol.x)**2 + (c.y - sol.y)**2 + (c.z - sol.z)**2)
            if r_sol < r_crit_min:
                finales_list.append({
                    "t_anos": round(t_curr, 4),
                    "clone_id": pid,
                    "motivo_eliminacao": "queda_solar",
                    "dist_au": round(r_sol, 4)
                })
            elif r_sol > r_crit_max:
                finales_list.append({
                    "t_anos": round(t_curr, 4),
                    "clone_id": pid,
                    "motivo_eliminacao": "ejecao_sistema",
                    "dist_au": round(r_sol, 4)
                })

            # Detecção de mínimo local de aproximação (< 0.05 AU)
            if pid in last_d and pid in prev_d:
                if last_d[pid] < prev_d[pid] and last_d[pid] < d_au and last_d[pid] < 0.05:
                    candidatos_enc.append((pid, t_curr - dt_step))

            prev_d[pid] = last_d[pid]
            last_d[pid] = d_au

    # Fecha arquivos de streaming temporal
    f_orb.close()
    f_xyz.close()
    f_con.close()
    f_ene.close()

    # Refinamento de precisão nos encontros detectados
    encontros_list = []
    if candidatos_enc:
        janelas = {}
        for pid, t_cand in candidatos_enc:
            j_key = round(t_cand, 1)
            if j_key not in janelas:
                janelas[j_key] = []
            janelas[j_key].append(pid)

        for j_key, pids in janelas.items():
            t_min_j = max(0.0, j_key - 0.03)
            t_max_j = j_key + 0.03
            pids_set = set(pids)

            r_temp = ReboundExperimentRunner(use_evorb_k=True)
            r_temp.setup_from_datos(datos)
            sim_t = r_temp.sim
            t0_t = sim_t.t

            sim_t.integrate(t0_t + t_min_j * two_pi)
            micro_steps = int((t_max_j - t_min_j) / 0.0001) + 1
            e_t = sim_t.particles[earth_idx]

            best_d = {pid: 999.0 for pid in pids_set}
            best_t = {pid: 0.0 for pid in pids_set}
            best_v = {pid: 0.0 for pid in pids_set}

            for _ in range(micro_steps):
                sim_t.integrate(sim_t.t + 0.0001 * two_pi)
                t_now = (sim_t.t - t0_t) / two_pi
                for pid in pids_set:
                    c = sim_t.particles[f"clone_{pid}"]
                    d = math.sqrt((e_t.x - c.x)**2 + (e_t.y - c.y)**2 + (e_t.z - c.z)**2)
                    if d < best_d[pid]:
                        best_d[pid] = d
                        best_t[pid] = t_now
                        vx = (e_t.vx - c.vx) * two_pi
                        vy = (e_t.vy - c.vy) * two_pi
                        vz = (e_t.vz - c.vz) * two_pi
                        best_v[pid] = (math.sqrt(vx**2 + vy**2 + vz**2) * AU_IN_KM) / (JULIAN_YEAR_DAYS * 86400.0)

            for pid in pids_set:
                d_km = best_d[pid] * AU_IN_KM
                d_radii = d_km / EARTH_RADIUS_KM
                jd_enc = jd_base + (best_t[pid] * JULIAN_YEAR_DAYS)
                encontros_list.append({
                    "t_anos": round(best_t[pid], 4),
                    "clone_id": pid,
                    "data_estimada": jd_to_gregorian(jd_enc),
                    "dist_raios_terra": round(d_radii, 1),
                    "dist_km": round(d_km, 1),
                    "vel_rel_kms": round(best_v[pid], 2)
                })

    elapsed = time.perf_counter() - start_time
    print(f"      -> Concluído em {elapsed:.2f} s ({t_total_anos/elapsed:.1f} anos/s)!")

    # 5. Salva os arquivos restantes de temporal/ e todos os de final/
    print("\n[4/4] Gerando os arquivos nas pastas final/ e temporal/...")

    # A) Encontros (temporal = cronológico; final = ordenado por dist_km)
    df_enc_cron = pd.DataFrame(encontros_list) if encontros_list else pd.DataFrame(
        columns=["t_anos", "clone_id", "data_estimada", "dist_raios_terra", "dist_km", "vel_rel_kms"]
    )
    if not df_enc_cron.empty:
        df_enc_cron = df_enc_cron.sort_values(by=["t_anos", "clone_id"]).reset_index(drop=True)
    df_enc_cron.to_csv(temporal_dir / "encuent.csv", index=False)

    df_enc_final = df_enc_cron.copy()
    if not df_enc_final.empty:
        df_enc_final = df_enc_final.sort_values(by="dist_km").reset_index(drop=True)
    df_enc_final.to_csv(final_dir / "encuent.csv", index=False)

    # B) Colisões (temporal = cronológico; final = consolidado)
    df_col_cron = pd.DataFrame(colisoes_list) if colisoes_list else pd.DataFrame(
        columns=["t_anos", "clone_id", "corpo_atingido", "dist_km", "vel_impacto_kms"]
    )
    df_col_cron.to_csv(temporal_dir / "colisio.csv", index=False)
    df_col_cron.to_csv(final_dir / "colisio.csv", index=False)

    # C) Finais / Ejeções (temporal = cronológico; final = consolidado)
    df_fin_cron = pd.DataFrame(finales_list) if finales_list else pd.DataFrame(
        columns=["t_anos", "clone_id", "motivo_eliminacao", "dist_au"]
    )
    df_fin_cron.to_csv(temporal_dir / "finales.csv", index=False)
    df_fin_cron.to_csv(final_dir / "finales.csv", index=False)

    # D) PASTA FINAL/ (Snapshots consolidados dos 408 corpos no instante final t = 20.0)
    sol_final = sim.particles[0]
    final_orbitas_records = []
    final_condin_records = []

    for i in range(1, sim.N):
        p = sim.particles[i]
        orb = p.orbit(primary=sol_final)
        nome = nomes_planetas[i - 1] if i < 9 else f"clone_{p.name.split('_')[-1]}"
        tipo = "planeta" if i < 9 else "clone"

        # Orbitais
        final_orbitas_records.append({
            "N": i,
            "nome": nome,
            "tipo": tipo,
            "a": round(float(orb.a), 6),
            "e": round(float(orb.e), 6),
            "i_graus": round(float(np.degrees(orb.inc)), 4),
            "nodo_graus": round(float(np.degrees(orb.Omega)), 4),
            "argper_graus": round(float(np.degrees(orb.omega)), 4),
            "anomedia_graus": round(float(np.degrees(orb.M)), 4)
        })

        # Cartesianos (X, Y, Z, Vx, Vy, Vz)
        xh = p.x - sol_final.x
        yh = p.y - sol_final.y
        zh = p.z - sol_final.z
        vxh = (p.vx - sol_final.vx) * two_pi
        vyh = (p.vy - sol_final.vy) * two_pi
        vzh = (p.vz - sol_final.vz) * two_pi
        final_condin_records.append({
            "N": i,
            "nome": nome,
            "tipo": tipo,
            "x_au": round(float(xh), 6),
            "y_au": round(float(yh), 6),
            "z_au": round(float(zh), 6),
            "vx_au_ano": round(float(vxh), 6),
            "vy_au_ano": round(float(vyh), 6),
            "vz_au_ano": round(float(vzh), 6)
        })

    pd.DataFrame(final_orbitas_records).to_csv(final_dir / "orbitas.csv", index=False)
    df_final_condin = pd.DataFrame(final_condin_records)
    df_final_condin.to_csv(final_dir / "condin.csv", index=False)
    df_final_condin.to_csv(final_dir / "xyzv.csv", index=False)

    # Energia final (resumo)
    e_end = sim.energy()
    t_kin_end = 0.5 * sum(sim.particles[k].m * (sim.particles[k].vx**2 + sim.particles[k].vy**2 + sim.particles[k].vz**2) for k in range(9))
    de_rel_end = abs((e_end - e_ini) / e_ini) if abs(e_ini) > 1e-15 else 0.0
    pd.DataFrame([{
        "t_final_anos": round(t_total_anos, 4),
        "energia_total": round(float(e_end * gm_factor), 12),
        "energia_cinetica": round(float(t_kin_end * gm_factor), 12),
        "energia_potencial": round(float((e_end - t_kin_end) * gm_factor), 12),
        "delta_e_relativo": float(f"{de_rel_end:.6e}")
    }]).to_csv(final_dir / "energia.csv", index=False)

    print("\n" + "=" * 80)
    print("✅ SIMULAÇÃO CONCLUÍDA COM SUCESSO! ESTRUTURA GERADA:")
    print("=" * 80)
    print(f"📁 Pasta de Saída Geral: {output_path}")
    print("\n   [📂 final/ -> Snapshots Consolidados no Tempo Final T=20]")
    print("      1. orbitas.csv  (elementos orbitais finais)")
    print("      2. condin.csv   (coordenadas cartesianas finais, espelho condin.xxx)")
    print("      3. encuent.csv  (encontros ordenados por menor distância)")
    print("      4. colisio.csv  (balanço de colisões)")
    print("      5. finales.csv  (balanço de ejeções)")
    print("      6. energia.csv  (balanço de conservação de energia)")
    print("      7. xyzv.csv     (snapshot cartesiano final)")
    print("\n   [📂 temporal/ -> Séries Temporais Completas T=0 até T=20]")
    print("      1. orbitas.csv  (evolução temporal orbital completa)")
    print("      2. condin.csv   (evolução cartesiana para re-arranque)")
    print("      3. encuent.csv  (cronograma de encontros)")
    print("      4. colisio.csv  (cronograma de impactos)")
    print("      5. finales.csv  (cronograma de ejeções)")
    print("      6. energia.csv  (curva temporal de energia)")
    print("      7. xyzv.csv     (histórico completo de posições e velocidades)")
    print("=" * 80)


if __name__ == "__main__":
    main()
