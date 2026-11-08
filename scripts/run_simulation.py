"""
===============================================================================
OrbitREBOUND: Simulador Orbital (Substituto Direto do EVORB15 em Python/C)
===============================================================================
Lê as condições iniciais de 'datos.dat' e executa a integração orbital de N-corpos
no REBOUND com o solver Gauss-Radau IAS15 de 15ª ordem.

Gera exatamente os 7 arquivos de saída equivalentes aos gerados pelo EVORB15,
todos em formato CSV limpo, padronizado e diretamente analisável:

1. 'orbitas.csv'   -> Substituto de 'orbitas.sal' (a, e, i, nodo, argper, anomedia)
2. 'condin.csv'    -> Substituto de 'condin.xxx'  (estado cartesiano final: X, Y, Z, Vx, Vy, Vz)
3. 'encuent.csv'   -> Substituto de 'encuent.sal' (detalhes de encontros próximos / rasantes)
4. 'colisio.csv'   -> Substituto de 'colisio.sal' (registro de colisões físicas com planetas)
5. 'finales.csv'   -> Substituto de 'finales.sal' (corpos ejetados ou eliminados por limites)
6. 'energia.csv'   -> Substituto de 'energia.sal' (auditoria de conservação de energia do sistema)
7. 'xyzv.csv'      -> Substituto de 'xyzv.sal'    (histórico temporal cartesiano X, Y, Z, Vx, Vy, Vz)
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
        help="Pasta onde os 7 arquivos CSV serão salvos (padrão: data/rebound_output/<pasta_entrada>)"
    )
    parser.add_argument(
        "--years", "-y",
        type=float,
        default=None,
        help="Anos de simulação (padrão: lê do datos.dat, tipicamente 20.0 anos)"
    )
    parser.add_argument(
        "--dt-xyzv",
        type=float,
        default=0.1,
        help="Intervalo de amostragem em anos para o arquivo xyzv.csv (padrão: 0.1 ano = ~36.5 dias)"
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
    output_path.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("🪐 OrbitREBOUND: SIMULADOR ORBITAL (SUBSTITUTO DIRETO DO EVORB15)")
    print("=" * 80)
    print(f"📂 Entrada: {input_path}")
    print(f"📁 Saída:   {output_path}")

    # 1. Carrega dados.dat
    print("\n[1/4] Lendo dados.dat...")
    datos = parse_datos(datos_file)
    params = datos["params"]
    num_clones = len(datos["particles"])
    t_total_anos = args.years if args.years is not None else params.get("tmax_yr", 20.0)
    r_crit_min = params.get("rcrit1_au", 0.1)
    r_crit_max = params.get("rcrit2_au", 100.0)

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

    # 3. Propagação e monitoramento de encontros, órbitas, energia e posições
    print(f"[3/4] Integrando {t_total_anos:.1f} anos no futuro com amostragem dos 7 arquivos...")
    start_time = time.perf_counter()

    dt_step = 0.01  # passo de detecção fina de ~3.65 dias
    dt_xyzv = max(0.001, args.dt_xyzv)
    total_steps = int(math.ceil(t_total_anos / dt_step))

    last_d = {p["id"]: 999.0 for p in datos["particles"]}
    prev_d = {p["id"]: 999.0 for p in datos["particles"]}
    candidatos_enc = []
    colisoes_list = []
    finales_list = []
    energia_records = []
    xyzv_records = []

    # Energia inicial para cálculo relativo delta E / E
    e_ini = sim.energy()
    next_xyzv_t = 0.0

    # Gravação do estado inicial (t = 0) em energia e xyzv
    sol0 = sim.particles[0]
    # Energia cinética (Sol + 8 planetas)
    t_kin0 = 0.5 * sum(sim.particles[i].m * (sim.particles[i].vx**2 + sim.particles[i].vy**2 + sim.particles[i].vz**2) for i in range(9))
    energia_records.append({
        "t_anos": 0.0,
        "energia_total": round(float(e_ini), 12),
        "energia_cinetica": round(float(t_kin0), 12),
        "energia_potencial": round(float(e_ini - t_kin0), 12),
        "delta_e_relativo": 0.0
    })

    # Snapshot inicial em xyzv
    for i in range(1, sim.N):
        p = sim.particles[i]
        nome = nomes_planetas[i - 1] if i < 9 else f"clone_{p.name.split('_')[-1]}"
        tipo = "planeta" if i < 9 else "clone"
        xyzv_records.append({
            "t_anos": 0.0,
            "N": i,
            "nome": nome,
            "tipo": tipo,
            "x_au": round(float(p.x - sol0.x), 6),
            "y_au": round(float(p.y - sol0.y), 6),
            "z_au": round(float(p.z - sol0.z), 6),
            "vx_au_ano": round(float((p.vx - sol0.vx) * two_pi), 6),
            "vy_au_ano": round(float((p.vy - sol0.vy) * two_pi), 6),
            "vz_au_ano": round(float((p.vz - sol0.vz) * two_pi), 6)
        })

    next_xyzv_t += dt_xyzv

    for s in range(1, total_steps + 1):
        t_curr = min(s * dt_step, t_total_anos)
        sim.integrate(t_start + t_curr * two_pi)
        sol = sim.particles[0]
        earth = sim.particles[earth_idx]

        # Amostragem periódica de energia e xyzv
        if t_curr >= (next_xyzv_t - 1e-6) or s == total_steps:
            e_now = sim.energy()
            t_kin = 0.5 * sum(sim.particles[i].m * (sim.particles[i].vx**2 + sim.particles[i].vy**2 + sim.particles[i].vz**2) for i in range(9))
            de_rel = abs((e_now - e_ini) / e_ini) if abs(e_ini) > 1e-15 else 0.0

            energia_records.append({
                "t_anos": round(t_curr, 4),
                "energia_total": round(float(e_now), 12),
                "energia_cinetica": round(float(t_kin), 12),
                "energia_potencial": round(float(e_now - t_kin), 12),
                "delta_e_relativo": float(f"{de_rel:.6e}")
            })

            for i in range(1, sim.N):
                p = sim.particles[i]
                nome = nomes_planetas[i - 1] if i < 9 else f"clone_{p.name.split('_')[-1]}"
                tipo = "planeta" if i < 9 else "clone"
                xyzv_records.append({
                    "t_anos": round(t_curr, 4),
                    "N": i,
                    "nome": nome,
                    "tipo": tipo,
                    "x_au": round(float(p.x - sol.x), 6),
                    "y_au": round(float(p.y - sol.y), 6),
                    "z_au": round(float(p.z - sol.z), 6),
                    "vx_au_ano": round(float((p.vx - sol.vx) * two_pi), 6),
                    "vy_au_ano": round(float((p.vy - sol.vy) * two_pi), 6),
                    "vz_au_ano": round(float((p.vz - sol.vz) * two_pi), 6)
                })

            next_xyzv_t += dt_xyzv

        # Monitoramento de encontros, limites e colisões para cada clone
        for p in datos["particles"]:
            pid = p["id"]
            c = sim.particles[f"clone_{pid}"]
            dx = earth.x - c.x
            dy = earth.y - c.y
            dz = earth.z - c.z
            d_au = math.sqrt(dx**2 + dy**2 + dz**2)

            # Colisão física com a Terra (raio < raio terrestre)
            if d_au <= (EARTH_RADIUS_KM / AU_IN_KM):
                colisoes_list.append({
                    "clone_id": pid,
                    "t_anos": round(t_curr, 4),
                    "corpo_atingido": "Terra",
                    "dist_km": round(d_au * AU_IN_KM, 1),
                    "vel_impacto_kms": round(math.sqrt((earth.vx - c.vx)**2 + (earth.vy - c.vy)**2 + (earth.vz - c.vz)**2) * two_pi * AU_IN_KM / (JULIAN_YEAR_DAYS * 86400.0), 2)
                })

            # Checagem de limites de ejeção / colisão solar (finales.sal)
            r_sol = math.sqrt((c.x - sol.x)**2 + (c.y - sol.y)**2 + (c.z - sol.z)**2)
            if r_sol < r_crit_min:
                finales_list.append({
                    "clone_id": pid,
                    "t_anos": round(t_curr, 4),
                    "motivo_eliminacao": "queda_solar",
                    "dist_au": round(r_sol, 4)
                })
            elif r_sol > r_crit_max:
                finales_list.append({
                    "clone_id": pid,
                    "t_anos": round(t_curr, 4),
                    "motivo_eliminacao": "ejecao_sistema",
                    "dist_au": round(r_sol, 4)
                })

            # Detecção de mínimo local de aproximação (< 0.05 AU)
            if pid in last_d and pid in prev_d:
                if last_d[pid] < prev_d[pid] and last_d[pid] < d_au and last_d[pid] < 0.05:
                    candidatos_enc.append((pid, t_curr - dt_step))

            prev_d[pid] = last_d[pid]
            last_d[pid] = d_au

    # Refinamento de precisão micrométrica nos encontros detectados
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
                    "clone_id": pid,
                    "t_anos": round(best_t[pid], 4),
                    "data_estimada": jd_to_gregorian(jd_enc),
                    "dist_raios_terra": round(d_radii, 1),
                    "dist_km": round(d_km, 1),
                    "vel_rel_kms": round(best_v[pid], 2)
                })

    elapsed = time.perf_counter() - start_time
    print(f"      -> Concluído em {elapsed:.2f} s ({t_total_anos/elapsed:.1f} anos/s)!")

    # 4. Geração dos 7 Arquivos Oficiais Padronizados em CSV
    print("\n[4/4] Gerando os 7 arquivos de saída oficiais em CSV...")
    sol = sim.particles[0]

    # 1. ORBITAS (substituto de orbitas.sal: a, e, i, nodo, argper, anomedia)
    orbitas_records = []
    for i in range(1, sim.N):
        p = sim.particles[i]
        orb = p.orbit(primary=sol)
        nome = nomes_planetas[i - 1] if i < 9 else f"clone_{p.name.split('_')[-1]}"
        tipo = "planeta" if i < 9 else "clone"

        orbitas_records.append({
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
    df_orbitas = pd.DataFrame(orbitas_records)
    df_orbitas.to_csv(output_path / "orbitas.csv", index=False)

    # 2. CONDIN (substituto direto do condin.xxx binário: X, Y, Z, Vx, Vy, Vz no instante final)
    condin_records = []
    for i in range(1, sim.N):
        p = sim.particles[i]
        nome = nomes_planetas[i - 1] if i < 9 else f"clone_{p.name.split('_')[-1]}"
        tipo = "planeta" if i < 9 else "clone"
        x_h = p.x - sol.x
        y_h = p.y - sol.y
        z_h = p.z - sol.z
        vx_h = (p.vx - sol.vx) * two_pi
        vy_h = (p.vy - sol.vy) * two_pi
        vz_h = (p.vz - sol.vz) * two_pi

        condin_records.append({
            "N": i,
            "nome": nome,
            "tipo": tipo,
            "x_au": round(float(x_h), 6),
            "y_au": round(float(y_h), 6),
            "z_au": round(float(z_h), 6),
            "vx_au_ano": round(float(vx_h), 6),
            "vy_au_ano": round(float(vy_h), 6),
            "vz_au_ano": round(float(vz_h), 6)
        })
    df_condin = pd.DataFrame(condin_records)
    df_condin.to_csv(output_path / "condin.csv", index=False)

    # 3. ENCUENT (substituto de encuent.sal)
    df_encuent = pd.DataFrame(encontros_list) if encontros_list else pd.DataFrame(
        columns=["clone_id", "t_anos", "data_estimada", "dist_raios_terra", "dist_km", "vel_rel_kms"]
    )
    if not df_encuent.empty:
        df_encuent = df_encuent.sort_values(by="dist_km").reset_index(drop=True)
    df_encuent.to_csv(output_path / "encuent.csv", index=False)

    # 4. COLISIO (substituto de colisio.sal)
    df_colisio = pd.DataFrame(colisoes_list) if colisoes_list else pd.DataFrame(
        columns=["clone_id", "t_anos", "corpo_atingido", "dist_km", "vel_impacto_kms"]
    )
    df_colisio.to_csv(output_path / "colisio.csv", index=False)

    # 5. FINALES (substituto de finales.sal)
    df_finales = pd.DataFrame(finales_list) if finales_list else pd.DataFrame(
        columns=["clone_id", "t_anos", "motivo_eliminacao", "dist_au"]
    )
    df_finales.to_csv(output_path / "finales.csv", index=False)

    # 6. ENERGIA (substituto de energia.sal)
    df_energia = pd.DataFrame(energia_records)
    df_energia.to_csv(output_path / "energia.csv", index=False)

    # 7. XYZV (substituto de xyzv.sal)
    df_xyzv = pd.DataFrame(xyzv_records)
    df_xyzv.to_csv(output_path / "xyzv.csv", index=False)

    print("\n" + "=" * 80)
    print("✅ SIMULAÇÃO CONCLUÍDA COM SUCESSO! OS 7 ARQUIVOS PADRONIZADOS:")
    print("=" * 80)
    print(f"📁 Pasta de Saída: {output_path}")
    print(f"   1. orbitas.csv  -> Elementos orbitais finais (a, e, i, nodo, argper, anomed) [{len(df_orbitas)} corpos]")
    print(f"   2. condin.csv   -> Coordenadas cartesianas finais (X, Y, Z, Vx, Vy, Vz)    [{len(df_condin)} corpos]")
    print(f"   3. encuent.csv  -> Detalhes de encontros próximos / rasantes               [{len(df_encuent)} registros]")
    print(f"   4. colisio.csv  -> Registro de impactos na superfície de planetas          [{len(df_colisio)} impactos]")
    print(f"   5. finales.csv  -> Registro de partículas eliminadas / ejetadas            [{len(df_finales)} ejetados]")
    print(f"   6. energia.csv  -> Auditoria temporal de conservação de energia            [{len(df_energia)} amostras]")
    print(f"   7. xyzv.csv     -> Histórico temporal das trajetórias cartesianas           [{len(df_xyzv)} pontos]")
    print("=" * 80)


if __name__ == "__main__":
    main()
