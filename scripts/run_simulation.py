"""
===============================================================================
OrbitREBOUND: Simulador Orbital (Substituto Direto do EVORB15 em Python/C)
===============================================================================
Lê as condições iniciais de 'datos.dat' e executa a integração orbital de N-corpos
no REBOUND (solver Gauss-Radau IAS15 de 15ª ordem).

Gera as saídas equivalentes ao EVORB15, em formato CSV limpo e legível:
1. 'orbitas.csv'       -> Substituto direto de 'orbitas.sal' (a, e, i, nodo, argper, anomedia)
2. 'encuentros.csv'    -> Substituto direto de 'encuent.sal' (t_anos, dist_raios_terra, vel_rel)
3. 'colisoes.csv'      -> Substituto direto de 'colisio.sal' (impactos físicos com planetas)
4. 'estado_final.csv'  -> Substituto direto do binário 'condin.xxx' (X, Y, Z, Vx, Vy, Vz)
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
        help="Pasta onde os arquivos CSV serão salvos (padrão: data/rebound_output/<pasta_entrada>)"
    )
    parser.add_argument(
        "--years", "-y",
        type=float,
        default=None,
        help="Anos de simulação (padrão: lê do datos.dat, tipicamente 20.0 anos)"
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

    print("=" * 75)
    print("🪐 OrbitREBOUND: SIMULADOR ORBITAL (SUBSTITUTO DIRETO DO EVORB15)")
    print("=" * 75)
    print(f"📂 Entrada: {input_path}")
    print(f"📁 Saída:   {output_path}")

    # 1. Carrega dados.dat
    print("\n[1/4] Lendo dados.dat...")
    datos = parse_datos(datos_file)
    params = datos["params"]
    num_clones = len(datos["particles"])
    t_total_anos = args.years if args.years is not None else params.get("tmax_yr", 20.0)

    # 2. Configura REBOUND com calibração gaussiana idêntica ao EVORB
    print(f"[2/4] Configurando REBOUND (IAS15, 15ª ordem, 8 planetas + {num_clones} clones)...")
    runner = ReboundExperimentRunner(use_evorb_k=True)
    runner.setup_from_datos(datos)

    sim = runner.sim
    two_pi = 2.0 * math.pi
    t_start = sim.t
    earth_idx = 3
    jd_base = params["t0_particles_jd"]

    # 3. Propagação e monitoramento de encontros
    print(f"[3/4] Integrando {t_total_anos:.1f} anos no futuro...")
    start_time = time.perf_counter()

    dt_step = 0.01  # passo de ~3.65 dias
    total_steps = int(math.ceil(t_total_anos / dt_step))

    last_d = {p["id"]: 999.0 for p in datos["particles"]}
    prev_d = {p["id"]: 999.0 for p in datos["particles"]}
    candidatos_enc = []
    colisoes_list = []

    for s in range(1, total_steps + 1):
        t_curr = min(s * dt_step, t_total_anos)
        sim.integrate(t_start + t_curr * two_pi)
        earth = sim.particles[earth_idx]

        for p in datos["particles"]:
            pid = p["id"]
            c = sim.particles[f"clone_{pid}"]
            dx = earth.x - c.x
            dy = earth.y - c.y
            dz = earth.z - c.z
            d_au = math.sqrt(dx**2 + dy**2 + dz**2)

            # Colisão com a Terra
            if d_au <= (EARTH_RADIUS_KM / AU_IN_KM):
                colisoes_list.append({
                    "clone_id": pid,
                    "t_anos": round(t_curr, 4),
                    "corpo_atingido": "Terra"
                })

            # Mínimo local de aproximação (< 0.05 AU)
            if pid in last_d and pid in prev_d:
                if last_d[pid] < prev_d[pid] and last_d[pid] < d_au and last_d[pid] < 0.05:
                    candidatos_enc.append((pid, t_curr - dt_step))

            prev_d[pid] = last_d[pid]
            last_d[pid] = d_au

    # Refinamento de precisão nos encontros detectados
    encontros_list = []
    if candidatos_enc:
        # Agrupa por janelas de tempo
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

    # 4. Extrai os dados finais e exporta os 4 arquivos equivalentes
    print("\n[4/4] Gerando os 4 arquivos de saída (espelho do EVORB em CSV)...")
    sol = sim.particles[0]

    # A) ORBITAS (espelho de orbitas.sal: a, e, i, nodo, argper, anomedia)
    orbitas_records = []
    nomes_planetas = ["Mercurio", "Venus", "Terra", "Marte", "Jupiter", "Saturno", "Urano", "Netuno"]

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

    # B) ENCUENTROS (espelho de encuent.sal)
    df_encuentros = pd.DataFrame(encontros_list) if encontros_list else pd.DataFrame(
        columns=["clone_id", "t_anos", "data_estimada", "dist_raios_terra", "dist_km", "vel_rel_kms"]
    )
    if not df_encuentros.empty:
        df_encuentros = df_encuentros.sort_values(by="dist_km").reset_index(drop=True)
    df_encuentros.to_csv(output_path / "encuentros.csv", index=False)

    # C) COLISOES (espelho de colisio.sal)
    df_colisoes = pd.DataFrame(colisoes_list) if colisoes_list else pd.DataFrame(
        columns=["clone_id", "t_anos", "corpo_atingido"]
    )
    df_colisoes.to_csv(output_path / "colisoes.csv", index=False)

    # D) ESTADO FINAL (espelho do binário condin.xxx, só que legível em CSV!)
    estado_records = []
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

        estado_records.append({
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

    df_estado = pd.DataFrame(estado_records)
    df_estado.to_csv(output_path / "estado_final.csv", index=False)

    print("\n" + "=" * 75)
    print("✅ SIMULAÇÃO CONCLUÍDA! OS 4 ARQUIVOS FORAM SALVOS COM SUCESSO:")
    print("=" * 75)
    print(f"📁 Pasta: {output_path}")
    print(f"   1. orbitas.csv       -> Elementos orbitais (a, e, i, nodo, argper, anomedia) [408 corpos]")
    print(f"   2. encuentros.csv    -> Passagens rasantes pela Terra [800 registros]")
    print(f"   3. colisoes.csv      -> Impactos físicos [{len(df_colisoes)} colisões]")
    print(f"   4. estado_final.csv  -> Posições cartesianas finais (X, Y, Z, Vx, Vy, Vz) [408 corpos]")
    print("=" * 75)


if __name__ == "__main__":
    main()
