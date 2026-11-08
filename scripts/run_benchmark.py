"""
Interface de Linha de Comando (CLI) para execução do Benchmarking Orbital 2024 YR4.
Confronta a integração N-corpos no REBOUND (IAS15) contra os gabaritos do EVORB15.
"""

import argparse
import sys
import time
from pathlib import Path
import numpy as np

# Adiciona a raiz do projeto ao sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.parsers.evorb_parser import parse_datos, parse_condin_binary, parse_encuent
from src.simulation.rebound_runner import ReboundExperimentRunner
from src.analysis.metrics import calculate_position_residual, calculate_velocity_residual


def parse_args():
    parser = argparse.ArgumentParser(
        description="2024 YR4 Benchmarking CLI: EVORB15 vs. REBOUND (IAS15)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default="data/evorb_raw/datos1",
        help="Caminho relativo ou absoluto da pasta contendo os dados brutos do EVORB."
    )
    parser.add_argument(
        "--years",
        type=float,
        default=20.0,
        help="Horizonte temporal de integração contínua em anos Julianos."
    )
    parser.add_argument(
        "--use-evorb-k",
        action="store_true",
        help="Ativa calibração com a constante gaussiana k do EVORB (GM = (k*365.25)^2)."
    )
    return parser.parse_args()


def main():
    args = parse_args()
    input_dir = PROJECT_ROOT / args.input_dir if not Path(args.input_dir).is_absolute() else Path(args.input_dir)

    print("=" * 70)
    print("🚀 2024 YR4 BENCHMARKING CLI (EVORB15 vs. REBOUND IAS15)")
    print("=" * 70)

    datos_path = input_dir / "datos.dat"
    condin_path = input_dir / "condin.xxx"
    encuent_path = input_dir / "encuent.sal"

    if not datos_path.exists():
        print(f"[-] Erro: Arquivo de configuração não encontrado: {datos_path}")
        sys.exit(1)

    print(f"\n[1/4] Ingerindo dados de entrada a partir de: {input_dir.name}/")
    datos = parse_datos(datos_path)
    print(f"      -> 8 planetas e {len(datos['particles'])} clones de 2024 YR4 decodificados.")

    print("\n[2/4] Configurando simulador REBOUND (IAS15)...")
    runner = ReboundExperimentRunner(use_evorb_k=args.use_evorb_k)
    runner.setup_from_datos(datos)
    calib_str = " (Calibração Gaussiana do EVORB ativada: GM=(k*365.25)^2)" if args.use_evorb_k else " (Modo Newtoniano padrão G=1.0)"
    print(f"      -> Sistema alinhado via LLEVA e centrado no baricentro inercial.{calib_str}")

    print(f"\n[3/4] Integrando {args.years:.1f} anos na janela crítica...")
    t_start = time.perf_counter()
    E_ini = runner.energy()
    L_ini = runner.angular_momentum()
    runner.integrate_years(args.years)
    dt_wall = time.perf_counter() - t_start
    E_fim = runner.energy()
    L_fim = runner.angular_momentum()

    dE = abs((E_fim - E_ini) / (abs(E_ini) + 1e-30))
    dL = np.linalg.norm(L_fim - L_ini) / (np.linalg.norm(L_ini) + 1e-30)

    print(f"      -> Concluído em {dt_wall:.2f} segundos!")
    print(f"      -> Conservação de Energia (dE/E):            {dE:.2e}")
    print(f"      -> Conservação de Momento Angular (dL/L):   {dL:.2e}")

    # 4. Confronto com gabaritos se disponíveis
    if condin_path.exists():
        print("\n[4/4] Executando benchmarking forense contra condin.xxx...")
        condin = parse_condin_binary(condin_path)
        states_reb = runner.get_heliocentric_states()

        r_reb = np.array([[s["x"], s["y"], s["z"]] for s in states_reb])
        v_reb = np.array([[s["vx"], s["vy"], s["vz"]] for s in states_reb])
        r_evo = condin["positions"]
        v_evo = condin["velocities"]

        delta_r = calculate_position_residual(r_reb, r_evo)
        delta_v = calculate_velocity_residual(v_reb, v_evo)

        print("\n" + "-" * 70)
        print(f"🪐 RESUMO PLANETÁRIO (8 Planetas em T={args.years:.1f} anos):")
        print(f"   Delta r Médio: {np.mean(delta_r[:8]):.4e} AU ({np.mean(delta_r[:8]) * 149597870.7:.1f} km)")
        print(f"   Delta v Médio: {np.mean(delta_v[:8]):.4e} AU/ano")

        print(f"\n☄️  RESUMO DOS CLONES ({len(delta_r[8:])} partículas de teste):")
        print(f"   Delta r Médio: {np.mean(delta_r[8:]):.4e} AU ({np.mean(delta_r[8:]) * 149597870.7 / 1e6:.1f} milhões de km)")
        print(f"   Delta v Médio: {np.mean(delta_v[8:]):.4e} AU/ano")
        print("-" * 70)
    else:
        print(f"\n[4/4] Arquivo {condin_path.name} não encontrado para confronto numérico.")

    print("\n✅ Execução finalizada com sucesso.")


if __name__ == "__main__":
    main()
