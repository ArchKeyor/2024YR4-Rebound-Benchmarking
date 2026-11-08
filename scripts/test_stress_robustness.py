import sys
from pathlib import Path
ROOT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT_DIR))

import numpy as np
from src.parsers.evorb_parser import parse_datos
from src.simulation.rebound_runner import ReboundExperimentRunner

def run_stress():
    print("=" * 65)
    print("🔬 BATERIA DE ESTRESSE E ROBUSTEZ TÉCNICA")
    print("=" * 65)

    datos = parse_datos(ROOT_DIR / "data/evorb_raw/datos1/datos.dat")

    print("\n[Teste 1] Idempotência do Setup (3 reinicializações consecutivas):")
    runner = ReboundExperimentRunner()
    for k in range(3):
        runner.setup_from_datos(datos)
        print(f"  -> Ciclo {k+1}: Total Corpos = {runner.sim.N}, Ativos = {runner.sim.N_active}, Alinhado = {runner.alinhado}")
        assert runner.sim.N == 409
        assert runner.sim.N_active == 9
    print("  ✅ APROVADO: Zero vazamento de corpos ou acumulação de memória.")

    print("\n[Teste 2] Integração por Fatias Incrementais (5 anos x 4 = 20 anos):")
    E0 = runner.energy()
    L0 = runner.angular_momentum()
    for step in range(4):
        runner.integrate_years(5.0)
        E_curr = runner.energy()
        L_curr = runner.angular_momentum()
        dE = abs((E_curr - E0) / (abs(E0) + 1e-30))
        dL = np.linalg.norm(L_curr - L0) / (np.linalg.norm(L0) + 1e-30)
        t_yr = runner.sim.t / (2.0 * np.pi)
        print(f"  -> Ano {int((step+1)*5):02d}: dE/E = {dE:.4e} | dL/L = {dL:.4e} | t_acumulado = {t_yr:.2f} anos")
        assert dE < 1e-12
        assert dL < 1e-12
    print("  ✅ APROVADO: Estabilidade mantida em integrações particionadas.")

    print("\n[Teste 3] Varredura de Integridade Numérica (NaNs / Infs nos 408 Corpos):")
    states = runner.get_heliocentric_states()
    assert len(states) == 408
    nan_count = 0
    inf_count = 0
    for s in states:
        vals = [s["x"], s["y"], s["z"], s["vx"], s["vy"], s["vz"]]
        if any(np.isnan(v) for v in vals):
            nan_count += 1
        if any(np.isinf(v) for v in vals):
            inf_count += 1
    print(f"  -> NaNs encontrados: {nan_count}")
    print(f"  -> Infs encontrados: {inf_count}")
    assert nan_count == 0 and inf_count == 0
    print("  ✅ APROVADO: 100% dos corpos com coordenadas e velocidades finitas.")

    print("\n" + "=" * 65)
    print("🎉 DIAGNÓSTICO: O SOFTWARE MANTEVE E ELEVOU SUA ROBUSTEZ!")
    print("=" * 65)

if __name__ == "__main__":
    run_stress()

