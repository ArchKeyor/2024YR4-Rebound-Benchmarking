"""
Testes automatizados para o módulo de simulação com REBOUND.
"""

import pytest
from pathlib import Path
import numpy as np
from src.parsers.evorb_parser import parse_datos
from src.simulation.rebound_runner import ReboundExperimentRunner

DATOS1_PATH = Path("C:/Users/Kauê/Documents/02_Programacao_e_Tech/Projetos/2024YR4-Rebound-Benchmarking/data/evorb_raw/datos1/datos.dat")


def test_rebound_runner_setup_and_short_integration():
    """Valida o setup completo e a integração de curto prazo (0.1 ano)."""
    assert DATOS1_PATH.exists(), f"Arquivo não encontrado: {DATOS1_PATH}"
    datos = parse_datos(str(DATOS1_PATH))

    runner = ReboundExperimentRunner()
    runner.setup_from_datos(datos)

    # Verifica se 409 partículas foram carregadas (Sol + 8 planetas + 400 clones)
    assert runner.sim.N == 409
    assert runner.sim.N_active == 9
    assert runner.alinhado is True

    # Conservação de energia após 0.1 ano
    E_ini = runner.energy()
    L_ini = runner.angular_momentum()
    runner.integrate_years(0.1)
    E_fim = runner.energy()
    L_fim = runner.angular_momentum()

    dE = abs((E_fim - E_ini) / E_ini)
    assert dE < 1e-12, f"Erro de conservação de energia muito alto: {dE}"

    dL = np.linalg.norm(L_fim - L_ini) / np.linalg.norm(L_ini)
    assert dL < 1e-12, f"Erro de conservação de momento angular muito alto: {dL}"

    # Checa estado da Terra
    terra_state = runner.get_body_state("Terra")
    r_terra = (terra_state["x"]**2 + terra_state["y"]**2 + terra_state["z"]**2)**0.5
    assert 0.98 < r_terra < 1.02, f"Raio orbital da Terra fora do esperado: {r_terra}"


def test_rebound_runner_heliocentric_states():
    """Valida a conversão heliocêntrica e as unidades de velocidade (AU/ano)."""
    datos = parse_datos(str(DATOS1_PATH))
    runner = ReboundExperimentRunner()
    runner.setup_from_datos(datos)

    states = runner.get_heliocentric_states()
    assert len(states) == 408  # 8 planetas + 400 clones

    # O terceiro corpo (índice 2) é a Terra
    terra = states[2]
    assert terra["name"] == "Terra"
    assert states[8]["name"] == "clone_9"

    r_terra = (terra["x"]**2 + terra["y"]**2 + terra["z"]**2)**0.5
    v_terra = (terra["vx"]**2 + terra["vy"]**2 + terra["vz"]**2)**0.5

    assert 0.98 < r_terra < 1.02, f"Raio heliocêntrico da Terra inválido: {r_terra}"
    # Velocidade orbital da Terra em AU/ano deve ser ~2*pi (aprox. 6.28 AU/ano)
    assert 6.0 < v_terra < 6.5, f"Velocidade heliocêntrica da Terra em AU/ano inválida: {v_terra}"


def test_rebound_runner_defensive_and_idempotency():
    """Valida proteção contra chamadas fora de ordem e chamadas repetidas (idempotência)."""
    runner = ReboundExperimentRunner()

    # 1. Chamar integrate sem setup deve falhar com RuntimeError
    with pytest.raises(RuntimeError):
        runner.integrate_years(1.0)

    # 2. Chamar estados sem setup deve falhar com RuntimeError
    with pytest.raises(RuntimeError):
        runner.get_heliocentric_states()

    # 3. Setup e Idempotência: rodar setup duas vezes NÃO deve duplicar corpos
    datos = parse_datos(str(DATOS1_PATH))
    runner.setup_from_datos(datos)
    assert runner.sim.N == 409

    # Segundo setup consecutivo na mesma instância
    runner.setup_from_datos(datos)
    assert runner.sim.N == 409  # Continua 409, não 818!

    # 4. Tempo negativo de integração deve falhar com ValueError
    with pytest.raises(ValueError):
        runner.integrate_years(-1.0)


def test_ephem1_analytical_equivalence_with_rebound():
    """
    Comprova que a matriz perifocal de rotação P e Q de Murray & Dermott
    usada na sub-rotina EPHEM1 do Fortran é matematicamente idêntica
    à conversão interna de elementos orbitais do REBOUND (tolerância < 1e-13 AU).
    """
    import rebound
    import numpy as np

    datos = parse_datos(str(DATOS1_PATH))
    planets = datos["planets"]

    for p in planets:
        inc = np.radians(p["inc_deg"])
        anode = np.radians(p["node_deg"])
        peri = np.radians(p["argper_deg"])
        M = np.radians(p["M_deg"])
        e = p["e"]
        a = p["a"]

        # Resolução de Kepler pura
        E = M
        for _ in range(100):
            f = E - e * np.sin(E) - M
            f_prime = 1.0 - e * np.cos(E)
            dE = -f / f_prime
            E += dE
            if abs(dE) < 1e-15:
                break

        x1 = a * (np.cos(E) - e)
        x2 = a * np.sqrt(1.0 - e**2) * np.sin(E)

        # Vetores P e Q da sub-rotina EPHEM1 (linhas 1833-1839 do evorb15corr.for)
        ck, sk = np.cos(anode), np.sin(anode)
        ci, si = np.cos(inc), np.sin(inc)
        co, so = np.cos(peri), np.sin(peri)

        P = np.array([ck * co - sk * ci * so, sk * co + ck * ci * so, si * so])
        Q = np.array([-ck * so - sk * ci * co, -sk * so + ck * ci * co, si * co])
        pos_ephem1 = x1 * P + x2 * Q

        sim = rebound.Simulation()
        sim.units = ("AU", "Msun", "yr2pi")
        sim.add(m=1.0)
        sim.add(
            m=p["mass_msun"],
            a=a,
            e=e,
            inc=inc,
            Omega=anode,
            omega=peri,
            M=M,
            primary=sim.particles[0]
        )
        pos_reb = np.array([sim.particles[1].x, sim.particles[1].y, sim.particles[1].z])

        diff = np.linalg.norm(pos_ephem1 - pos_reb)
        assert diff < 1e-13, f"Discrepância analítica em T=0 no planeta {p['id']}: {diff} AU"


def test_rebound_runner_gaussian_calibration():
    """Valida a calibração com a constante gaussiana k do EVORB (dr da Terra cai para < 1e-4 AU)."""
    from src.parsers.evorb_parser import parse_condin_binary
    from src.analysis.metrics import calculate_position_residual

    condin_path = Path("C:/Users/Kauê/Documents/02_Programacao_e_Tech/Projetos/2024YR4-Rebound-Benchmarking/data/evorb_raw/datos1/condin.xxx")
    if not condin_path.exists():
        pytest.skip("condin.xxx não encontrado.")

    datos = parse_datos(str(DATOS1_PATH))
    condin = parse_condin_binary(str(condin_path))

    runner = ReboundExperimentRunner(use_evorb_k=True)
    runner.setup_from_datos(datos)
    runner.integrate_years(20.0)

    states = runner.get_heliocentric_states()
    r_reb = np.array([[s["x"], s["y"], s["z"]] for s in states[:8]])
    r_evo = condin["positions"][:8]

    dr = calculate_position_residual(r_reb, r_evo)
    # Terra (índice 2) deve ter dr < 1e-4 AU (~8.500 km em vez de 813.000 km)
    dr_terra = dr[2]
    assert dr_terra < 1e-4, f"Resíduo da Terra após calibração de G muito alto: {dr_terra} AU"
    assert dr_terra * 149597870.7 < 10000.0, "Resíduo da Terra deve ser inferior a 10.000 km"
