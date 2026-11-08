import pytest
import sys
from pathlib import Path
import numpy as np

# Add src to sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.parsers.evorb_parser import (
    clean_fortran_float,
    parse_datos,
    parse_condin_binary,
    parse_encuent
)
from src.analysis.metrics import calculate_position_residual


def test_clean_fortran_float():
    assert clean_fortran_float("1.2345D-02") == 0.012345
    assert clean_fortran_float(" -5.678D+01 ") == -56.78
    assert clean_fortran_float("1.0E+00") == 1.0
    # Notação Fortran compacta sem a letra 'E' ou 'D'
    assert clean_fortran_float("1.2345-05") == 1.2345e-05
    assert clean_fortran_float("1.2345+05") == 123450.0


def test_calculate_position_residual():
    r1 = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
    r2 = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    res = calculate_position_residual(r1, r2)
    assert np.isclose(res[0], 0.0)
    assert np.isclose(res[1], 2.0)

    # Validação defensiva de dimensão incompatível
    with pytest.raises(ValueError):
        calculate_position_residual(r1, np.array([[1.0, 0.0, 0.0]]))


def test_parse_datos():
    datos_path = PROJECT_ROOT / "data" / "evorb_raw" / "datos1" / "datos.dat"
    if not datos_path.exists():
        pytest.skip("Arquivo datos.dat não encontrado para teste de integração.")

    data = parse_datos(datos_path)
    assert data["params"]["tmax_yr"] == 20.0
    assert len(data["planets"]) == 8
    assert len(data["particles"]) == 400
    # Checar Terra (planeta 3)
    earth = data["planets"][2]
    assert np.isclose(earth["a"], 0.99998836)
    assert earth["radius_km"] == 6371.0


def test_parse_condin_binary():
    condin_path = PROJECT_ROOT / "data" / "evorb_raw" / "datos1" / "condin.xxx"
    if not condin_path.exists():
        pytest.skip("Arquivo condin.xxx não encontrado.")

    data = parse_condin_binary(condin_path)
    assert data["t_final_yr"] == 20.0
    assert data["nbod"] == 408
    assert data["npla"] == 8
    assert data["positions"].shape == (408, 3)
    assert data["velocities"].shape == (408, 3)


def test_parse_encuent():
    encuent_path = PROJECT_ROOT / "data" / "evorb_raw" / "datos1" / "encuent.sal"
    if not encuent_path.exists():
        pytest.skip("Arquivo encuent.sal não encontrado.")

    encounters = parse_encuent(encuent_path)
    assert len(encounters) > 0
    # Checar que o primeiro encontro é com a Terra (Pla = 3)
    assert encounters[0]["planet_id"] == 3
    assert np.isclose(encounters[0]["dist_radii"], 130.4)
