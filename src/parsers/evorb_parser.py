"""Parser module for legacy EVORB15 Fortran files."""

import re
import struct
from pathlib import Path
from typing import Dict, Any, Generator, List, Tuple
import numpy as np


def clean_fortran_float(val_str: str) -> float:
    """
    Converts Fortran float string to Python float.
    Handles 'D' / 'd' double-precision exponents and tight-field notations
    where the 'E'/'D' letter was omitted (e.g. '1.2345-05' -> '1.2345E-05').
    """
    clean = val_str.strip().upper().replace("D", "E")
    if "E" not in clean:
        clean = re.sub(r"(?<=\d)(?=[+-]\d)", "E", clean)
    return float(clean)


def parse_datos(file_path: Path | str) -> Dict[str, Any]:
    """
    Parses full simulation setup from datos.dat:
    - Integrator options (relativistic corrections, step size, limits)
    - 8 major planets orbital elements and masses
    - 400 asteroid test particles orbital elements (2024 YR4 clones)
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Arquivo datos.dat não encontrado em: {path}")

    with open(path, "r", encoding="latin1") as f:
        lines = [line.strip() for line in f if line.strip()]

    # Filter out pure comment lines starting with '*'
    data_lines = [l for l in lines if not l.startswith("*")]

    iefrel = int(data_lines[0])
    ixyz = int(data_lines[1])
    rcrit1, rcrit2 = [clean_fortran_float(x) for x in data_lines[2].split()]
    cmass, rstar = [clean_fortran_float(x) for x in data_lines[3].split()]
    tmax = clean_fortran_float(data_lines[4])
    tsal = clean_fortran_float(data_lines[5])
    h_step = clean_fortran_float(data_lines[6])
    istat = int(data_lines[7])
    npla = int(data_lines[8])
    t0_planets = clean_fortran_float(data_lines[9])

    # Parse Planets: L, a, e, i, Omega, omega, M, mass, radius
    planets = []
    idx = 10
    for _ in range(npla):
        tokens = data_lines[idx].split()
        planets.append({
            "id": int(tokens[0]),
            "a": clean_fortran_float(tokens[1]),
            "e": clean_fortran_float(tokens[2]),
            "inc_deg": clean_fortran_float(tokens[3]),
            "node_deg": clean_fortran_float(tokens[4]),
            "argper_deg": clean_fortran_float(tokens[5]),
            "M_deg": clean_fortran_float(tokens[6]),
            "mass_msun": clean_fortran_float(tokens[7]),
            "radius_km": clean_fortran_float(tokens[8])
        })
        idx += 1

    nparti = int(data_lines[idx])
    idx += 1
    t0_particles = clean_fortran_float(data_lines[idx])
    idx += 1

    # Parse Particles: L, a, e, i, Omega, omega, M
    particles = []
    for _ in range(nparti):
        tokens = data_lines[idx].split()
        particles.append({
            "id": int(tokens[0]),
            "a": clean_fortran_float(tokens[1]),
            "e": clean_fortran_float(tokens[2]),
            "inc_deg": clean_fortran_float(tokens[3]),
            "node_deg": clean_fortran_float(tokens[4]),
            "argper_deg": clean_fortran_float(tokens[5]),
            "M_deg": clean_fortran_float(tokens[6]),
        })
        idx += 1

    return {
        "params": {
            "iefrel": iefrel,
            "ixyz": ixyz,
            "rcrit1_au": rcrit1,
            "rcrit2_au": rcrit2,
            "cmass_msun": cmass,
            "rstar_km": rstar,
            "tmax_yr": tmax,
            "tsal_yr": tsal,
            "h_step_yr": h_step,
            "istat": istat,
            "t0_planets_jd": t0_planets,
            "t0_particles_jd": t0_particles,
        },
        "planets": planets,
        "particles": particles
    }


def parse_condin_binary(file_path: Path | str) -> Dict[str, Any]:
    """
    Parses Fortran unformatted binary checkpoint file (condin.xxx).
    Unpacks positions, velocities, and physics metadata for all bodies.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Arquivo condin.xxx não encontrado em: {path}")

    with open(path, "rb") as f:
        # Header: T (real*8), NBOD (int*4), NPLA (int*4)
        m_start = struct.unpack("<i", f.read(4))[0]
        if m_start != 16:
            raise ValueError(f"Marcador inicial de cabeçalho Fortran inválido: {m_start} (esperado 16)")
        t = struct.unpack("<d", f.read(8))[0]
        nbod, npla = struct.unpack("<ii", f.read(8))
        m_end = struct.unpack("<i", f.read(4))[0]
        if m_end != 16:
            raise ValueError(f"Marcador final de cabeçalho Fortran inválido: {m_end} (esperado 16)")

        positions = np.zeros((nbod, 3), dtype=np.float64)
        velocities = np.zeros((nbod, 3), dtype=np.float64)
        masses = np.zeros(nbod, dtype=np.float64)
        radii_au = np.zeros(nbod, dtype=np.float64)
        names = []

        for i in range(nbod):
            # Pos (X, Y, Z) - 24 bytes
            p1 = struct.unpack("<i", f.read(4))[0]
            positions[i] = struct.unpack("<ddd", f.read(24))
            p2 = struct.unpack("<i", f.read(4))[0]
            if p1 != 24 or p2 != 24:
                raise ValueError(f"Marcadores de posição do corpo {i} inválidos: {p1}, {p2}")

            # Vel (Vx, Vy, Vz) - 24 bytes
            v1 = struct.unpack("<i", f.read(4))[0]
            velocities[i] = struct.unpack("<ddd", f.read(24))
            v2 = struct.unpack("<i", f.read(4))[0]
            if v1 != 24 or v2 != 24:
                raise ValueError(f"Marcadores de velocidade do corpo {i} inválidos: {v1}, {v2}")

            # Metadata: IPRIM, RPH, MPLA, RPLA, NOMBRE - 32 bytes
            m1 = struct.unpack("<i", f.read(4))[0]
            iprim = struct.unpack("<i", f.read(4))[0]
            rph = struct.unpack("<d", f.read(8))[0]
            mpla = struct.unpack("<d", f.read(8))[0]
            rpla = struct.unpack("<d", f.read(8))[0]
            nombre = struct.unpack("<i", f.read(4))[0]
            m2 = struct.unpack("<i", f.read(4))[0]
            if m1 != 32 or m2 != 32:
                raise ValueError(f"Marcadores de metadados do corpo {i} inválidos: {m1}, {m2}")

            radii_au[i] = rph
            masses[i] = mpla
            names.append(nombre)

        # Garantia de consumo de 100% dos bytes
        sobra = f.read()
        if len(sobra) > 0:
            raise ValueError(f"Bytes órfãos detectados no final do condin.xxx: {len(sobra)} bytes não lidos")

    return {
        "t_final_yr": t,
        "nbod": nbod,
        "npla": npla,
        "positions": positions,
        "velocities": velocities,
        "masses": masses,
        "radii_au": radii_au,
        "names": names
    }


def parse_encuent(file_path: Path | str) -> List[Dict[str, Any]]:
    """
    Parses encounter log file (encuent.sal).
    Columns: t enc, Dm/R, Vp(K/s), energ, Pla, Par.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Arquivo encuent.sal não encontrado em: {path}")

    encounters = []
    with open(path, "r", encoding="latin1") as f:
        lines = f.readlines()

    # Skip header
    for line in lines[1:]:
        tokens = line.strip().split()
        if len(tokens) >= 6:
            encounters.append({
                "t_enc_yr": float(tokens[0]),
                "dist_radii": float(tokens[1]),
                "v_rel_kms": float(tokens[2]),
                "energy": float(tokens[3]),
                "planet_id": int(tokens[4]),
                "particle_id": int(tokens[5])
            })

    return encounters


def parse_sal_stream(file_path: Path | str, chunk_size: int = 10000) -> Generator[np.ndarray, None, None]:
    """Memory-efficient streaming generator to read large .sal files."""
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Arquivo .sal não encontrado em: {path}")

    batch = []
    with open(path, "r", encoding="latin1") as f:
        for line in f:
            tokens = line.strip().split()
            if not tokens or not tokens[0].replace(".", "", 1).replace("-", "", 1).isdigit():
                continue
            try:
                row = [clean_fortran_float(t) for t in tokens]
                batch.append(row)
                if len(batch) >= chunk_size:
                    yield np.array(batch, dtype=np.float64)
                    batch = []
            except ValueError:
                continue

    if batch:
        yield np.array(batch, dtype=np.float64)
