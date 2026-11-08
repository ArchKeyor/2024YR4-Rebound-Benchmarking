import sys
from pathlib import Path
ROOT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT_DIR))

import numpy as np
import rebound
from src.parsers.evorb_parser import parse_datos

datos = parse_datos(ROOT_DIR / "data/evorb_raw/datos1/datos.dat")
planets = datos["planets"]

def ephem1_pure(a, e, inc_deg, node_deg, peri_deg, M_deg):
    inc = np.radians(inc_deg)
    anode = np.radians(node_deg)
    peri = np.radians(peri_deg)
    M = np.radians(M_deg)
    
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
    
    ck, sk = np.cos(anode), np.sin(anode)
    ci, si = np.cos(inc), np.sin(inc)
    co, so = np.cos(peri), np.sin(peri)
    
    P = np.array([ck*co - sk*ci*so, sk*co + ck*ci*so, si*so])
    Q = np.array([-ck*so - sk*ci*co, -sk*so + ck*ci*co, si*co])
    
    return x1 * P + x2 * Q

diffs = []
for p in planets:
    pos_ephem1 = ephem1_pure(p['a'], p['e'], p['inc_deg'], p['node_deg'], p['argper_deg'], p['M_deg'])
    
    sim = rebound.Simulation()
    sim.units = ('AU', 'Msun', 'yr2pi')
    sim.add(m=1.0)
    sim.add(
        m=p['mass_msun'],
        a=p['a'],
        e=p['e'],
        inc=np.radians(p['inc_deg']),
        Omega=np.radians(p['node_deg']),
        omega=np.radians(p['argper_deg']),
        M=np.radians(p['M_deg']),
        primary=sim.particles[0]
    )
    pos_reb = np.array([sim.particles[1].x, sim.particles[1].y, sim.particles[1].z])
    
    diff = np.linalg.norm(pos_ephem1 - pos_reb)
    diffs.append(diff)
    print(f"Planeta {p['id']}: Diferenca EPHEM1 vs REBOUND = {diff:.4e} AU")

print(f"\nDiferenca maxima nos 8 planetas em T=0: {max(diffs):.4e} AU")
