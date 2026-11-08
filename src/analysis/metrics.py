"""Mathematical and statistical metrics for EVORB vs REBOUND benchmarking."""

import numpy as np
from typing import Dict, Any, Tuple


def calculate_position_residual(r_rebound: np.ndarray, r_evorb: np.ndarray) -> np.ndarray:
    """
    Calculates absolute Euclidean distance residual:
    Delta r = || r_rebound - r_evorb ||
    """
    r_reb = np.asarray(r_rebound, dtype=np.float64)
    r_evo = np.asarray(r_evorb, dtype=np.float64)
    if r_reb.shape != r_evo.shape:
        raise ValueError(f"Dimensões incompatíveis para resíduo de posição: {r_reb.shape} vs {r_evo.shape}")
    diff = r_reb - r_evo
    return np.linalg.norm(diff, axis=-1)


def calculate_velocity_residual(v_rebound: np.ndarray, v_evorb: np.ndarray) -> np.ndarray:
    """
    Calculates absolute velocity residual:
    Delta v = || v_rebound - v_evorb ||
    """
    v_reb = np.asarray(v_rebound, dtype=np.float64)
    v_evo = np.asarray(v_evorb, dtype=np.float64)
    if v_reb.shape != v_evo.shape:
        raise ValueError(f"Dimensões incompatíveis para resíduo de velocidade: {v_reb.shape} vs {v_evo.shape}")
    diff = v_reb - v_evo
    return np.linalg.norm(diff, axis=-1)


def find_closest_approach(times: np.ndarray, r_target: np.ndarray, r_earth: np.ndarray) -> Tuple[float, float]:
    """
    Determines minimum approach distance and corresponding timestamp:
    Returns (t_min, min_distance).
    """
    rel_dist = np.linalg.norm(r_target - r_earth, axis=-1)
    min_idx = np.argmin(rel_dist)
    return float(times[min_idx]), float(rel_dist[min_idx])
