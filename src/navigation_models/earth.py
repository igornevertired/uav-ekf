"""Упрощённая модель Земли: радиусы кривизны и вектор g в локальной НСК [N, U, E]."""

from __future__ import annotations

import numpy as np


def earth_model(
    h: float,
    fi: float,
    lamd: float,
) -> tuple[np.ndarray, float, float, float]:
    """
    Возвращает (g_vector_NUE, R_N, R_M, unused_placeholder).

    g_vector: приближённо [0, -g, 0] в НСК North-Up-East, м/с².
    R_N, R_M: радиусы кривизны + высота, м (для интегрирования координат).
    """
    a = 6378245.0
    b = 6356856.0
    e2 = (a**2 - b**2) / a**2
    sf = np.sin(fi)
    cf = np.cos(fi)
    w = np.sqrt(1.0 - e2 * sf * sf)
    rn = a / w
    rm = a * (1.0 - e2) / (w**3)
    r_n = rn + h
    r_m = rm + h
    g0 = 9.80665 * (6371000.0 / (6371000.0 + float(h))) ** 2
    g_g = np.array([0.0, -g0, 0.0], dtype=float)
    return g_g, float(r_n), float(r_m), 0.0
