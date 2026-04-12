from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class GnssNoise:
    wn_fi: float
    wn_lm: float
    wn_vn: float
    wn_ve: float


def gnss_reference_sigmas() -> dict[str, float]:
    """Табличные 1σ параметры GPS из Table 2."""
    pseudorange = 12.0
    pseudorange_rate = 0.5
    pos_std_rad = pseudorange / 111_000.0 * np.pi / 180.0
    return {
        "pseudorange": pseudorange,
        "pseudorange_rate": pseudorange_rate,
        "phi": pos_std_rad,
        "lambda": pos_std_rad,
        "vn": pseudorange_rate,
        "ve": pseudorange_rate,
    }


def gnss_measurement(
    fi: float,
    lm: float,
    vn: float,
    ve: float,
    noise: GnssNoise,
    scheme: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    z = np.array(
        [
            fi + noise.wn_fi,
            lm + noise.wn_lm,
            vn + noise.wn_vn,
            ve + noise.wn_ve,
        ],
        dtype=float,
    )
    if scheme == 1:
        ref = gnss_reference_sigmas()
        sigma = np.array(
            [ref["phi"], ref["lambda"], ref["vn"], ref["ve"]],
            dtype=float,
        )
    elif scheme == 2:
        sigma = np.array([6.27e-7, 6.27e-7, 0.1e-7, 0.1e-7], dtype=float)
    else:
        raise ValueError("scheme должен быть 1 или 2")
    return z, sigma
