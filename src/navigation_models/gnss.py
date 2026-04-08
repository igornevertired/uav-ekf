from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class GnssNoise:
    wn_fi: float
    wn_lm: float
    wn_vn: float
    wn_ve: float


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
        sko1 = 0.5
        sko2 = 0.1
        sigma = np.array(
            [sko1 / 111_000.0 * np.pi / 180.0, sko1 / 111_000.0 * np.pi / 180.0, sko2, sko2],
            dtype=float,
        )
    elif scheme == 2:
        sigma = np.array([6.27e-7, 6.27e-7, 0.1e-7, 0.1e-7], dtype=float)
    else:
        raise ValueError("scheme должен быть 1 или 2")
    return z, sigma
