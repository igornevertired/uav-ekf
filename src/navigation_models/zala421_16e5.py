from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .bins import Bins2Params
from .noise import imu_reference_sigmas


def initial_attitude_nue(theta: float) -> tuple[np.ndarray, np.ndarray]:
    """Связанная СК: [Forward, Right, Down]. Навигация: [North, Up, East]."""
    ct, st = np.cos(theta), np.sin(theta)
    c_bn = np.array(
        [
            [ct, 0.0, st],
            [st, 0.0, -ct],
            [0.0, 1.0, 0.0],
        ],
        dtype=float,
    )
    ct2 = np.cos(theta / 2.0)
    st2 = np.sin(theta / 2.0)
    sq2 = np.sqrt(2.0)
    q = np.array([ct2 / sq2, ct2 / sq2, st2 / sq2, st2 / sq2], dtype=float)
    return c_bn, q


@dataclass(frozen=True)
class Zala421Config:
    np0: np.ndarray
    q0: np.ndarray
    c_bn0: np.ndarray
    a_true_body: np.ndarray
    w_true_body: np.ndarray
    bins_params: Bins2Params
    gnss_scheme: int
    dt: float


def zala421_default_config() -> Zala421Config:
    fi0 = np.deg2rad(55.75)
    lm0 = np.deg2rad(37.62)
    h0 = 300.0
    v_cruise = 30.0
    np0 = np.array([v_cruise, 0.0, 0.0, h0, fi0, lm0], dtype=float)

    theta0 = 0.07
    c_bn0, q0 = initial_attitude_nue(theta0)

    a_true_body = np.array([0.0, 0.0, 0.0], dtype=float)
    w_true_body = np.array([0.0, 0.0, 0.0], dtype=float)

    imu_sig = imu_reference_sigmas()
    bins_params = Bins2Params(
        dw=imu_sig["gyro_bias"].copy(),
        da=imu_sig["accel_bias"].copy(),
        kmw=np.zeros((3, 3), dtype=float),
        kma=np.zeros((3, 3), dtype=float),
        fiw=np.eye(3, dtype=float),
        fia=np.eye(3, dtype=float),
    )
    return Zala421Config(
        np0=np0,
        q0=q0,
        c_bn0=c_bn0,
        a_true_body=a_true_body,
        w_true_body=w_true_body,
        bins_params=bins_params,
        gnss_scheme=1,
        dt=0.01,
    )
