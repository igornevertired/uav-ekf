from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .earth import earth_model
from .noise import ImuNoise


@dataclass(frozen=True)
class Bins2Params:
    dw: np.ndarray
    da: np.ndarray
    kmw: np.ndarray
    kma: np.ndarray
    fiw: np.ndarray
    fia: np.ndarray


def _dus(w_true: np.ndarray, params: Bins2Params, noise: ImuNoise) -> np.ndarray:
    w_p = params.fiw @ w_true
    return (np.eye(3) + params.kmw) @ w_p + params.dw + noise.wn_w


def _dlu(a_true: np.ndarray, params: Bins2Params, noise: ImuNoise) -> np.ndarray:
    a_p = params.fia @ a_true
    return (np.eye(3) + params.kma) @ a_p + params.da + noise.wn_a


def gyro_measurement(w_true: np.ndarray, params: Bins2Params, noise: ImuNoise) -> np.ndarray:
    return _dus(w_true=w_true, params=params, noise=noise)


def accel_measurement(a_true: np.ndarray, params: Bins2Params, noise: ImuNoise) -> np.ndarray:
    return _dlu(a_true=a_true, params=params, noise=noise)


def _navigation(
    a_dlu: np.ndarray,
    c_bn: np.ndarray,
    np_state: np.ndarray,
    dt: float,
) -> tuple[np.ndarray, np.ndarray]:
    ue = 7292115.0e-11
    vbe_n = np_state[:3].copy()
    h = np_state[3]
    fi = np_state[4]
    lamd = np_state[5]

    wei_n = np.array([ue * np.cos(fi), ue * np.sin(fi), 0.0], dtype=float)
    _, r1, r2, _ = earth_model(h=h, fi=fi, lamd=lamd)

    a = 6378245.0
    b = 6356856.0
    e2 = (a**2 - b**2) / a**2
    rn = np.array(
        [
            -r1 * e2 * np.sin(fi) * np.cos(fi),
            r1 * (1.0 - e2 * (np.sin(fi)) ** 2),
            0.0,
        ],
        dtype=float,
    )

    ap_n = -np.cross(wei_n, np.cross(wei_n, rn))
    ak_n = -2.0 * np.cross(wei_n, vbe_n)
    wge_n = np.array([vbe_n[2] / r1, vbe_n[2] * np.tan(fi) / r1, -vbe_n[0] / r2], dtype=float)
    at_n = -np.cross(wge_n, vbe_n)

    g_g, _, _, _ = earth_model(h=h, fi=fi, lamd=lamd)
    a_n = c_bn @ a_dlu

    dvbe_n = a_n + ap_n + ak_n + at_n + g_g
    dh = vbe_n[1]
    dlamd = vbe_n[2] / r1 / max(np.cos(fi), 1e-6)
    dfi = vbe_n[0] / r2
    dnp = np.array([dvbe_n[0], dvbe_n[1], dvbe_n[2], dh, dfi, dlamd], dtype=float)

    np_new = np_state + dnp * dt
    return np_new, dvbe_n


def _orientation(
    w_dus: np.ndarray,
    q: np.ndarray,
    c_bn_prev: np.ndarray,
    np_state: np.ndarray,
    dt: float,
) -> tuple[np.ndarray, np.ndarray]:
    ue = 7292115.0e-11
    vbe_n = np_state[:3].copy()
    h = np_state[3]
    fi = np_state[4]
    lamd = np_state[5]

    wei_n = np.array([ue * np.cos(fi), ue * np.sin(fi), 0.0], dtype=float)
    _, r1, r2, _ = earth_model(h=h, fi=fi, lamd=lamd)
    wne_n = np.array([vbe_n[2] / r1, vbe_n[2] * np.tan(fi) / r1, -vbe_n[0] / r2], dtype=float)
    wni_n = wne_n + wei_n

    a_nb = c_bn_prev.T
    wni_b = a_nb @ wni_n
    wbnb = w_dus - wni_b

    qqq = np.dot(q, q)
    dq0 = 0.5 * (-q[1] * wbnb[0] - q[2] * wbnb[1] - q[3] * wbnb[2] + q[0] * (1.0 - qqq))
    dq1 = 0.5 * (q[0] * wbnb[0] - q[3] * wbnb[1] + q[2] * wbnb[2] + q[1] * (1.0 - qqq))
    dq2 = 0.5 * (q[3] * wbnb[0] + q[0] * wbnb[1] - q[1] * wbnb[2] + q[2] * (1.0 - qqq))
    dq3 = 0.5 * (-q[2] * wbnb[0] + q[1] * wbnb[1] + q[0] * wbnb[2] + q[3] * (1.0 - qqq))
    q_new = np.array(
        [q[0] + dq0 * dt, q[1] + dq1 * dt, q[2] + dq2 * dt, q[3] + dq3 * dt],
        dtype=float,
    )

    c_nb = np.zeros((3, 3), dtype=float)
    c_nb[0, 0] = 1.0 - 2.0 * (q_new[2] ** 2 + q_new[3] ** 2)
    c_nb[0, 1] = 2.0 * q_new[1] * q_new[2] + 2.0 * q_new[3] * q_new[0]
    c_nb[0, 2] = 2.0 * q_new[1] * q_new[3] - 2.0 * q_new[2] * q_new[0]
    c_nb[1, 0] = 2.0 * q_new[1] * q_new[2] - 2.0 * q_new[3] * q_new[0]
    c_nb[1, 1] = 1.0 - 2.0 * (q_new[1] ** 2 + q_new[3] ** 2)
    c_nb[1, 2] = 2.0 * q_new[2] * q_new[3] + 2.0 * q_new[1] * q_new[0]
    c_nb[2, 0] = 2.0 * q_new[1] * q_new[3] + 2.0 * q_new[2] * q_new[0]
    c_nb[2, 1] = 2.0 * q_new[2] * q_new[3] - 2.0 * q_new[1] * q_new[0]
    c_nb[2, 2] = 1.0 - 2.0 * (q_new[1] ** 2 + q_new[2] ** 2)
    c_bn_new = c_nb.T
    return c_bn_new, q_new


def bins2_step(
    a_true: np.ndarray,
    w_true: np.ndarray,
    np_state: np.ndarray,
    c_bn: np.ndarray,
    q: np.ndarray,
    dt: float,
    params: Bins2Params,
    noise: ImuNoise,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    w_dus = gyro_measurement(w_true=w_true, params=params, noise=noise)
    a_dlu = accel_measurement(a_true=a_true, params=params, noise=noise)
    np_new, _ = _navigation(
        a_dlu=a_dlu,
        c_bn=c_bn,
        np_state=np_state,
        dt=dt,
    )
    c_bn_new, q_new = _orientation(
        w_dus=w_dus,
        q=q,
        c_bn_prev=c_bn,
        np_state=np_state,
        dt=dt,
    )
    return q_new, c_bn_new, np_new, a_dlu
