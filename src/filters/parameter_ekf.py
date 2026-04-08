"""
Parameter EKF: совместная оценка [V,α,q,θ] и [CLα, CLδe, Cmα, Cmq, Cmδe].
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from models.aircraft_dynamics import NonlinearLongitudinalParams, longitudinal_dynamics

STATE_DIM = 9


def _params_from_state(x: np.ndarray, base: NonlinearLongitudinalParams) -> NonlinearLongitudinalParams:
    return replace(
        base,
        cla=float(x[4]),
        cl_de=float(x[5]),
        cma=float(x[6]),
        cmq=float(x[7]),
        cmde=float(x[8]),
    )


def _rk4_step(
    x: np.ndarray, delta_e: float, dt: float, base: NonlinearLongitudinalParams
) -> np.ndarray:
    def f_full(xv: np.ndarray) -> np.ndarray:
        m = _params_from_state(xv, base)
        d4 = longitudinal_dynamics(xv[:4], delta_e, m)
        return np.concatenate([d4, np.zeros(5, dtype=float)])

    k1 = f_full(x)
    k2 = f_full(x + 0.5 * dt * k1)
    k3 = f_full(x + 0.5 * dt * k2)
    k4 = f_full(x + dt * k3)
    x_new = x + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    x_new[0] = float(np.clip(x_new[0], base.v_clip_min, base.v_clip_max))
    return x_new


def measurement_h(
    x: np.ndarray,
    delta_e: float,
    base: NonlinearLongitudinalParams,
    *,
    include_nz: bool = False,
) -> np.ndarray:
    v, alpha, q, theta = float(x[0]), float(x[1]), float(x[2]), float(x[3])
    if not include_nz:
        return np.array([v, alpha, q, theta], dtype=float)
    mdl = _params_from_state(x, base)
    rho, S, g, m = mdl.rho, mdl.S, mdl.g, mdl.m
    qb = 0.5 * rho * max(v, 1.0) ** 2
    cl = mdl.cl0 + mdl.cla * alpha + mdl.cla2 * alpha**2 + mdl.cl_de * delta_e
    lift = qb * S * cl
    nz = lift / (m * g)
    return np.array([v, alpha, q, theta, nz], dtype=float)


def _jacobian_F(
    x: np.ndarray,
    delta_e: float,
    dt: float,
    base: NonlinearLongitudinalParams,
    eps: float = 1e-5,
) -> np.ndarray:
    def Phi(xv: np.ndarray) -> np.ndarray:
        return _rk4_step(xv, delta_e, dt, base)

    f0 = Phi(x)
    jac = np.zeros((STATE_DIM, STATE_DIM), dtype=float)
    for j in range(STATE_DIM):
        dx = np.zeros(STATE_DIM, dtype=float)
        dx[j] = eps
        jac[:, j] = (Phi(x + dx) - f0) / eps
    return jac


def _jacobian_H(
    x: np.ndarray,
    delta_e: float,
    base: NonlinearLongitudinalParams,
    *,
    include_nz: bool,
    meas_dim: int,
    eps: float = 1e-5,
) -> np.ndarray:
    h0 = measurement_h(x, delta_e, base, include_nz=include_nz)
    jac = np.zeros((meas_dim, STATE_DIM), dtype=float)
    for j in range(STATE_DIM):
        dx = np.zeros(STATE_DIM, dtype=float)
        dx[j] = eps
        jac[:, j] = (measurement_h(x + dx, delta_e, base, include_nz=include_nz) - h0) / eps
    return jac


@dataclass(frozen=True)
class ParameterEkfConfig:
    q_diag: np.ndarray
    r_diag: np.ndarray
    include_nz: bool = False


def parameter_ekf_step(
    x: np.ndarray,
    p_mat: np.ndarray,
    y: np.ndarray,
    delta_e: float,
    dt: float,
    base: NonlinearLongitudinalParams,
    cfg: ParameterEkfConfig,
) -> tuple[np.ndarray, np.ndarray]:
    qd = np.asarray(cfg.q_diag, dtype=float)
    q_mat = np.diag(np.where(qd > 0.0, np.maximum(qd, 1e-24), 0.0))
    r_mat = np.diag(np.maximum(cfg.r_diag, 1e-20))
    incl = bool(cfg.include_nz)
    meas_dim = int(y.shape[0])
    if meas_dim != (5 if incl else 4):
        raise ValueError("Длина y и include_nz/r_diag несогласованы")

    x_pred = _rk4_step(x, delta_e, dt, base)
    f_disc = _jacobian_F(x, delta_e, dt, base)
    p_pred = f_disc @ p_mat @ f_disc.T + q_mat

    h_pred = measurement_h(x_pred, delta_e, base, include_nz=incl)
    innov = y - h_pred
    h_mat = _jacobian_H(x_pred, delta_e, base, include_nz=incl, meas_dim=meas_dim)
    s_mat = h_mat @ p_pred @ h_mat.T + r_mat
    ph = p_pred @ h_mat.T
    k_gain = np.linalg.solve(s_mat, ph.T).T
    x_new = x_pred + k_gain @ innov
    i_kh = np.eye(STATE_DIM, dtype=float) - k_gain @ h_mat
    p_new = i_kh @ p_pred @ i_kh.T + k_gain @ r_mat @ k_gain.T
    x_new[0] = float(np.clip(x_new[0], base.v_clip_min, base.v_clip_max))
    return x_new, p_new
