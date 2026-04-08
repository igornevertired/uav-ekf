from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class LongitudinalState:
    """x_true = [V, alpha, q, theta]^T"""

    v: float
    alpha: float
    q: float
    theta: float

    def as_vector(self) -> np.ndarray:
        return np.array([self.v, self.alpha, self.q, self.theta], dtype=float)

    @staticmethod
    def from_vector(x: np.ndarray) -> "LongitudinalState":
        return LongitudinalState(
            v=float(x[0]), alpha=float(x[1]), q=float(x[2]), theta=float(x[3])
        )


@dataclass(frozen=True)
class NonlinearLongitudinalParams:
    m: float
    g: float
    rho: float
    S: float
    cbar: float
    Iyy: float
    thrust: float
    auto_trim_thrust_along_path: bool
    cl0: float
    cla: float
    cla2: float
    cl_de: float
    cd0: float
    k_induced: float
    cm0: float
    cma: float
    cmq: float
    cmde: float
    q_damping: float
    v_clip_min: float
    v_clip_max: float

    @staticmethod
    def default() -> "NonlinearLongitudinalParams":
        return NonlinearLongitudinalParams(
            m=18.0,
            g=9.80665,
            rho=1.225,
            S=0.55,
            cbar=0.18,
            Iyy=10.0,
            thrust=95.0,
            auto_trim_thrust_along_path=True,
            cl0=0.28,
            cla=4.5,
            cla2=-1.2,
            cl_de=0.0,
            cd0=0.028,
            k_induced=0.045,
            cm0=0.015,
            cma=-0.48,
            cmq=-12.0,
            cmde=-0.55,
            q_damping=0.55,
            v_clip_min=12.0,
            v_clip_max=90.0,
        )

    @staticmethod
    def default_initial_state() -> np.ndarray:
        return np.array([30.0, 0.07, 0.0, 0.07], dtype=float)


def _qbar(rho: float, v: float) -> float:
    v_safe = max(float(v), 1.0)
    return 0.5 * rho * v_safe**2


def _cl(alpha: float, delta_e: float, p: NonlinearLongitudinalParams) -> float:
    return p.cl0 + p.cla * alpha + p.cla2 * alpha**2 + p.cl_de * delta_e


def _cd(cl: float, p: NonlinearLongitudinalParams) -> float:
    return p.cd0 + p.k_induced * cl**2


def _cm(alpha: float, q_hat: float, delta_e: float, p: NonlinearLongitudinalParams) -> float:
    return p.cm0 + p.cma * alpha + p.cmq * q_hat + p.cmde * delta_e


def equilibrium_thrust_for_vdot_zero(
    x_true: np.ndarray, model: NonlinearLongitudinalParams
) -> float:
    v, alpha, _, theta = float(x_true[0]), float(x_true[1]), float(x_true[2]), float(x_true[3])
    gamma = theta - alpha
    qb = _qbar(model.rho, v)
    cl = _cl(alpha, 0.0, model)
    cd = _cd(cl, model)
    drag = qb * model.S * cd
    return float((drag + model.m * model.g * np.sin(gamma)) / max(np.cos(alpha), 0.2))


def longitudinal_dynamics(
    x_true: np.ndarray, delta_e: float, model: NonlinearLongitudinalParams
) -> np.ndarray:
    v, alpha, q, theta = float(x_true[0]), float(x_true[1]), float(x_true[2]), float(x_true[3])
    gamma = theta - alpha
    v_path = max(v, 12.0)

    qb = _qbar(model.rho, v)
    cl = _cl(alpha, delta_e, model)
    cd = _cd(cl, model)
    lift = qb * model.S * cl
    drag = qb * model.S * cd

    if model.auto_trim_thrust_along_path:
        thrust = (drag + model.m * model.g * np.sin(gamma)) / max(np.cos(alpha), 0.2)
    else:
        thrust = model.thrust

    q_hat = model.cbar * q / (2.0 * max(v, 1.0))
    cm = _cm(alpha, q_hat, delta_e, model)
    moment = qb * model.S * model.cbar * cm

    v_dot = (thrust * np.cos(alpha) - drag - model.m * model.g * np.sin(gamma)) / model.m
    gamma_dot = (thrust * np.sin(alpha) + lift - model.m * model.g * np.cos(gamma)) / (
        model.m * v_path
    )
    alpha_dot = q - gamma_dot
    theta_dot = q
    q_dot = moment / model.Iyy - model.q_damping * q

    return np.array([v_dot, alpha_dot, q_dot, theta_dot], dtype=float)


def aircraft_longitudinal_step(
    x_true: np.ndarray,
    delta_e: float,
    dt: float,
    model: NonlinearLongitudinalParams | None = None,
) -> np.ndarray:
    m = model or NonlinearLongitudinalParams.default()
    k1 = longitudinal_dynamics(x_true, delta_e, m)
    k2 = longitudinal_dynamics(x_true + 0.5 * dt * k1, delta_e, m)
    k3 = longitudinal_dynamics(x_true + 0.5 * dt * k2, delta_e, m)
    k4 = longitudinal_dynamics(x_true + dt * k3, delta_e, m)
    x_new = x_true + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    x_new[0] = float(np.clip(x_new[0], m.v_clip_min, m.v_clip_max))
    return x_new


def load_factor_nz(
    x_true: np.ndarray, delta_e: float, model: NonlinearLongitudinalParams
) -> float:
    """n_z ≈ L/(m g)."""
    v, alpha = float(x_true[0]), float(x_true[1])
    qb = _qbar(model.rho, v)
    cl = _cl(alpha, delta_e, model)
    lift = qb * model.S * cl
    return float(lift / (model.m * model.g))


def compute_specific_force_body(
    x_true: np.ndarray, delta_e: float, model: NonlinearLongitudinalParams
) -> np.ndarray:
    """
    Удельная сила в связанной СК [Forward, Right, Down] (без гравитации в «чистом» акселерометре
    здесь — аэродинамика+тяга)/m.
    """
    v, alpha, _, theta = float(x_true[0]), float(x_true[1]), float(x_true[2]), float(x_true[3])
    gamma = theta - alpha
    qb = _qbar(model.rho, v)
    cl = _cl(alpha, delta_e, model)
    cd = _cd(cl, model)
    lift = qb * model.S * cl
    drag = qb * model.S * cd
    if model.auto_trim_thrust_along_path:
        thrust = (drag + model.m * model.g * np.sin(gamma)) / max(np.cos(alpha), 0.2)
    else:
        thrust = model.thrust
    fx = (thrust - drag * np.cos(alpha) + lift * np.sin(alpha)) / model.m
    fy = 0.0
    fz = (-lift * np.cos(alpha) - drag * np.sin(alpha)) / model.m
    return np.array([fx, fy, fz], dtype=float)
