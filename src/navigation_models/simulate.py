"""Ядро: симуляция, нав. EKF, кинематика fusion — импортируется из simulate.py."""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from models.aircraft_dynamics import (
    NonlinearLongitudinalParams,
    aircraft_longitudinal_step,
    compute_specific_force_body,
    equilibrium_thrust_for_vdot_zero,
    load_factor_nz,
)
from src.filters.parameter_ekf import ParameterEkfConfig, parameter_ekf_step
from src.navigation_models.altimeter import altimeter_measurement, altimeter_reference_sigma
from src.navigation_models.bins import bins2_step, gyro_measurement
from src.navigation_models.earth import earth_model
from src.navigation_models.gnss import GnssNoise, gnss_measurement, gnss_reference_sigmas
from src.navigation_models.noise import NoiseStreams
from src.navigation_models.noise import q_gyro_measurement_sigma, theta_ins_measurement_sigma
from src.navigation_models.zala421_16e5 import Zala421Config, zala421_default_config


@dataclass(frozen=True)
class NavigationAccuracyProfile:
    name: str = "baseline"
    gnss_pos_scale: float = 1.0
    gnss_vel_scale: float = 1.0
    altimeter_scale: float = 1.0
    theta_scale: float = 1.0
    q_scale: float = 1.0


@dataclass(frozen=True)
class IdentificationWindow:
    t_start: float = 5.0
    t_end: float = 30.0


def pitch_from_c_bn(c_bn: np.ndarray) -> float:
    """Тангаж θ из DCM тело→НСК [N,U,E]."""
    return float(np.arctan2(c_bn[1, 0], c_bn[0, 0]))


def _delta_angle(a0: float, a1: float) -> float:
    d = float(a1 - a0)
    return float((d + np.pi) % (2.0 * np.pi) - np.pi)


def elevator_step_command_deg(t: float, t_model: float) -> float:
    schedule = [
        (0.0, 5.0, 0.0),
        (5.0, 8.0, -1.0),
        (8.0, 11.0, 1.0),
        (11.0, 13.0, 0.0),
        (13.0, 15.0, -1.0),
        (15.0, 17.0, 1.0),
        (17.0, 18.0, 0.0),
        (18.0, 19.0, -1.0),
        (19.0, 20.0, 1.0),
    ]
    for t0, t1, level in schedule:
        if t0 <= t < t1:
            return float(level)
    return 0.0


def q_excitation_command_deg(t: float) -> float:
    """Дополнительный короткий doublet для возбуждения канала q."""
    schedule = [
        (5.0, 6.2, -2.5),
        (6.2, 7.4, 2.5),
        (13.0, 14.2, -2.5),
        (14.2, 15.4, 2.5),
        (18.0, 18.8, -2.0),
        (18.8, 19.6, 2.0),
    ]
    for t0, t1, level in schedule:
        if t0 <= t < t1:
            return float(level)
    return 0.0


def nav_ekf_initial_sigmas() -> dict[str, float]:
    """Начальные 1σ ошибки из таблицы на скриншоте, в проекции на наш nav-state."""
    return {
        "phi": 1.57e-7,
        "lambda": 1.57e-7,
        "vn": 0.02,
        "ve": 0.02,
        "vh": 0.02,
        "h": 1.0,
        # В 6D модели оставлен только тангаж; берём малую стартовую ошибку БИНС по pitch.
        "theta": 1.0e-4,
    }


def scaled_gnss_sigmas(profile: NavigationAccuracyProfile) -> dict[str, float]:
    ref = gnss_reference_sigmas()
    return {
        "pseudorange": ref["pseudorange"] * float(profile.gnss_pos_scale),
        "pseudorange_rate": ref["pseudorange_rate"] * float(profile.gnss_vel_scale),
        "phi": ref["phi"] * float(profile.gnss_pos_scale),
        "lambda": ref["lambda"] * float(profile.gnss_pos_scale),
        "vn": ref["vn"] * float(profile.gnss_vel_scale),
        "ve": ref["ve"] * float(profile.gnss_vel_scale),
    }


def scaled_altimeter_sigma(profile: NavigationAccuracyProfile) -> float:
    return float(altimeter_reference_sigma() * float(profile.altimeter_scale))


def scaled_theta_ins_sigma(profile: NavigationAccuracyProfile) -> float:
    return float(theta_ins_measurement_sigma() * float(profile.theta_scale))


def scaled_q_gyro_sigma(profile: NavigationAccuracyProfile) -> float:
    return float(q_gyro_measurement_sigma() * float(profile.q_scale))


def nz_from_specific_force_measurement(a_body_hist: np.ndarray, g: float = 9.80665) -> np.ndarray:
    acc = np.asarray(a_body_hist, dtype=float)
    if acc.ndim != 2 or acc.shape[1] != 3:
        raise ValueError("a_body_hist должен иметь размерность [N,3]")
    return -acc[:, 2] / float(g)


def identification_mask(time_hist: np.ndarray, window: IdentificationWindow) -> np.ndarray:
    t = np.asarray(time_hist, dtype=float)
    return (t >= float(window.t_start)) & (t <= float(window.t_end))


def longitudinal_measurement_from_nav_state(x_nav: np.ndarray) -> np.ndarray:
    vn = float(x_nav[2])
    ve = float(x_nav[3])
    vh = float(x_nav[4])
    theta = float(x_nav[6])
    q = float(x_nav[7])
    r_h = float(np.sqrt(max(vn * vn + ve * ve, 1e-12)))
    v_mag = float(np.sqrt(max(vn * vn + ve * ve + vh * vh, 1e-12)))
    gamma = float(np.arctan2(vh, r_h))
    alpha = float(theta - gamma)
    return np.array([v_mag, alpha, q, theta], dtype=float)


def longitudinal_measurement_jacobian_from_nav(x_nav: np.ndarray) -> np.ndarray:
    vn = float(x_nav[2])
    ve = float(x_nav[3])
    vh = float(x_nav[4])
    v_mag = float(np.sqrt(max(vn * vn + ve * ve + vh * vh, 1e-12)))
    r_h = float(np.sqrt(max(vn * vn + ve * ve, 1e-12)))
    den = float(max(vh * vh + r_h * r_h, 1e-12))

    jac = np.zeros((4, 8), dtype=float)
    jac[0, 2] = vn / v_mag
    jac[0, 3] = ve / v_mag
    jac[0, 4] = vh / v_mag

    dgamma_dvn = -vh * vn / (r_h * den)
    dgamma_dve = -vh * ve / (r_h * den)
    dgamma_dvh = r_h / den
    jac[1, 2] = -dgamma_dvn
    jac[1, 3] = -dgamma_dve
    jac[1, 4] = -dgamma_dvh
    jac[1, 6] = 1.0

    jac[2, 7] = 1.0
    jac[3, 6] = 1.0
    return jac


def longitudinal_measurement_covariance_from_nav(x_nav: np.ndarray, p_nav: np.ndarray) -> np.ndarray:
    jac = longitudinal_measurement_jacobian_from_nav(np.asarray(x_nav, dtype=float))
    p = np.asarray(p_nav, dtype=float)
    r = jac @ p @ jac.T
    r = 0.5 * (r + r.T)
    diag = np.maximum(np.diag(r), 1e-20)
    r[np.diag_indices_from(r)] = diag
    return r


def ekf_fuse_bins_gnss(
    bins_hist: np.ndarray,
    gnss_hist: np.ndarray,
    altimeter_hist: np.ndarray,
    theta_ins_hist: np.ndarray,
    q_gyro_hist: np.ndarray,
    c_bn_hist: np.ndarray,
    dt: float,
    profile: NavigationAccuracyProfile | None = None,
) -> dict[str, np.ndarray]:
    n = bins_hist.shape[0]
    state_dim = 8
    meas_dim = 7
    x_hat = np.zeros((n, state_dim), dtype=float)
    p_hist = np.zeros((n, state_dim, state_dim), dtype=float)
    sig3_hist = np.zeros((n, state_dim), dtype=float)

    nav_profile = profile or NavigationAccuracyProfile()
    gnss_sig = scaled_gnss_sigmas(nav_profile)
    gnss_pos_std = gnss_sig["phi"]
    gnss_vel_std = gnss_sig["vn"]
    alt_std = scaled_altimeter_sigma(nav_profile)
    init_sig = nav_ekf_initial_sigmas()
    theta_std = scaled_theta_ins_sigma(nav_profile)
    q_std = scaled_q_gyro_sigma(nav_profile)

    bins_nav = np.column_stack(
        [bins_hist[:, 4], bins_hist[:, 5], bins_hist[:, 0], bins_hist[:, 2]]
    )
    vh_bins = bins_hist[:, 1].astype(float)
    h_bins = bins_hist[:, 3].astype(float)

    H = np.zeros((meas_dim, state_dim), dtype=float)
    for i in range(4):
        H[i, i] = 1.0
    H[4, 5] = 1.0
    H[5, 6] = 1.0
    H[6, 7] = 1.0

    x = np.concatenate(
        [
            gnss_hist[0].copy(),
            np.array([vh_bins[0], h_bins[0], theta_ins_hist[0], q_gyro_hist[0]], dtype=float),
        ]
    )
    p = np.diag(
        [
            init_sig["phi"] ** 2,
            init_sig["lambda"] ** 2,
            init_sig["vn"] ** 2,
            init_sig["ve"] ** 2,
            init_sig["vh"] ** 2,
            init_sig["h"] ** 2,
            init_sig["theta"] ** 2,
            q_std**2,
        ]
    )

    q_vel = 1e-5
    q_pos = 1e-14
    q_vh = 2e-6
    q_h = 1e-3
    q_th = 3e-10
    q_q = 5e-7
    q_mat = np.diag([q_pos, q_pos, q_vel, q_vel, q_vh, q_h, q_th, q_q])

    r_mat = np.diag(
        [gnss_pos_std**2, gnss_pos_std**2, gnss_vel_std**2, gnss_vel_std**2, alt_std**2, theta_std**2, q_std**2]
    )

    i_state = np.eye(state_dim, dtype=float)

    for k in range(n):
        if k == 0:
            x_pred = x.copy()
            f_mat = np.eye(state_dim, dtype=float)
            f_mat[5, 4] = dt
        else:
            dphi, dlm, dvn, dve = (bins_nav[k] - bins_nav[k - 1]).tolist()
            dvh = float(vh_bins[k] - vh_bins[k - 1])
            dth = _delta_angle(theta_ins_hist[k - 1], theta_ins_hist[k])
            dq = float(q_gyro_hist[k] - q_gyro_hist[k - 1])
            x_pred = x.copy()
            x_pred[:4] = x[:4] + np.array([dphi, dlm, dvn, dve], dtype=float)
            x_pred[4] = x[4] + dvh
            x_pred[5] = x[5] + x[4] * dt
            x_pred[6] = x[6] + dth
            x_pred[7] = x[7] + dq
            f_mat = np.eye(state_dim, dtype=float)
            f_mat[5, 4] = dt

        p_pred = f_mat @ p @ f_mat.T + q_mat

        y = np.concatenate(
            [
                gnss_hist[k],
                np.array([altimeter_hist[k], theta_ins_hist[k], q_gyro_hist[k]], dtype=float),
            ]
        )
        innov = y - H @ x_pred
        s = H @ p_pred @ H.T + r_mat
        k_gain = p_pred @ H.T @ np.linalg.inv(s)
        x = x_pred + k_gain @ innov
        p = (i_state - k_gain @ H) @ p_pred

        x_hat[k] = x
        p_hist[k] = p
        sig3_hist[k] = 3.0 * np.sqrt(np.maximum(np.diag(p), 0.0))

    return {
        "x_hat": x_hat,
        "p": p_hist,
        "sig3": sig3_hist,
        "nav_state_dim": state_dim,
    }


def nominal_aero_params_for_nz(
    plant_params: NonlinearLongitudinalParams,
) -> NonlinearLongitudinalParams:
    nom = NonlinearLongitudinalParams.default()
    return replace(
        plant_params,
        cl0=nom.cl0,
        cla=nom.cla,
        cla2=nom.cla2,
        cl_de=nom.cl_de,
        cd0=nom.cd0,
        k_induced=nom.k_induced,
        cm0=nom.cm0,
        cma=nom.cma,
        cmq=nom.cmq,
        cmde=nom.cmde,
    )


def longitudinal_state_from_fusion(
    x_hat_nav: np.ndarray,
    bins_hist: np.ndarray,
    c_bn_hist: np.ndarray,
    dt: float,
) -> np.ndarray:
    n = x_hat_nav.shape[0]
    vn = x_hat_nav[:, 2]
    ve = x_hat_nav[:, 3]
    if x_hat_nav.shape[1] >= 7:
        vh = x_hat_nav[:, 4]
        theta_hat = x_hat_nav[:, 6]
    elif x_hat_nav.shape[1] >= 6:
        vh = x_hat_nav[:, 4]
        theta_hat = x_hat_nav[:, 5]
    else:
        vh = bins_hist[:, 1]
        theta_hat = np.array([pitch_from_c_bn(c_bn_hist[k]) for k in range(n)], dtype=float)
    r_h = np.sqrt(np.maximum(vn**2 + ve**2, 1e-12))
    v_mag = np.sqrt(np.maximum(vn**2 + ve**2 + vh**2, 1e-12))
    gamma_hat = np.arctan2(vh, r_h)
    alpha_hat = theta_hat - gamma_hat
    if x_hat_nav.shape[1] >= 8:
        q_hat = x_hat_nav[:, 7]
    else:
        q_hat = np.gradient(theta_hat, dt)
    return np.column_stack([v_mag, alpha_hat, q_hat, theta_hat])


upper_stage_kinematics_for_param_ekf = longitudinal_state_from_fusion


def aero_param_initial_guess(p_true: np.ndarray) -> np.ndarray:
    factors = np.array([1.6, 0.4, 1.6, 0.4, 1.6], dtype=float)
    return np.asarray(p_true, dtype=float) * factors


def state5_from_nav_solution(
    x_hat_nav: np.ndarray,
    bins_hist: np.ndarray,
    c_bn_hist: np.ndarray,
    dt: float,
) -> np.ndarray:
    long4 = np.vstack([longitudinal_measurement_from_nav_state(x_hat_nav[k]) for k in range(x_hat_nav.shape[0])])
    if x_hat_nav.shape[1] >= 8:
        h_hat = x_hat_nav[:, 5]
    elif x_hat_nav.shape[1] >= 7:
        h_hat = x_hat_nav[:, 5]
    else:
        h_hat = bins_hist[:, 3]
    return np.column_stack([long4[:, 0], long4[:, 1], long4[:, 3], long4[:, 2], h_hat])


def state5_from_sensor_measurements(
    gnss_hist: np.ndarray,
    bins_hist: np.ndarray,
    theta_ins_hist: np.ndarray,
    q_gyro_hist: np.ndarray,
    altimeter_hist: np.ndarray,
    dt: float,
) -> np.ndarray:
    vn = gnss_hist[:, 2]
    ve = gnss_hist[:, 3]
    vh = bins_hist[:, 1]
    r_h = np.sqrt(np.maximum(vn**2 + ve**2, 1e-12))
    v_mag = np.sqrt(np.maximum(vn**2 + ve**2 + vh**2, 1e-12))
    gamma = np.arctan2(vh, r_h)
    alpha = theta_ins_hist - gamma
    return np.column_stack([v_mag, alpha, theta_ins_hist, q_gyro_hist, altimeter_hist])


def sigma3_state5_from_nav_covariance(
    x_hat_nav: np.ndarray,
    bins_hist: np.ndarray,
    p_hist: np.ndarray,
    dt: float,
    profile: NavigationAccuracyProfile | None = None,
) -> np.ndarray:
    nav_profile = profile or NavigationAccuracyProfile()
    gnss_sig = scaled_gnss_sigmas(nav_profile)
    long_sig3 = sigma3_longitudinal_state(
        x_hat_nav,
        bins_hist,
        p_hist,
        gnss_vel_std=gnss_sig["vn"],
        sigma_theta_ins=scaled_theta_ins_sigma(nav_profile),
        dt=dt,
        sigma_vh=nav_ekf_initial_sigmas()["vh"],
        sigma_q_gyro=scaled_q_gyro_sigma(nav_profile),
    )
    h_sig3 = 3.0 * np.sqrt(np.maximum(p_hist[:, 5, 5], 0.0))
    th_sig3 = 3.0 * np.sqrt(np.maximum(p_hist[:, 6, 6], 0.0))
    if p_hist.shape[1] >= 8:
        q_sig3 = 3.0 * np.sqrt(np.maximum(p_hist[:, 7, 7], 0.0))
    else:
        q_sig3 = long_sig3[:, 2]
    return np.column_stack([long_sig3[:, 0], long_sig3[:, 1], th_sig3, q_sig3, h_sig3])


def state_error_series(truth: np.ndarray, estimate: np.ndarray) -> np.ndarray:
    return np.asarray(estimate, dtype=float) - np.asarray(truth, dtype=float)


def q_dot_from_q_series(q_hist: np.ndarray, dt: float) -> np.ndarray:
    q = np.asarray(q_hist, dtype=float)
    if q.size == 0:
        return np.zeros(0, dtype=float)
    if q.size < 5:
        return np.gradient(q, max(float(dt), 1e-9))
    kernel = np.array([1.0, 2.0, 3.0, 2.0, 1.0], dtype=float)
    kernel /= kernel.sum()
    q_smooth = np.convolve(q, kernel, mode="same")
    return np.gradient(q_smooth, max(float(dt), 1e-9))


def navigation_outputs_to_parameter_measurements(
    x_hat_nav: np.ndarray,
    p_nav_hist: np.ndarray,
    *,
    include_nz: bool,
    include_qdot: bool,
    nz_hist: np.ndarray | None,
    delta_e_hist: np.ndarray,
    dt: float,
    nz_var: float = (5e-3) ** 2,
) -> tuple[np.ndarray, np.ndarray]:
    n = x_hat_nav.shape[0]
    meas_dim = 4 + int(include_nz) + int(include_qdot)
    y_hist = np.zeros((n, meas_dim), dtype=float)
    r_hist = np.zeros((n, meas_dim, meas_dim), dtype=float)
    q_dot_hist = q_dot_from_q_series(x_hat_nav[:, 7], dt)
    q_var = np.maximum(p_nav_hist[:, 7, 7], 1e-20)
    dt_eff = max(5.0 * float(dt), 1e-9)
    q_dot_var = np.maximum((2.0 * q_var) / dt_eff**2, 0.12**2)
    for k in range(n):
        y_base = longitudinal_measurement_from_nav_state(x_hat_nav[k])
        r_base = longitudinal_measurement_covariance_from_nav(x_hat_nav[k], p_nav_hist[k])
        y_parts = [y_base]
        r_blocks = [r_base]
        if include_nz:
            if nz_hist is None:
                raise ValueError("Для include_nz=True требуется nz_hist")
            y_parts.append(np.array([float(nz_hist[k])], dtype=float))
            r_blocks.append(np.array([[nz_var]], dtype=float))
        if include_qdot:
            y_parts.append(np.array([q_dot_hist[k]], dtype=float))
            r_blocks.append(np.array([[q_dot_var[k]]], dtype=float))
        y_hist[k] = np.concatenate(y_parts)
        r_hist[k] = np.zeros((meas_dim, meas_dim), dtype=float)
        idx = 0
        for block in r_blocks:
            block_dim = block.shape[0]
            r_hist[k, idx : idx + block_dim, idx : idx + block_dim] = block
            idx += block_dim
    return y_hist, r_hist


def sigma3_longitudinal_state(
    x_hat_nav: np.ndarray,
    bins_hist: np.ndarray,
    p_hist: np.ndarray,
    gnss_vel_std: float,
    sigma_theta_ins: float,
    dt: float,
    *,
    sigma_vh: float | None = None,
    sigma_q_gyro: float | None = None,
    inflation: float = 1.08,
) -> np.ndarray:
    n = x_hat_nav.shape[0]
    sig3 = np.zeros((n, 4), dtype=float)
    sig_vh = float(gnss_vel_std if sigma_vh is None else sigma_vh)
    sigma_q_val = float(q_gyro_measurement_sigma() if sigma_q_gyro is None else sigma_q_gyro)
    nav7 = x_hat_nav.shape[1] >= 7 and p_hist.shape[1] >= 7
    nav6 = x_hat_nav.shape[1] >= 6 and p_hist.shape[1] >= 6

    for k in range(n):
        vn, ve = float(x_hat_nav[k, 2]), float(x_hat_nav[k, 3])
        vh = float(x_hat_nav[k, 4]) if (nav7 or nav6) else float(bins_hist[k, 1])
        p = p_hist[k]
        p22, p33 = float(p[2, 2]), float(p[3, 3])
        if nav7:
            p44, p66 = float(p[4, 4]), float(p[6, 6])
            sig_th_k = np.sqrt(max(p66, sigma_theta_ins**2))
        elif nav6:
            p44, p55 = float(p[4, 4]), float(p[5, 5])
            sig_th_k = np.sqrt(max(p55, sigma_theta_ins**2))
        else:
            p44, p55 = sig_vh**2, sigma_theta_ins**2
            sig_th_k = sigma_theta_ins
        v = float(np.sqrt(max(vn * vn + ve * ve + vh * vh, 1e-12)))
        dvn = vn / v
        dve = ve / v
        dvh = vh / v
        var_v = dvn**2 * p22 + dve**2 * p33 + dvh**2 * p44
        sig3[k, 0] = 3.0 * np.sqrt(max(var_v, 0.0))

        r_h = float(np.sqrt(max(vn**2 + ve**2, 1e-12)))
        den = vh * vh + r_h * r_h
        dg_dvn = -vh * vn / (r_h * den) if den > 0 else 0.0
        dg_dve = -vh * ve / (r_h * den) if den > 0 else 0.0
        dg_dvh = r_h / den if den > 0 else 0.0
        var_gamma = dg_dvn**2 * p22 + dg_dve**2 * p33 + dg_dvh**2 * p44
        if nav7:
            var_alpha = sig_th_k**2 + var_gamma
        elif nav6:
            var_alpha = sig_th_k**2 + var_gamma
        else:
            var_alpha = sigma_theta_ins**2 + var_gamma
        sig3[k, 1] = 3.0 * np.sqrt(max(var_alpha, 0.0))

        sig3[k, 3] = 3.0 * sig_th_k
        sig_q_k = np.sqrt((np.sqrt(2.0) * sig_th_k / max(dt, 1e-9)) ** 2 + sigma_q_val**2)
        sig3[k, 2] = 3.0 * sig_q_k

    sig3 *= float(inflation)
    return sig3

def run_simulation(
    config: Zala421Config | None = None,
    t_model: float = 120.0,
    *,
    auto_trim_thrust: bool = False,
    thrust_newtons: float | None = None,
    param_ekf_with_nz: bool = True,
    profile: NavigationAccuracyProfile | None = None,
    ident_window: IdentificationWindow | None = None,
) -> dict[str, np.ndarray | dict]:
    cfg = config or zala421_default_config()
    nav_profile = profile or NavigationAccuracyProfile()
    ident_cfg = ident_window or IdentificationWindow()
    steps = int(t_model / cfg.dt)
    base_dyn = NonlinearLongitudinalParams.default()
    plant_params = replace(
        base_dyn,
        cla=5.08,
        cl_de=0.022,
        cma=-0.49,
        cmq=-12.6,
        cmde=-0.51,
    )
    x_true = NonlinearLongitudinalParams.default_initial_state().copy()
    if auto_trim_thrust:
        dyn_params = plant_params
        thrust_used = float("nan")
    else:
        t_fixed = (
            float(thrust_newtons)
            if thrust_newtons is not None
            else equilibrium_thrust_for_vdot_zero(x_true, plant_params)
        )
        dyn_params = replace(plant_params, auto_trim_thrust_along_path=False, thrust=t_fixed)
        thrust_used = t_fixed

    np_state = cfg.np0.copy()
    q = cfg.q0.copy()
    c_bn = cfg.c_bn0.copy()
    noise_streams = NoiseStreams.from_seed(42116)

    bins_hist = np.zeros((steps, 6), dtype=float)
    c_bn_hist = np.zeros((steps, 3, 3), dtype=float)
    gnss_hist = np.zeros((steps, 4), dtype=float)
    altimeter_hist = np.zeros(steps, dtype=float)
    theta_ins_hist = np.zeros(steps, dtype=float)
    q_gyro_hist = np.zeros(steps, dtype=float)
    accel_meas_hist = np.zeros((steps, 3), dtype=float)
    nz_meas_hist = np.zeros(steps, dtype=float)
    t_hist = np.arange(steps, dtype=float) * cfg.dt
    x_true_hist = np.zeros((steps, 4), dtype=float)
    delta_e_ref_hist = np.zeros(steps, dtype=float)
    delta_e_cmd_hist = np.zeros(steps, dtype=float)
    truth_nav_hist = np.zeros((steps, 8), dtype=float)

    fi_true = float(cfg.np0[4])
    lm_true = float(cfg.np0[5])
    h_true = float(cfg.np0[3])
    h_ref = h_true
    theta_ref = float(x_true[3])

    for i in range(steps):
        t = t_hist[i]
        delta_e_pulse_deg = elevator_step_command_deg(float(t), float(t_model))
        delta_e_q_exc_deg = q_excitation_command_deg(float(t))
        v_ctrl = float(x_true[0])
        alpha_ctrl = float(x_true[1])
        q_ctrl = float(x_true[2])
        theta_ctrl = float(x_true[3])
        vh_ctrl = v_ctrl * np.sin(theta_ctrl - alpha_ctrl)
        h_err = h_ref - h_true
        excitation_active = abs(delta_e_q_exc_deg) > 0.0
        q_hold_gain = 0.0 if excitation_active else 2.0e-1
        theta_hold_gain = 2.0e-1 if excitation_active else 8.0e-1
        delta_e_hold = (
            1.0e-3 * h_err
            + 5.0e-2 * vh_ctrl
            + theta_hold_gain * (theta_ctrl - theta_ref)
            + q_hold_gain * q_ctrl
        )
        delta_e_ref_deg = delta_e_pulse_deg + delta_e_q_exc_deg
        delta_e = np.deg2rad(delta_e_ref_deg) + delta_e_hold
        delta_e = float(np.clip(delta_e, np.deg2rad(-8.0), np.deg2rad(8.0)))
        delta_e_ref_hist[i] = np.deg2rad(delta_e_ref_deg)
        delta_e_cmd_hist[i] = delta_e

        a_true_body = compute_specific_force_body(x_true, delta_e, dyn_params)
        x_true = aircraft_longitudinal_step(x_true=x_true, delta_e=delta_e, dt=cfg.dt, model=dyn_params)
        x_true_hist[i] = x_true

        v_true = float(x_true[0])
        alpha_true = float(x_true[1])
        q_pitch = float(x_true[2])
        theta_true = float(x_true[3])
        gamma_true = theta_true - alpha_true

        vn_true = v_true * np.cos(gamma_true)
        vh_true = v_true * np.sin(gamma_true)
        ve_true = 0.0

        _, r1, r2, _ = earth_model(h=h_true, fi=fi_true, lamd=lm_true)
        fi_true += (vn_true / r2) * cfg.dt
        lm_true += (ve_true / (r1 * max(np.cos(fi_true), 1e-6))) * cfg.dt
        h_true += vh_true * cfg.dt

        truth_nav_hist[i] = np.array(
            [fi_true, lm_true, vn_true, ve_true, vh_true, h_true, theta_true, q_pitch],
            dtype=float,
        )

        w_true_body = np.array([0.0, q_pitch, 0.0], dtype=float)
        imu_noise = noise_streams.sample_imu(cfg.dt)
        w_gyro = gyro_measurement(w_true=w_true_body, params=cfg.bins_params, noise=imu_noise)
        q, c_bn, np_state, a_dlu = bins2_step(
            a_true=a_true_body,
            w_true=w_true_body,
            np_state=np_state,
            c_bn=c_bn,
            q=q,
            dt=cfg.dt,
            params=cfg.bins_params,
            noise=imu_noise,
        )
        bins_hist[i] = np_state
        c_bn_hist[i] = c_bn
        theta_noise_gain = max(float(nav_profile.theta_scale) - 1.0, 0.0)
        theta_ins_hist[i] = float(pitch_from_c_bn(c_bn) + theta_noise_gain * noise_streams.sample_theta_ins())
        q_gyro_hist[i] = float(w_gyro[1] + nav_profile.q_scale * noise_streams.sample_q_gyro())
        accel_meas_hist[i] = a_dlu
        nz_meas_hist[i] = float(nz_from_specific_force_measurement(a_dlu[np.newaxis, :])[0])

        wn_fi, wn_lm, wn_vn, wn_ve = noise_streams.sample_gnss()
        z, _ = gnss_measurement(
            fi=fi_true,
            lm=lm_true,
            vn=vn_true,
            ve=ve_true,
            noise=GnssNoise(
                wn_fi=wn_fi * nav_profile.gnss_pos_scale,
                wn_lm=wn_lm * nav_profile.gnss_pos_scale,
                wn_vn=wn_vn * nav_profile.gnss_vel_scale,
                wn_ve=wn_ve * nav_profile.gnss_vel_scale,
            ),
            scheme=cfg.gnss_scheme,
        )
        gnss_hist[i] = z
        h_meas, _ = altimeter_measurement(h_true, nav_profile.altimeter_scale * noise_streams.sample_altimeter())
        altimeter_hist[i] = h_meas

    ekf_out = ekf_fuse_bins_gnss(
        bins_hist=bins_hist,
        gnss_hist=gnss_hist,
        altimeter_hist=altimeter_hist,
        theta_ins_hist=theta_ins_hist,
        q_gyro_hist=q_gyro_hist,
        c_bn_hist=c_bn_hist,
        dt=cfg.dt,
        profile=nav_profile,
    )
    ekf_out["delta_truth_minus_hat"] = truth_nav_hist - ekf_out["x_hat"]

    pe_hist = np.zeros((steps, 9), dtype=float)
    include_qdot = False
    meas_dim = 4 + int(param_ekf_with_nz) + int(include_qdot)
    p_true = np.array(
        [
            plant_params.cla,
            plant_params.cl_de,
            plant_params.cma,
            plant_params.cmq,
            plant_params.cmde,
        ],
        dtype=float,
    )
    x_pe = np.zeros(9, dtype=float)
    p_init = aero_param_initial_guess(p_true)
    x_pe[4:9] = p_init.copy()
    x0_pe = np.concatenate([np.zeros(4, dtype=float), p_init.copy()])
    p_pe = np.diag(
        [
            0.5**2,
            0.02**2,
            0.02**2,
            0.02**2,
            0.6**2,
            0.0025**2,
            0.15**2,
            2.0**2,
            0.15**2,
        ]
    )
    pe_cfg = ParameterEkfConfig(
        q_diag=np.array(
            [8e-5, 2e-6, 2e-6, 2e-6, 0.0, 0.0, 0.0, 0.0, 0.0],
            dtype=float,
        ),
        r_diag=np.zeros(meas_dim, dtype=float),
        include_nz=bool(param_ekf_with_nz),
        include_qdot=bool(include_qdot),
    )
    x_hat_nav = ekf_out["x_hat"]
    y_pe_hist, r_meas_hist = navigation_outputs_to_parameter_measurements(
        x_hat_nav,
        np.asarray(ekf_out["p"], dtype=float),
        include_nz=bool(param_ekf_with_nz),
        include_qdot=bool(include_qdot),
        nz_hist=nz_meas_hist,
        delta_e_hist=delta_e_cmd_hist,
        dt=cfg.dt,
    )
    x_pe[:4] = y_pe_hist[0, :4].copy()
    x0_pe[:4] = x_pe[:4].copy()
    ident_mask = identification_mask(t_hist, ident_cfg)
    for k in range(steps):
        if ident_mask[k]:
            de = float(delta_e_cmd_hist[k])
            yk = y_pe_hist[k].copy()
            x_pe, p_pe = parameter_ekf_step(
                x_pe,
                p_pe,
                yk,
                de,
                cfg.dt,
                plant_params,
                pe_cfg,
                r_override=r_meas_hist[k],
            )
        pe_hist[k] = x_pe

    p_hat_hist = pe_hist[:, 4:9]
    param_error = p_hat_hist - p_true[np.newaxis, :]
    param_rmse = np.sqrt(np.mean(param_error**2, axis=0))
    param_mae = np.mean(np.abs(param_error), axis=0)
    param_bias_end = param_error[-1, :]
    p_hat_final = pe_hist[-1, 4:9].copy()
    with np.errstate(divide="ignore", invalid="ignore"):
        param_rel_rmse = param_rmse / (np.abs(p_true) + 1e-12)
        param_rel_err_pct_initial = 100.0 * np.abs(p_init - p_true) / (np.abs(p_true) + 1e-12)
        param_rel_err_pct_final = 100.0 * np.abs(p_hat_final - p_true) / (np.abs(p_true) + 1e-12)

    return {
        "time": t_hist,
        "x_true": x_true_hist,
        "delta_e_ref": delta_e_ref_hist,
        "delta_e_cmd": delta_e_cmd_hist,
        "truth_nav": truth_nav_hist,
        "bins": bins_hist,
        "c_bn": c_bn_hist,
        "gnss": gnss_hist,
        "altimeter": altimeter_hist,
        "theta_ins": theta_ins_hist,
        "q_gyro": q_gyro_hist,
        "nz_meas": nz_meas_hist,
        "accel_meas": accel_meas_hist,
        "ekf": ekf_out,
        "plant_params": plant_params,
        "p_true": p_true,
        "param_ekf": {
            "x_hat": pe_hist,
            "y": y_pe_hist,
            "p_cov": p_pe,
            "x0": x0_pe,
            "r_meas": r_meas_hist,
            "ident_mask": ident_mask,
        },
        "param_error": param_error,
        "param_metrics": {
            "rmse": param_rmse,
            "mae": param_mae,
            "bias_final": param_bias_end,
            "rel_rmse": param_rel_rmse,
            "p_initial": p_init,
            "p_hat_final": p_hat_final,
            "rel_error_pct_initial": param_rel_err_pct_initial,
            "rel_error_pct_final": param_rel_err_pct_final,
        },
        "thrust_fixed_n": thrust_used,
        "auto_trim_thrust": auto_trim_thrust,
        "param_ekf_with_nz": param_ekf_with_nz,
        "nav_profile": nav_profile,
        "ident_window": ident_cfg,
    }


def _fill_sigma_bands_gradient(
    ax,
    t: np.ndarray,
    y: np.ndarray,
    sig: np.ndarray,
    *,
    color: str = "C0",
    n_levels: int = 3,
) -> None:
    """Полосы ±k·σ/3, k=1..n_levels — заливка с убывающей непрозрачностью (градиент по σ)."""
    ax.plot(t, y, color=color, lw=1.0)
    base_alpha = 0.22
    for k in range(1, n_levels + 1):
        frac = k / float(n_levels)
        lo = y - frac * sig
        hi = y + frac * sig
        a = base_alpha * (1.0 - (k - 1) * 0.25)
        ax.fill_between(t, lo, hi, color=color, alpha=max(a, 0.04), linewidth=0)


def state_plot_style(name: str) -> dict[str, float | bool]:
    highlight = name in {"θ", "q"}
    return {
        "model_on_top": highlight,
        "sigma_on_top": highlight,
        "sigma_alpha": 0.2 if highlight else 0.12,
        "error_smoothing_window": 3 if name == "q" else 1,
    }


def smooth_series_for_display(series: np.ndarray, window: int) -> np.ndarray:
    data = np.asarray(series, dtype=float)
    win = max(int(window), 1)
    if win <= 1 or data.size < 3:
        return data.copy()
    if win % 2 == 0:
        win += 1
    if data.size < win:
        return data.copy()
    kernel = np.ones(win, dtype=float) / float(win)
    pad = win // 2
    padded = np.pad(data, (pad, pad), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def save_simulation_plots(result: dict, out_dir: str | Path) -> list[Path]:
    """Сохраняет набор PNG по результату run_simulation."""
    import matplotlib.pyplot as plt

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for old_png in out.glob("*.png"):
        old_png.unlink()
    saved: list[Path] = []

    t = np.asarray(result["time"], dtype=float)
    x_true = np.asarray(result["x_true"], dtype=float)
    delta_e_ref = np.asarray(result["delta_e_ref"], dtype=float)
    delta_e_cmd = np.asarray(result["delta_e_cmd"], dtype=float)
    truth_nav = np.asarray(result["truth_nav"], dtype=float)
    bins_hist = np.asarray(result["bins"], dtype=float)
    gnss_hist = np.asarray(result["gnss"], dtype=float)
    altimeter_hist = np.asarray(result["altimeter"], dtype=float)
    theta_ins_hist = np.asarray(result["theta_ins"], dtype=float)
    q_gyro_hist = np.asarray(result["q_gyro"], dtype=float)
    ekf = result["ekf"]
    x_hat = np.asarray(ekf["x_hat"], dtype=float)
    p_hist = np.asarray(ekf["p"], dtype=float)
    c_bn_hist = np.asarray(result["c_bn"], dtype=float)
    p_true = np.asarray(result["p_true"], dtype=float)
    p_init = np.asarray(result["param_metrics"]["p_initial"], dtype=float)
    p_final = np.asarray(result["param_metrics"]["p_hat_final"], dtype=float)
    rel_init = np.asarray(result["param_metrics"]["rel_error_pct_initial"], dtype=float)
    rel_final = np.asarray(result["param_metrics"]["rel_error_pct_final"], dtype=float)
    pe = np.asarray(result["param_ekf"]["x_hat"], dtype=float)
    nav_profile = result.get("nav_profile", NavigationAccuracyProfile())

    cfg_dt = float(t[1] - t[0]) if len(t) > 1 else 0.01
    truth_state5 = np.column_stack([x_true[:, 0], x_true[:, 1], x_true[:, 3], x_true[:, 2], truth_nav[:, 5]])
    sensor_state5 = state5_from_sensor_measurements(
        gnss_hist, bins_hist, theta_ins_hist, q_gyro_hist, altimeter_hist, cfg_dt
    )
    ekf_state5 = state5_from_nav_solution(x_hat, bins_hist, c_bn_hist, cfg_dt)
    state5_delta = state_error_series(truth_state5, ekf_state5)
    state5_sig3 = sigma3_state5_from_nav_covariance(x_hat, bins_hist, p_hist, cfg_dt, profile=nav_profile)
    state_specs = [
        ("V", "m/s", 1.0, "state_v.png"),
        ("α", "deg", 180.0 / np.pi, "state_alpha.png"),
        ("θ", "deg", 180.0 / np.pi, "state_theta.png"),
        ("q", "deg/s", 180.0 / np.pi, "state_q.png"),
        ("h", "m", 1.0, "state_h.png"),
    ]
    t_mask = t <= min(40.0, float(t[-1])) if len(t) else np.array([], dtype=bool)

    for i, (name, unit, scale, filename) in enumerate(state_specs):
        style = state_plot_style(name)
        truth_plot = truth_state5[t_mask, i] * scale
        sensor_plot = sensor_state5[t_mask, i] * scale
        ekf_plot = ekf_state5[t_mask, i] * scale
        err_plot = smooth_series_for_display(
            state5_delta[t_mask, i] * scale,
            int(style["error_smoothing_window"]),
        )
        sig3_plot = state5_sig3[t_mask, i] * scale

        fig, axs = plt.subplots(1, 2, figsize=(11, 4), sharex=True)
        if style["model_on_top"]:
            axs[0].plot(t[t_mask], ekf_plot, color="r", lw=1.0, label="Оценка", zorder=2)
            axs[0].plot(t[t_mask], truth_plot, color="k", lw=1.1, label="Модель", zorder=3)
        else:
            axs[0].plot(t[t_mask], truth_plot, color="k", lw=1.0, label="Модель", zorder=2)
            axs[0].plot(t[t_mask], ekf_plot, color="r", lw=1.0, label="Оценка", zorder=3)
        axs[0].set_ylabel(f"{name} ({unit})")
        axs[0].set_xlabel("Время (с)")
        axs[0].grid(True, alpha=0.3)
        axs[0].legend(fontsize=8)

        if style["sigma_on_top"]:
            axs[1].plot(t[t_mask], err_plot, color="r", lw=1.0, label="Ошибка оценки", zorder=2)
            axs[1].fill_between(
                t[t_mask],
                -sig3_plot,
                sig3_plot,
                color="C0",
                alpha=float(style["sigma_alpha"]),
                label="±3σ",
                zorder=3,
            )
        else:
            axs[1].fill_between(
                t[t_mask],
                -sig3_plot,
                sig3_plot,
                color="C0",
                alpha=float(style["sigma_alpha"]),
                label="±3σ",
                zorder=1,
            )
            axs[1].plot(t[t_mask], err_plot, color="r", lw=1.0, label="Ошибка оценки", zorder=2)
        axs[1].axhline(0.0, color="k", lw=0.4)
        axs[1].set_ylabel(f"Δ{name} ({unit})")
        axs[1].set_xlabel("Время (с)")
        axs[1].grid(True, alpha=0.3)
        axs[1].legend(fontsize=8)

        fig.suptitle(name)
        fig.tight_layout()
        p = out / filename
        fig.savefig(p, dpi=140)
        plt.close(fig)
        saved.append(p)

    fig, ax = plt.subplots(figsize=(10, 4), sharex=True)
    ax.step(t, np.rad2deg(delta_e_ref), where="post", color="C2", lw=1.0)
    ax.axhline(0.0, color="k", lw=0.4)
    ax.set_xlabel("Время (с)")
    ax.set_ylabel("δe, deg")
    ax.set_title("Ступенчатая команда руля высоты")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    p = out / "elevator_deflection.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    saved.append(p)

    p_names = ["C_Lα", "C_Lδe", "C_mα", "C_mq", "C_mδe"]
    fig, axs = plt.subplots(5, 1, figsize=(10, 11), sharex=True)
    for j in range(5):
        trace = pe[:, 4 + j]
        axs[j].plot(t, trace, color="r", lw=1.1, label="Оценка", zorder=3)
        axs[j].axhline(p_true[j], color="k", ls="--", lw=0.9, label="Истина")
        axs[j].axhline(p_init[j], color="C1", ls=":", lw=0.9, label="Начальное")
        axs[j].plot(t[0], p_init[j], marker="o", color="C1", ms=4)
        axs[j].plot(t[-1], p_final[j], marker="o", color="C0", ms=4)
        axs[j].set_ylabel(p_names[j])
        y_all = np.array([trace.min(), trace.max(), p_true[j], p_init[j], p_final[j]], dtype=float)
        y_span = float(max(y_all.max() - y_all.min(), 1e-6))
        y_pad = 0.12 * y_span
        axs[j].set_ylim(float(y_all.min() - y_pad), float(y_all.max() + y_pad))
        axs[j].grid(True, alpha=0.3)
        axs[j].text(
            0.01,
            0.08,
            f"нач. ошибка={rel_init[j]:.1f}% | итог={rel_final[j]:.1f}%",
            transform=axs[j].transAxes,
            fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.7, "edgecolor": "none"},
        )
    axs[0].legend(fontsize=8)
    axs[-1].set_xlabel("Время (с)")
    fig.suptitle("Идентификация аэродинамических параметров: старт ±60% и итог относительно истины")
    fig.tight_layout()
    p = out / "aero_param_identification.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    saved.append(p)

    return saved


def scenario_slug(name: str) -> str:
    slug = re.sub(r"[^a-zA-Z0-9_-]+", "_", name.strip().lower()).strip("_")
    return slug or "scenario"


def write_scenario_results_md(result: dict, out_dir: str | Path) -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    profile = result.get("nav_profile", NavigationAccuracyProfile())
    ident_cfg = result.get("ident_window", IdentificationWindow())
    metrics = result["param_metrics"]
    p_true = np.asarray(result["p_true"], dtype=float)
    p_final = np.asarray(metrics["p_hat_final"], dtype=float)
    rel_final = np.asarray(metrics["rel_error_pct_final"], dtype=float)
    param_names = ["CLα", "CLδe", "Cmα", "Cmq", "Cmδe"]
    lines = [
        f"# Сценарий `{profile.name}`",
        "",
        "## Параметры точности навигации",
        f"- GNSS положение: x{profile.gnss_pos_scale:.2f}",
        f"- GNSS скорость: x{profile.gnss_vel_scale:.2f}",
        f"- Высотомер: x{profile.altimeter_scale:.2f}",
        f"- `θ_ins`: x{profile.theta_scale:.2f}",
        f"- `q_gyro`: x{profile.q_scale:.2f}",
        "",
        "## Окно идентификации",
        f"- начало: {ident_cfg.t_start:.2f} с",
        f"- конец: {ident_cfg.t_end:.2f} с",
        "",
        "## Итог параметров",
    ]
    for name, true_v, est_v, rel_v in zip(param_names, p_true, p_final, rel_final, strict=True):
        lines.append(f"- `{name}`: истина = {true_v:.6g}, оценка = {est_v:.6g}, ошибка = {rel_v:.3f}%")
    lines.extend(
        [
            "",
            "## Файлы графиков",
            "- `state_v.png`",
            "- `state_alpha.png`",
            "- `state_theta.png`",
            "- `state_q.png`",
            "- `state_h.png`",
            "- `elevator_deflection.png`",
            "- `aero_param_identification.png`",
            "",
        ]
    )
    path = out / "RESULTS.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def run_navigation_accuracy_sweep(
    out_root: str | Path,
    *,
    t_model: float = 30.0,
    scenarios: list[NavigationAccuracyProfile] | None = None,
    ident_window: IdentificationWindow | None = None,
) -> list[Path]:
    root = Path(out_root)
    root.mkdir(parents=True, exist_ok=True)
    scenario_list = scenarios or [
        NavigationAccuracyProfile(name="good_nav", gnss_pos_scale=0.5, gnss_vel_scale=0.5, altimeter_scale=0.5, theta_scale=0.7, q_scale=0.7),
        NavigationAccuracyProfile(name="baseline_nav"),
        NavigationAccuracyProfile(name="poor_nav", gnss_pos_scale=2.0, gnss_vel_scale=2.0, altimeter_scale=2.0, theta_scale=1.5, q_scale=1.5),
    ]
    ident_cfg = ident_window or IdentificationWindow()
    summary_lines = [
        "# Sweep по точности навигации",
        "",
        "| Сценарий | GNSS pos | GNSS vel | Alt | Theta | q | Cmq err % | RESULTS |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    written: list[Path] = []
    for profile in scenario_list:
        scenario_dir = root / scenario_slug(profile.name)
        result = run_simulation(t_model=t_model, profile=profile, ident_window=ident_cfg)
        save_simulation_plots(result, scenario_dir)
        md_path = write_scenario_results_md(result, scenario_dir)
        written.append(md_path)
        cmq_err = float(np.asarray(result["param_metrics"]["rel_error_pct_final"], dtype=float)[3])
        summary_lines.append(
            f"| `{profile.name}` | {profile.gnss_pos_scale:.2f} | {profile.gnss_vel_scale:.2f} | {profile.altimeter_scale:.2f} | {profile.theta_scale:.2f} | {profile.q_scale:.2f} | {cmq_err:.3f} | `{scenario_dir.name}/RESULTS.md` |"
        )
    summary_path = root / "SUMMARY.md"
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    written.append(summary_path)
    return written


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Симуляция БИНС/GNSS + нав. EKF + parameter EKF")
    p.add_argument("--t-model", type=float, default=120.0, help="Длительность модели, с")
    p.add_argument("--auto-trim-thrust", action="store_true", help="Режим с автотримом тяги")
    p.add_argument("--thrust", type=float, default=None, help="Фиксированная тяга, Н (иначе равновесие)")
    p.add_argument("--without-nz", action="store_true", help="Отключить n_z в измерениях parameter EKF")
    p.add_argument("--run-sweep", action="store_true", help="Запустить sweep по точности навигации")
    p.add_argument("--ident-start", type=float, default=5.0, help="Начало окна идентификации, с")
    p.add_argument("--ident-end", type=float, default=30.0, help="Конец окна идентификации, с")
    p.add_argument(
        "--out-dir",
        type=str,
        default="sim_plots",
        help="Каталог для PNG",
    )
    p.add_argument("--no-plots", action="store_true", help="Не сохранять графики")
    return p.parse_args()


if __name__ == "__main__":
    if __package__ is None:
        sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    args = _parse_args()
    thrust = args.thrust
    ident_cfg = IdentificationWindow(t_start=float(args.ident_start), t_end=float(args.ident_end))
    if bool(args.run_sweep):
        written = run_navigation_accuracy_sweep(args.out_dir, t_model=args.t_model, ident_window=ident_cfg)
        print(f"Сохранено {len(written)} markdown-файлов в {args.out_dir}")
        raise SystemExit(0)
    res = run_simulation(
        t_model=args.t_model,
        auto_trim_thrust=bool(args.auto_trim_thrust),
        thrust_newtons=thrust,
        param_ekf_with_nz=not bool(args.without_nz),
        ident_window=ident_cfg,
    )
    if not args.no_plots:
        paths = save_simulation_plots(res, args.out_dir)
        print(f"Сохранено {len(paths)} файлов в {args.out_dir}")
    print("param metrics RMSE:", res["param_metrics"]["rmse"])
