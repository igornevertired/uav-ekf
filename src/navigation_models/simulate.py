"""Ядро: симуляция, нав. EKF, кинематика fusion — импортируется из simulate.py."""

from __future__ import annotations

import argparse
import sys
from dataclasses import replace
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
from src.navigation_models.bins import bins2_step
from src.navigation_models.earth import earth_model
from src.navigation_models.gnss import GnssNoise, gnss_measurement
from src.navigation_models.noise import NoiseStreams
from src.navigation_models.zala421_16e5 import Zala421Config, zala421_default_config


def pitch_from_c_bn(c_bn: np.ndarray) -> float:
    """Тангаж θ из DCM тело→НСК [N,U,E]."""
    return float(np.arctan2(c_bn[1, 0], c_bn[0, 0]))


def _delta_angle(a0: float, a1: float) -> float:
    d = float(a1 - a0)
    return float((d + np.pi) % (2.0 * np.pi) - np.pi)


def ekf_fuse_bins_gnss(
    bins_hist: np.ndarray,
    gnss_hist: np.ndarray,
    truth_hist: np.ndarray,
    c_bn_hist: np.ndarray,
    dt: float,
) -> dict[str, np.ndarray]:
    n = bins_hist.shape[0]
    state_dim = 6
    meas_dim = 4
    x_hat = np.zeros((n, state_dim), dtype=float)
    p_hist = np.zeros((n, state_dim, state_dim), dtype=float)
    delta_truth_minus_hat = np.zeros((n, state_dim), dtype=float)
    sig3_hist = np.zeros((n, state_dim), dtype=float)

    gnss_pos_std = (1.0 / 111_000.0) * np.pi / 180.0
    gnss_vel_std = 1e-2

    bins_nav = np.column_stack(
        [bins_hist[:, 4], bins_hist[:, 5], bins_hist[:, 0], bins_hist[:, 2]]
    )
    vh_bins = bins_hist[:, 1].astype(float)
    theta_bins = np.array([pitch_from_c_bn(c_bn_hist[k]) for k in range(n)], dtype=float)

    H = np.zeros((meas_dim, state_dim), dtype=float)
    for i in range(meas_dim):
        H[i, i] = 1.0

    p0_vh = 0.5**2
    p0_th = (0.05) ** 2
    x = np.concatenate([gnss_hist[0].copy(), np.array([vh_bins[0], theta_bins[0]], dtype=float)])
    p = np.diag(
        [
            gnss_pos_std**2,
            gnss_pos_std**2,
            gnss_vel_std**2,
            gnss_vel_std**2,
            p0_vh,
            p0_th,
        ]
    )

    q_vel = 1e-5
    q_pos = 1e-14
    q_vh = 2e-6
    q_th = 3e-10
    q_mat = np.diag([q_pos, q_pos, q_vel, q_vel, q_vh, q_th])

    r_mat = np.diag([gnss_pos_std**2, gnss_pos_std**2, gnss_vel_std**2, gnss_vel_std**2])

    i6 = np.eye(state_dim, dtype=float)

    for k in range(n):
        if k == 0:
            x_pred = x.copy()
        else:
            dphi, dlm, dvn, dve = (bins_nav[k] - bins_nav[k - 1]).tolist()
            dvh = float(vh_bins[k] - vh_bins[k - 1])
            dth = _delta_angle(theta_bins[k - 1], theta_bins[k])
            x_pred = x + np.array([dphi, dlm, dvn, dve, dvh, dth], dtype=float)

        p_pred = p + q_mat

        innov = gnss_hist[k] - x_pred[:meas_dim]
        s = H @ p_pred @ H.T + r_mat
        k_gain = p_pred @ H.T @ np.linalg.inv(s)
        x = x_pred + k_gain @ innov
        p = (i6 - k_gain @ H) @ p_pred

        x_hat[k] = x
        p_hist[k] = p
        delta_truth_minus_hat[k] = truth_hist[k] - x
        sig3_hist[k] = 3.0 * np.sqrt(np.maximum(np.diag(p), 0.0))

    return {
        "x_hat": x_hat,
        "p": p_hist,
        "delta_truth_minus_hat": delta_truth_minus_hat,
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
    if x_hat_nav.shape[1] >= 6:
        vh = x_hat_nav[:, 4]
        theta_hat = x_hat_nav[:, 5]
    else:
        vh = bins_hist[:, 1]
        theta_hat = np.array([pitch_from_c_bn(c_bn_hist[k]) for k in range(n)], dtype=float)
    r_h = np.sqrt(np.maximum(vn**2 + ve**2, 1e-12))
    v_mag = np.sqrt(np.maximum(vn**2 + ve**2 + vh**2, 1e-12))
    gamma_hat = np.arctan2(vh, r_h)
    alpha_hat = theta_hat - gamma_hat
    q_hat = np.gradient(theta_hat, dt)
    return np.column_stack([v_mag, alpha_hat, q_hat, theta_hat])


upper_stage_kinematics_for_param_ekf = longitudinal_state_from_fusion


def sigma3_longitudinal_state(
    x_hat_nav: np.ndarray,
    bins_hist: np.ndarray,
    p_hist: np.ndarray,
    gnss_vel_std: float,
    sigma_theta_ins: float,
    dt: float,
    *,
    sigma_vh: float | None = None,
    inflation: float = 1.08,
) -> np.ndarray:
    n = x_hat_nav.shape[0]
    sig3 = np.zeros((n, 4), dtype=float)
    sig_vh = float(gnss_vel_std if sigma_vh is None else sigma_vh)
    sigma_q_gyro = 0.02
    nav6 = x_hat_nav.shape[1] >= 6 and p_hist.shape[1] >= 6

    for k in range(n):
        vn, ve = float(x_hat_nav[k, 2]), float(x_hat_nav[k, 3])
        vh = float(x_hat_nav[k, 4]) if nav6 else float(bins_hist[k, 1])
        p = p_hist[k]
        p22, p33 = float(p[2, 2]), float(p[3, 3])
        if nav6:
            p44, p55 = float(p[4, 4]), float(p[5, 5])
            sig_th_k = np.sqrt(max(p55, 0.0))
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
        var_alpha = p55 + var_gamma if nav6 else sigma_theta_ins**2 + var_gamma
        sig3[k, 1] = 3.0 * np.sqrt(max(var_alpha, 0.0))

        sig3[k, 3] = 3.0 * sig_th_k
        sig_q_k = np.sqrt((np.sqrt(2.0) * sig_th_k / max(dt, 1e-9)) ** 2 + sigma_q_gyro**2)
        sig3[k, 2] = 3.0 * sig_q_k

    sig3 *= float(inflation)
    return sig3

def run_simulation(
    config: Zala421Config | None = None,
    t_model: float = 120.0,
    *,
    auto_trim_thrust: bool = False,
    thrust_newtons: float | None = None,
    param_ekf_with_nz: bool = False,
) -> dict[str, np.ndarray | dict]:
    cfg = config or zala421_default_config()
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
    t_hist = np.arange(steps, dtype=float) * cfg.dt
    x_true_hist = np.zeros((steps, 4), dtype=float)
    delta_e_cmd_hist = np.zeros(steps, dtype=float)
    truth_nav_hist = np.zeros((steps, 6), dtype=float)

    fi_true = float(cfg.np0[4])
    lm_true = float(cfg.np0[5])
    h_true = float(cfg.np0[3])

    n_pulses = 4
    pulse_width = 3.0
    pulse_centers = np.linspace(
        0.5 * pulse_width + 3.0,
        max(t_model - 0.5 * pulse_width - 3.0, 0.5 * pulse_width + 3.0),
        n_pulses,
    )
    pulse_amp_deg = np.array([-1.25, 1.25, -1.25, 1.25], dtype=float)

    for i in range(steps):
        t = t_hist[i]
        delta_e_cmd_deg = 0.0
        for center, amp in zip(pulse_centers, pulse_amp_deg):
            if abs(t - center) <= pulse_width / 2.0:
                delta_e_cmd_deg = float(amp)
                break
        delta_e = np.deg2rad(delta_e_cmd_deg)
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
            [fi_true, lm_true, vn_true, ve_true, vh_true, theta_true],
            dtype=float,
        )

        w_true_body = np.array([0.0, q_pitch, 0.0], dtype=float)
        imu_noise = noise_streams.sample_imu()
        q, c_bn, np_state, _ = bins2_step(
            a_true=a_true_body,
            w_true=w_true_body,
            np_state=np_state,
            c_bn=c_bn,
            q=q,
            dt=cfg.dt,
            h_et=h_true,
            dh_et=vh_true,
            params=cfg.bins_params,
            noise=imu_noise,
        )
        bins_hist[i] = np_state
        c_bn_hist[i] = c_bn

        wn_fi, wn_lm, wn_vn, wn_ve = noise_streams.sample_gnss()
        z, _ = gnss_measurement(
            fi=fi_true,
            lm=lm_true,
            vn=vn_true,
            ve=ve_true,
            noise=GnssNoise(wn_fi=wn_fi, wn_lm=wn_lm, wn_vn=wn_vn, wn_ve=wn_ve),
            scheme=cfg.gnss_scheme,
        )
        gnss_hist[i] = z

    ekf_out = ekf_fuse_bins_gnss(
        bins_hist=bins_hist,
        gnss_hist=gnss_hist,
        truth_hist=truth_nav_hist,
        c_bn_hist=c_bn_hist,
        dt=cfg.dt,
    )

    pe_hist = np.zeros((steps, 9), dtype=float)
    meas_dim = 5 if param_ekf_with_nz else 4
    y_pe_hist = np.zeros((steps, meas_dim), dtype=float)
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
    x_pe[:4] = x_true_hist[0].copy()
    p_init = p_true * np.array([0.92, 1.08, 0.94, 1.04, 0.96], dtype=float)
    x_pe[4:9] = p_init.copy()
    p_pe = np.diag(
        [
            0.5**2,
            0.02**2,
            0.02**2,
            0.02**2,
            0.6**2,
            0.12**2,
            0.15**2,
            2.0**2,
            0.15**2,
        ]
    )
    r_list = [
        (5e-3) ** 2,
        (5e-5) ** 2,
        (5e-4) ** 2,
        (5e-5) ** 2,
    ]
    if param_ekf_with_nz:
        r_list.append((5e-3) ** 2)
    pe_cfg = ParameterEkfConfig(
        q_diag=np.array(
            [8e-5, 2e-6, 2e-6, 2e-6, 0.0, 0.0, 0.0, 0.0, 0.0],
            dtype=float,
        ),
        r_diag=np.array(r_list, dtype=float),
        include_nz=bool(param_ekf_with_nz),
    )
    x_hat_nav = ekf_out["x_hat"]
    long_nav = upper_stage_kinematics_for_param_ekf(x_hat_nav, bins_hist, c_bn_hist, cfg.dt)
    for k in range(steps):
        de = float(delta_e_cmd_hist[k])
        yk = np.asarray(long_nav[k], dtype=float).copy()
        if param_ekf_with_nz:
            nz_m = load_factor_nz(long_nav[k], de, plant_params)
            yk = np.concatenate([yk, np.array([nz_m], dtype=float)])
        y_pe_hist[k] = yk
        x_pe, p_pe = parameter_ekf_step(x_pe, p_pe, yk, de, cfg.dt, plant_params, pe_cfg)
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
        "delta_e_cmd": delta_e_cmd_hist,
        "truth_nav": truth_nav_hist,
        "bins": bins_hist,
        "c_bn": c_bn_hist,
        "gnss": gnss_hist,
        "ekf": ekf_out,
        "plant_params": plant_params,
        "p_true": p_true,
        "param_ekf": {"x_hat": pe_hist, "y": y_pe_hist, "p_cov": p_pe},
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


def save_simulation_plots(result: dict, out_dir: str | Path) -> list[Path]:
    """Сохраняет набор PNG по результату run_simulation."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    saved: list[Path] = []

    t = np.asarray(result["time"], dtype=float)
    x_true = np.asarray(result["x_true"], dtype=float)
    truth_nav = np.asarray(result["truth_nav"], dtype=float)
    bins_hist = np.asarray(result["bins"], dtype=float)
    gnss_hist = np.asarray(result["gnss"], dtype=float)
    ekf = result["ekf"]
    x_hat = np.asarray(ekf["x_hat"], dtype=float)
    delta = np.asarray(ekf["delta_truth_minus_hat"], dtype=float)
    sig3_nav = np.asarray(ekf["sig3"], dtype=float)
    c_bn_hist = np.asarray(result["c_bn"], dtype=float)
    cfg_dt = float(t[1] - t[0]) if len(t) > 1 else 0.01

    labels_nav = ["φ", "λ", "Vn", "Ve", "Vh", "θ"]
    fig, axs = plt.subplots(3, 2, figsize=(11, 8), sharex=True)
    axs = axs.ravel()
    for i in range(6):
        axs[i].plot(t, truth_nav[:, i], "k-", lw=1.0, label="истина")
        if i < 4:
            axs[i].plot(t, gnss_hist[:, i], "--", lw=0.8, alpha=0.85, label="GNSS")
        if i == 0:
            axs[i].plot(t, bins_hist[:, 4], ":", lw=0.9, alpha=0.9, label="БИНС φ")
        elif i == 1:
            axs[i].plot(t, bins_hist[:, 5], ":", lw=0.9, alpha=0.9, label="БИНС λ")
        elif i == 2:
            axs[i].plot(t, bins_hist[:, 0], ":", lw=0.9, alpha=0.9, label="БИНС Vn")
        elif i == 3:
            axs[i].plot(t, bins_hist[:, 2], ":", lw=0.9, alpha=0.9, label="БИНС Ve")
        elif i == 4:
            axs[i].plot(t, bins_hist[:, 1], ":", lw=0.9, alpha=0.9, label="БИНС Vh")
        elif i == 5:
            th_bins = np.array([pitch_from_c_bn(c_bn_hist[k]) for k in range(len(t))])
            axs[i].plot(t, th_bins, ":", lw=0.9, alpha=0.9, label="БИНС θ")
        axs[i].set_ylabel(labels_nav[i])
        axs[i].grid(True, alpha=0.3)
    axs[0].legend(fontsize=7, loc="upper right")
    axs[-2].set_xlabel("t, с")
    axs[-1].set_xlabel("t, с")
    fig.suptitle("Навигация: истина, GNSS, БИНС")
    fig.tight_layout()
    p = out / "nav_sensors.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    saved.append(p)

    fig, axs = plt.subplots(3, 2, figsize=(11, 8), sharex=True)
    axs = axs.ravel()
    for i in range(6):
        axs[i].plot(t, truth_nav[:, i], "k-", lw=1.0, label="истина")
        _fill_sigma_bands_gradient(axs[i], t, x_hat[:, i], sig3_nav[:, i], color="C0")
        axs[i].set_ylabel(labels_nav[i])
        axs[i].grid(True, alpha=0.3)
    h = [
        Line2D([0], [0], color="k", lw=1.0, label="истина"),
        Line2D([0], [0], color="C0", lw=1.0, label="EKF"),
        Patch(facecolor="C0", edgecolor="none", alpha=0.15, label="±3σ"),
    ]
    axs[0].legend(handles=h, fontsize=8)
    axs[-2].set_xlabel("t, с")
    axs[-1].set_xlabel("t, с")
    fig.suptitle("Навигационный EKF: оценка и ±3σ")
    fig.tight_layout()
    p = out / "ekf_nav.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    saved.append(p)

    fig, axs = plt.subplots(3, 2, figsize=(11, 8), sharex=True)
    axs = axs.ravel()
    for i in range(6):
        axs[i].plot(t, delta[:, i], color="C3", lw=0.9)
        axs[i].fill_between(t, -sig3_nav[:, i], sig3_nav[:, i], color="C0", alpha=0.15)
        axs[i].axhline(0.0, color="k", lw=0.4)
        axs[i].set_ylabel(f"Δ {labels_nav[i]}")
        axs[i].grid(True, alpha=0.3)
    axs[-2].set_xlabel("t, с")
    axs[-1].set_xlabel("t, с")
    fig.suptitle("Ошибка навигационного EKF (истина − оценка), полоса ±3σ")
    fig.tight_layout()
    p = out / "ekf_nav_error.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    saved.append(p)

    fig, axs = plt.subplots(2, 2, figsize=(10, 6), sharex=True)
    names = ["V", "α", "q", "θ"]
    for j in range(4):
        ax = axs.ravel()[j]
        ax.plot(t, x_true[:, j], "k-", lw=1.0, label="истина")
        ax.set_ylabel(names[j])
        ax.grid(True, alpha=0.3)
    axs[0, 0].legend(fontsize=8)
    fig.suptitle("Продольное движение (истина)")
    fig.tight_layout()
    p = out / "longitudinal_truth.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    saved.append(p)

    long_nav = longitudinal_state_from_fusion(x_hat, bins_hist, c_bn_hist, cfg_dt)
    long_sig3 = sigma3_longitudinal_state(
        x_hat,
        bins_hist,
        ekf["p"],
        gnss_vel_std=1e-2,
        sigma_theta_ins=0.05,
        dt=cfg_dt,
    )
    err_long = long_nav - x_true
    fig, axs = plt.subplots(2, 2, figsize=(10, 6), sharex=True)
    for j in range(4):
        ax = axs.ravel()[j]
        ax.plot(t, long_nav[:, j], color="C0", lw=0.9, label="из нав. EKF")
        ax.plot(t, x_true[:, j], "k--", lw=0.8, alpha=0.85, label="истина")
        ax.set_ylabel(names[j])
        ax.grid(True, alpha=0.3)
    axs[0, 0].legend(fontsize=7)
    fig.suptitle("Продольное состояние: кинематика из нав. EKF vs истина")
    fig.tight_layout()
    p = out / "longitudinal_nav_vs_truth.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    saved.append(p)

    fig, axs = plt.subplots(2, 2, figsize=(10, 6), sharex=True)
    for j in range(4):
        ax = axs.ravel()[j]
        ax.plot(t, err_long[:, j], color="C3", lw=0.8)
        ax.fill_between(t, -long_sig3[:, j], long_sig3[:, j], color="C0", alpha=0.12)
        ax.axhline(0.0, color="k", lw=0.4)
        ax.set_ylabel(f"Δ {names[j]}")
        ax.grid(True, alpha=0.3)
    fig.suptitle("Ошибка продольной кинематики ±3σ (аппрокс.)")
    fig.tight_layout()
    p = out / "longitudinal_error_sig3.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    saved.append(p)

    pe = np.asarray(result["param_ekf"]["x_hat"], dtype=float)
    p_true = np.asarray(result["p_true"], dtype=float)
    p_names = ["C_Lα", "C_Lδe", "C_mα", "C_mq", "C_mδe"]
    fig, axs = plt.subplots(5, 1, figsize=(10, 10), sharex=True)
    for j in range(5):
        axs[j].plot(t, pe[:, 4 + j], color="C0", lw=0.9, label="оценка")
        axs[j].axhline(p_true[j], color="k", ls="--", lw=0.9, label="истина")
        axs[j].set_ylabel(p_names[j])
        axs[j].grid(True, alpha=0.3)
    axs[0].legend(fontsize=8)
    axs[-1].set_xlabel("t, с")
    fig.suptitle("Parameter EKF: аэродинамические коэффициенты")
    fig.tight_layout()
    p = out / "param_ekf.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    saved.append(p)

    param_err = np.asarray(result["param_error"], dtype=float)
    fig, axs = plt.subplots(5, 1, figsize=(10, 10), sharex=True)
    for j in range(5):
        axs[j].plot(t, param_err[:, j], color="C3", lw=0.8)
        axs[j].axhline(0.0, color="k", lw=0.4)
        axs[j].set_ylabel(f"Δ {p_names[j]}")
        axs[j].grid(True, alpha=0.3)
    axs[-1].set_xlabel("t, с")
    fig.suptitle("Ошибка идентификации параметров")
    fig.tight_layout()
    p = out / "param_ekf_error.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    saved.append(p)

    return saved


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Симуляция БИНС/GNSS + нав. EKF + parameter EKF")
    p.add_argument("--t-model", type=float, default=120.0, help="Длительность модели, с")
    p.add_argument("--auto-trim-thrust", action="store_true", help="Режим с автотримом тяги")
    p.add_argument("--thrust", type=float, default=None, help="Фиксированная тяга, Н (иначе равновесие)")
    p.add_argument("--with-nz", action="store_true", help="Включить n_z в измерения parameter EKF")
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
    res = run_simulation(
        t_model=args.t_model,
        auto_trim_thrust=bool(args.auto_trim_thrust),
        thrust_newtons=thrust,
        param_ekf_with_nz=bool(args.with_nz),
    )
    if not args.no_plots:
        paths = save_simulation_plots(res, args.out_dir)
        print(f"Сохранено {len(paths)} файлов в {args.out_dir}")
    print("param metrics RMSE:", res["param_metrics"]["rmse"])
