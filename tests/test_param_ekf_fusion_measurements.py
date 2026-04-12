"""Измерения Parameter EKF: кинематика из навигации."""

from pathlib import Path

import numpy as np
from dataclasses import replace

from models.aircraft_dynamics import NonlinearLongitudinalParams
from src.navigation_models.altimeter import altimeter_reference_sigma
from src.navigation_models.bins import Bins2Params, bins2_step
from src.navigation_models.gnss import gnss_reference_sigmas
from src.navigation_models.noise import (
    ImuNoise,
    imu_reference_sigmas,
    q_gyro_measurement_sigma,
    theta_ins_measurement_sigma,
)
from src.navigation_models.simulate import (
    elevator_step_command_deg,
    longitudinal_measurement_covariance_from_nav,
    longitudinal_measurement_from_nav_state,
    navigation_outputs_to_parameter_measurements,
    q_excitation_command_deg,
    run_simulation,
    longitudinal_state_from_fusion,
    nav_ekf_initial_sigmas,
    nominal_aero_params_for_nz,
    pitch_from_c_bn,
    save_simulation_plots,
    sigma3_state5_from_nav_covariance,
    smooth_series_for_display,
    state_error_series,
    state_plot_style,
    state5_from_nav_solution,
)
from src.navigation_models.zala421_16e5 import initial_attitude_nue, zala421_default_config


def test_nominal_aero_keeps_mass_geometry_changes_lift_model():
    base = NonlinearLongitudinalParams.default()
    plant = replace(base, cla=99.0, cl_de=0.99)
    nom = nominal_aero_params_for_nz(plant)
    assert nom.m == plant.m and nom.S == plant.S and nom.rho == plant.rho
    assert nom.cla == base.cla and nom.cl_de == base.cl_de
    assert nom.cla != plant.cla


def test_longitudinal_state_from_fusion_matches_manual_first_row():
    dt = 0.01
    n = 3
    x_hat = np.zeros((n, 4))
    x_hat[:, 2] = [3.0, 3.0, 3.0]
    x_hat[:, 3] = [0.1, 0.1, 0.1]
    bins = np.zeros((n, 6))
    bins[:, 1] = [0.2, 0.2, 0.2]
    c_bn = np.tile(np.eye(3), (n, 1, 1))
    out = longitudinal_state_from_fusion(x_hat, bins, c_bn, dt)
    vn, ve, vh = 3.0, 0.1, 0.2
    v_exp = float(np.sqrt(vn * vn + ve * ve + vh * vh))
    assert abs(out[0, 0] - v_exp) < 1e-9
    assert abs(out[0, 3] - pitch_from_c_bn(c_bn[0])) < 1e-9


def test_longitudinal_state_from_fusion_7d_uses_state_vh_theta():
    dt = 0.01
    n = 2
    th = 0.05
    h = 123.0
    x_hat = np.zeros((n, 7))
    x_hat[:, 2] = 3.0
    x_hat[:, 3] = 0.1
    x_hat[:, 4] = 0.2
    x_hat[:, 5] = h
    x_hat[:, 6] = th
    bins = np.zeros((n, 6))
    c_bn = np.tile(np.eye(3), (n, 1, 1))
    out = longitudinal_state_from_fusion(x_hat, bins, c_bn, dt)
    assert abs(out[0, 3] - th) < 1e-9


def test_altimeter_reference_sigma_is_one_meter():
    assert np.isclose(altimeter_reference_sigma(), 1.0)


def test_q_gyro_measurement_sigma_is_large_enough_for_visible_sensor_error():
    assert np.isclose(q_gyro_measurement_sigma(), 0.01)


def test_theta_ins_measurement_sigma_matches_visible_theta_error_level():
    assert np.isclose(theta_ins_measurement_sigma(), np.deg2rad(0.05))


def test_nav_ekf_initial_sigmas_match_reference_table():
    sig = nav_ekf_initial_sigmas()
    assert np.isclose(sig["phi"], 1.57e-7)
    assert np.isclose(sig["lambda"], 1.57e-7)
    assert np.isclose(sig["vn"], 0.02)
    assert np.isclose(sig["ve"], 0.02)
    assert np.isclose(sig["vh"], 0.02)
    assert np.isclose(sig["theta"], 1.0e-4)


def test_gnss_reference_sigmas_match_table_2():
    sig = gnss_reference_sigmas()
    pos_std_rad = 12.0 / 111_000.0 * np.pi / 180.0
    assert np.isclose(sig["phi"], pos_std_rad)
    assert np.isclose(sig["lambda"], pos_std_rad)
    assert np.isclose(sig["vn"], 0.5)
    assert np.isclose(sig["ve"], 0.5)
    assert np.isclose(sig["pseudorange"], 12.0)
    assert np.isclose(sig["pseudorange_rate"], 0.5)


def test_imu_reference_sigmas_match_table_3():
    sig = imu_reference_sigmas()
    gyro_bias_rad_s = np.deg2rad(0.003) / 3600.0
    gyro_rw_rad_sqrt_s = np.deg2rad(0.0015) / 3600.0
    accel_bias_m_s2 = 25e-6 * 9.80665
    accel_rw_m_s2_sqrt_hz = 50e-6 * 9.80665
    assert np.allclose(sig["gyro_bias"], np.full(3, gyro_bias_rad_s))
    assert np.allclose(sig["gyro_rw"], np.full(3, gyro_rw_rad_sqrt_s))
    assert np.allclose(sig["accel_bias"], np.full(3, accel_bias_m_s2))
    assert np.allclose(sig["accel_rw"], np.full(3, accel_rw_m_s2_sqrt_hz))


def test_zala_default_config_uses_reference_imu_biases():
    cfg = zala421_default_config()
    sig = imu_reference_sigmas()
    assert np.allclose(cfg.bins_params.dw, sig["gyro_bias"])
    assert np.allclose(cfg.bins_params.da, sig["accel_bias"])


def test_bins_vertical_channel_integrates_internal_state_without_truth_injection():
    params = Bins2Params(
        dw=np.zeros(3, dtype=float),
        da=np.zeros(3, dtype=float),
        kmw=np.zeros((3, 3), dtype=float),
        kma=np.zeros((3, 3), dtype=float),
        fiw=np.eye(3, dtype=float),
        fia=np.eye(3, dtype=float),
    )
    noise = ImuNoise(wn_w=np.zeros(3, dtype=float), wn_a=np.zeros(3, dtype=float))
    c_bn, q = initial_attitude_nue(0.0)
    np_state = np.array([10.0, 3.0, 0.0, 100.0, np.deg2rad(55.0), np.deg2rad(37.0)], dtype=float)

    _, _, np_new, _ = bins2_step(
        a_true=np.zeros(3, dtype=float),
        w_true=np.zeros(3, dtype=float),
        np_state=np_state,
        c_bn=c_bn,
        q=q,
        dt=0.2,
        params=params,
        noise=noise,
    )

    assert np.isclose(np_new[3], np_state[3] + np_state[1] * 0.2)


def test_run_simulation_returns_altimeter_and_7d_nav_state():
    res = run_simulation(t_model=0.05)
    assert np.asarray(res["truth_nav"]).shape[1] == 8
    assert np.asarray(res["altimeter"]).ndim == 1
    assert np.asarray(res["theta_ins"]).ndim == 1
    assert np.asarray(res["q_gyro"]).ndim == 1
    assert int(res["ekf"]["nav_state_dim"]) == 8


def test_longitudinal_measurement_covariance_is_built_from_nav_covariance():
    x_nav = np.array([1.0, 2.0, 30.0, 0.2, 1.5, 300.0, 0.08, 0.01], dtype=float)
    p_small = np.diag([1e-10, 1e-10, 0.01, 0.01, 0.01, 0.04, 1e-4, 4e-4])
    p_big = 9.0 * p_small
    r_small = longitudinal_measurement_covariance_from_nav(x_nav, p_small)
    r_big = longitudinal_measurement_covariance_from_nav(x_nav, p_big)
    assert r_small.shape == (4, 4)
    assert np.all(np.diag(r_big) > np.diag(r_small))


def test_navigation_outputs_to_parameter_measurements_include_qdot():
    plant = NonlinearLongitudinalParams.default()
    x_hat_nav = np.array(
        [
            [0.0, 0.0, 30.0, 0.0, 0.5, 300.0, 0.08, 0.00],
            [0.0, 0.0, 30.0, 0.0, 0.5, 300.0, 0.08, 0.03],
            [0.0, 0.0, 30.0, 0.0, 0.5, 300.0, 0.08, -0.01],
        ],
        dtype=float,
    )
    p_nav = np.tile(np.diag([1e-8, 1e-8, 0.05, 0.05, 0.05, 1.0, 1e-3, 2e-3]), (3, 1, 1))
    delta_e = np.array([0.0, 0.02, -0.02], dtype=float)
    y_hist, r_hist = navigation_outputs_to_parameter_measurements(
        x_hat_nav,
        p_nav,
        include_nz=True,
        include_qdot=True,
        delta_e_hist=delta_e,
        plant_params=plant,
        dt=0.01,
    )
    assert y_hist.shape == (3, 6)
    assert r_hist.shape == (3, 6, 6)
    assert abs(float(y_hist[1, -1])) > 0.1


def test_state_error_series_returns_ekf_minus_truth():
    truth = np.array([[1.0, 2.0], [3.0, 4.0]])
    est = np.array([[0.7, 2.5], [3.5, 3.0]])
    err = state_error_series(truth, est)
    assert np.allclose(err, np.array([[-0.3, 0.5], [0.5, -1.0]]))


def test_state_plot_style_highlights_theta_and_q():
    style_theta = state_plot_style("θ")
    style_q = state_plot_style("q")
    style_v = state_plot_style("V")
    assert style_theta["model_on_top"] is True
    assert style_q["sigma_on_top"] is True
    assert style_theta["sigma_alpha"] > style_v["sigma_alpha"]
    assert style_q["error_smoothing_window"] > style_v["error_smoothing_window"]
    assert style_v["model_on_top"] is False


def test_smooth_series_for_display_preserves_length_and_softens_spikes():
    src = np.array([0.0, 0.0, 9.0, 0.0, 0.0], dtype=float)
    smoothed = smooth_series_for_display(src, 3)
    assert smoothed.shape == src.shape
    assert float(smoothed[2]) < float(src[2])
    assert np.isclose(smoothed.sum(), src.sum())


def test_run_simulation_initializes_parameter_filter_from_nav_solution():
    res = run_simulation(t_model=0.05)
    x0 = np.asarray(res["param_ekf"]["x0"], dtype=float)
    nav0 = np.asarray(res["ekf"]["x_hat"], dtype=float)[0]
    y0 = longitudinal_measurement_from_nav_state(nav0)
    assert np.allclose(x0[:4], y0)
    assert np.asarray(res["param_ekf"]["r_meas"]).shape[1:] == (5, 5)


def test_run_simulation_enables_nz_in_parameter_filter_by_default():
    res = run_simulation(t_model=0.05)
    assert bool(res["param_ekf_with_nz"]) is True
    assert np.asarray(res["param_ekf"]["y"]).shape[1] == 5


def test_q_excitation_command_adds_extra_pitch_rate_excitation():
    assert np.isclose(q_excitation_command_deg(4.0), 0.0)
    assert np.isclose(q_excitation_command_deg(5.5), -2.5)
    assert np.isclose(q_excitation_command_deg(6.5), 2.5)
    assert np.isclose(q_excitation_command_deg(14.5), 2.5)
    assert np.isclose(q_excitation_command_deg(21.0), 0.0)


def test_run_simulation_produces_larger_q_excursions_for_identification():
    res = run_simulation(t_model=25.0)
    q_hist = np.asarray(res["x_true"], dtype=float)[:, 2]
    assert float(np.max(np.abs(q_hist))) > 0.05


def test_q_gyro_measurement_is_not_identical_to_truth():
    res = run_simulation(t_model=10.0)
    q_true = np.asarray(res["truth_nav"], dtype=float)[:, 7]
    q_meas = np.asarray(res["q_gyro"], dtype=float)
    assert float(np.std(q_meas - q_true)) > 0.005


def test_run_simulation_keeps_q_rmse_below_point_seven_deg_per_sec():
    res = run_simulation(t_model=30.0)
    truth_nav = np.asarray(res["truth_nav"], dtype=float)
    x_true = np.asarray(res["x_true"], dtype=float)
    truth_state5 = np.column_stack([x_true[:, 0], x_true[:, 1], x_true[:, 3], x_true[:, 2], truth_nav[:, 5]])
    ekf_state5 = state5_from_nav_solution(
        np.asarray(res["ekf"]["x_hat"], dtype=float),
        np.asarray(res["bins"], dtype=float),
        np.asarray(res["c_bn"], dtype=float),
        0.01,
    )
    q_err = ekf_state5[:, 3] - truth_state5[:, 3]
    q_rmse_deg = float(np.rad2deg(np.sqrt(np.mean(q_err**2))))
    assert q_rmse_deg < 0.7


def test_run_simulation_keeps_cmq_closer_to_truth():
    res = run_simulation(t_model=30.0)
    rel_final = np.asarray(res["param_metrics"]["rel_error_pct_final"], dtype=float)
    assert float(rel_final[3]) < 45.0


def test_alpha_error_is_mostly_inside_three_sigma():
    res = run_simulation(t_model=30.0)
    truth_nav = np.asarray(res["truth_nav"], dtype=float)
    x_true = np.asarray(res["x_true"], dtype=float)
    truth_state5 = np.column_stack([x_true[:, 0], x_true[:, 1], x_true[:, 3], x_true[:, 2], truth_nav[:, 5]])
    ekf_state5 = state5_from_nav_solution(
        np.asarray(res["ekf"]["x_hat"], dtype=float),
        np.asarray(res["bins"], dtype=float),
        np.asarray(res["c_bn"], dtype=float),
        0.01,
    )
    sig3 = sigma3_state5_from_nav_covariance(
        np.asarray(res["ekf"]["x_hat"], dtype=float),
        np.asarray(res["bins"], dtype=float),
        np.asarray(res["ekf"]["p"], dtype=float),
        0.01,
    )
    alpha_err = np.abs(ekf_state5[:, 1] - truth_state5[:, 1])
    inside = alpha_err <= sig3[:, 1]
    assert float(np.mean(inside)) > 0.9


def test_run_simulation_holds_altitude_reasonably():
    res = run_simulation(t_model=30.0)
    h = np.asarray(res["truth_nav"], dtype=float)[:, 5]
    assert abs(float(h[-1] - h[0])) < 15.0


def test_elevator_step_command_is_piecewise_constant_for_two_seconds():
    assert np.isclose(elevator_step_command_deg(0.0, 120.0), 0.0)
    assert np.isclose(elevator_step_command_deg(4.0, 120.0), 0.0)
    assert np.isclose(elevator_step_command_deg(6.0, 120.0), -1.0)
    assert np.isclose(elevator_step_command_deg(9.0, 120.0), 1.0)
    assert np.isclose(elevator_step_command_deg(12.0, 120.0), 0.0)
    assert np.isclose(elevator_step_command_deg(14.0, 120.0), -1.0)
    assert np.isclose(elevator_step_command_deg(16.0, 120.0), 1.0)
    assert np.isclose(elevator_step_command_deg(18.2, 120.0), -1.0)
    assert np.isclose(elevator_step_command_deg(19.4, 120.0), 1.0)
    assert np.isclose(elevator_step_command_deg(22.0, 120.0), 0.0)


def test_aero_param_initial_errors_are_sixty_percent():
    res = run_simulation(t_model=0.05)
    initial_pct = np.asarray(res["param_metrics"]["rel_error_pct_initial"], dtype=float)
    assert np.allclose(initial_pct, np.full(5, 60.0))


def test_save_simulation_plots_creates_only_requested_files(tmp_path):
    res = run_simulation(t_model=0.05)
    paths = save_simulation_plots(res, tmp_path)
    names = sorted(path.name for path in paths)
    assert names == [
        "aero_param_identification.png",
        "elevator_deflection.png",
        "state_alpha.png",
        "state_h.png",
        "state_q.png",
        "state_theta.png",
        "state_v.png",
    ]


def test_sensor_models_document_exists():
    doc = Path("SENSOR_MODELS.md")
    assert doc.exists()
    text = doc.read_text(encoding="utf-8")
    assert "БИНС" in text and "GNSS" in text
    assert "Первый EKF" in text
    assert "Второй EKF" in text
    assert "Что подается на первый EKF" in text
    assert "Что выдает первый EKF" in text
    assert "Что подается на второй EKF" in text
    assert "Что выдает второй EKF" in text
    assert "Как формируются графики" in text
    assert "aero_param_identification.png" in text
