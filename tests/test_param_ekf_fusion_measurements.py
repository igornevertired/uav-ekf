"""Измерения Parameter EKF: кинематика из навигации."""

import numpy as np
from dataclasses import replace

from models.aircraft_dynamics import NonlinearLongitudinalParams
from src.navigation_models.simulate import longitudinal_state_from_fusion, nominal_aero_params_for_nz, pitch_from_c_bn


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


def test_longitudinal_state_from_fusion_6d_uses_state_vh_theta():
    dt = 0.01
    n = 2
    th = 0.05
    x_hat = np.zeros((n, 6))
    x_hat[:, 2] = 3.0
    x_hat[:, 3] = 0.1
    x_hat[:, 4] = 0.2
    x_hat[:, 5] = th
    bins = np.zeros((n, 6))
    c_bn = np.tile(np.eye(3), (n, 1, 1))
    out = longitudinal_state_from_fusion(x_hat, bins, c_bn, dt)
    assert abs(out[0, 3] - th) < 1e-9
