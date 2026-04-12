from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .altimeter import altimeter_reference_sigma
from .gnss import gnss_reference_sigmas


@dataclass(frozen=True)
class ImuNoise:
    wn_w: np.ndarray
    wn_a: np.ndarray


def imu_reference_sigmas() -> dict[str, np.ndarray]:
    """Табличные 1σ параметры INS из Table 3."""
    gyro_bias = np.full(3, np.deg2rad(0.003) / 3600.0, dtype=float)
    gyro_rw = np.full(3, np.deg2rad(0.0015) / 3600.0, dtype=float)
    accel_bias = np.full(3, 25e-6 * 9.80665, dtype=float)
    accel_rw = np.full(3, 50e-6 * 9.80665, dtype=float)
    return {
        "gyro_bias": gyro_bias,
        "gyro_rw": gyro_rw,
        "accel_bias": accel_bias,
        "accel_rw": accel_rw,
    }


@dataclass
class NoiseStreams:
    gyro_rng: tuple[np.random.Generator, np.random.Generator, np.random.Generator]
    accel_rng: tuple[np.random.Generator, np.random.Generator, np.random.Generator]
    gnss_rng: tuple[np.random.Generator, np.random.Generator, np.random.Generator, np.random.Generator]
    altimeter_rng: np.random.Generator
    q_gyro_rng: np.random.Generator

    @staticmethod
    def from_seed(seed: int = 12345) -> "NoiseStreams":
        root = np.random.default_rng(seed)
        seeds = root.integers(1, 100_000, size=12)
        gyro_rng = tuple(np.random.default_rng(int(s)) for s in seeds[:3])
        accel_rng = tuple(np.random.default_rng(int(s)) for s in seeds[3:6])
        gnss_rng = tuple(np.random.default_rng(int(s)) for s in seeds[6:10])
        altimeter_rng = np.random.default_rng(int(seeds[10]))
        q_gyro_rng = np.random.default_rng(int(seeds[11]))
        return NoiseStreams(
            gyro_rng=gyro_rng,
            accel_rng=accel_rng,
            gnss_rng=gnss_rng,
            altimeter_rng=altimeter_rng,
            q_gyro_rng=q_gyro_rng,
        )

    def sample_imu(self, dt: float = 1.0) -> ImuNoise:
        ref = imu_reference_sigmas()
        dt_safe = max(float(dt), 1e-12)
        gyro_std = ref["gyro_rw"] / np.sqrt(dt_safe)
        accel_std = ref["accel_rw"] / np.sqrt(dt_safe)
        wn_w = np.array(
            [rng.normal(0.0, gyro_std[i]) for i, rng in enumerate(self.gyro_rng)],
            dtype=float,
        )
        wn_a = np.array(
            [rng.normal(0.0, accel_std[i]) for i, rng in enumerate(self.accel_rng)],
            dtype=float,
        )
        return ImuNoise(wn_w=wn_w, wn_a=wn_a)

    def sample_gnss(self) -> tuple[float, float, float, float]:
        ref = gnss_reference_sigmas()
        fi = self.gnss_rng[0].normal(0.0, ref["phi"])
        lm = self.gnss_rng[1].normal(0.0, ref["lambda"])
        vn = self.gnss_rng[2].normal(0.0, ref["vn"])
        ve = self.gnss_rng[3].normal(0.0, ref["ve"])
        return float(fi), float(lm), float(vn), float(ve)

    def sample_altimeter(self) -> float:
        return float(self.altimeter_rng.normal(0.0, altimeter_reference_sigma()))

    def sample_q_gyro(self) -> float:
        return float(self.q_gyro_rng.normal(0.0, q_gyro_measurement_sigma()))


def q_gyro_measurement_sigma() -> float:
    """1σ шум прямого измерения угловой скорости q, рад/с."""
    return 0.01


def theta_ins_measurement_sigma() -> float:
    """1σ шум оценки тангажа θ_ins, рад."""
    return float(np.deg2rad(0.05))
