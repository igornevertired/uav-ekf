from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ImuNoise:
    wn_w: np.ndarray
    wn_a: np.ndarray


@dataclass
class NoiseStreams:
    gyro_rng: tuple[np.random.Generator, np.random.Generator, np.random.Generator]
    accel_rng: tuple[np.random.Generator, np.random.Generator, np.random.Generator]
    gnss_rng: tuple[np.random.Generator, np.random.Generator, np.random.Generator, np.random.Generator]

    @staticmethod
    def from_seed(seed: int = 12345) -> "NoiseStreams":
        root = np.random.default_rng(seed)
        seeds = root.integers(1, 100_000, size=10)
        gyro_rng = tuple(np.random.default_rng(int(s)) for s in seeds[:3])
        accel_rng = tuple(np.random.default_rng(int(s)) for s in seeds[3:6])
        gnss_rng = tuple(np.random.default_rng(int(s)) for s in seeds[6:10])
        return NoiseStreams(gyro_rng=gyro_rng, accel_rng=accel_rng, gnss_rng=gnss_rng)

    def sample_imu(self) -> ImuNoise:
        wn_w = np.array([rng.normal(0.0, 1e-4) for rng in self.gyro_rng], dtype=float)
        wn_a = np.array([rng.normal(0.0, 1e-4) for rng in self.accel_rng], dtype=float)
        return ImuNoise(wn_w=wn_w, wn_a=wn_a)

    def sample_gnss(self) -> tuple[float, float, float, float]:
        ang_std = (1.0 / 111_000.0) * np.pi / 180.0
        fi = self.gnss_rng[0].normal(0.0, ang_std)
        lm = self.gnss_rng[1].normal(0.0, ang_std)
        vn = self.gnss_rng[2].normal(0.0, 1e-2)
        ve = self.gnss_rng[3].normal(0.0, 1e-2)
        return float(fi), float(lm), float(vn), float(ve)
