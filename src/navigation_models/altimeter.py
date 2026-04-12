from __future__ import annotations


def altimeter_reference_sigma() -> float:
    """1σ высотомера, м."""
    return 1.0


def altimeter_measurement(h: float, noise: float) -> tuple[float, float]:
    return float(h + noise), float(altimeter_reference_sigma())
