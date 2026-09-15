"""Regression checks for the zero-deformation Kerr limit."""

import os
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "python"))
from pyblackray import trace_grid


def kerr_isco(spin: float) -> float:
    z1 = 1.0 + (1.0 - spin * spin) ** (1.0 / 3.0) * (
        (1.0 + spin) ** (1.0 / 3.0) + (1.0 - spin) ** (1.0 / 3.0)
    )
    z2 = np.sqrt(3.0 * spin * spin + z1 * z1)
    return 3.0 + z2 - np.sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2))


@pytest.mark.skipif(
    not os.environ.get("BLACKRAY_LIBRARY") and not (Path(__file__).parents[1] / "build").exists(),
    reason="build the shared library or set BLACKRAY_LIBRARY",
)
def test_kerr_isco_and_transfer_profile():
    spin = 0.99
    expected_isco = kerr_isco(spin)
    assert abs(expected_isco - 1.4544979381) < 1.0e-4

    screen = np.linspace(-8.0, 8.0, 9)
    x_obs, y_obs = np.meshgrid(screen, screen, indexing="ij")
    records, status = trace_grid(
        x_obs.ravel(), y_obs.ravel(), observer_radius=100.0,
        inclination=1.0, spin=spin, r_isco=expected_isco, r_out=30.0,
    )
    hits = records[status == 0]
    assert len(hits) > 0, "Kerr ray trace produced no disk hits"
    g_min = float(np.min(hits[:, 1]))
    g_max = float(np.max(hits[:, 1]))

    # Fixed-grid reference envelope for this observer and screen geometry.
    assert 0.15 - 1.0e-5 <= g_min <= 0.35 + 1.0e-5
    assert 1.0 - 1.0e-5 <= g_max <= 1.8 + 1.0e-5

    energies = 6.4 * hits[:, 1]
    histogram, edges = np.histogram(energies, bins=32, range=(0.5, 8.0))
    profile = histogram / max(float(histogram.sum()), 1.0)
    assert np.isclose(profile.sum(), 1.0, atol=1.0e-5)
    assert np.all(profile >= 0.0)
    assert np.count_nonzero(profile) >= 2