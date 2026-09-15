"""Production-style validation checks for the deformed-Kerr solver."""

import os
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "python"))
from pyblackray import find_isco, metric_diagnostics, trace_grid


@pytest.mark.skipif(
    not os.environ.get("BLACKRAY_LIBRARY") and not (Path(__file__).parents[1] / "build").exists(),
    reason="build the shared library or set BLACKRAY_LIBRARY",
)
def test_zero_deformation_matches_kerr_across_equatorial_grid():
    radii = np.array([1.5, 2.0, 3.0, 5.0, 10.0, 20.0], dtype=float)
    spin_values = [0.2, 0.5, 0.9, 0.99]
    theta = np.pi / 2.0
    max_det_error = 0.0
    max_g_error = 0.0

    for spin in spin_values:
        for radius in radii:
            det_model, g_model = metric_diagnostics(
                radius,
                theta,
                spin=spin,
                alpha13=0.0,
                alpha22=0.0,
                alpha52=0.0,
                epsilon3=0.0,
            )
            sigma = radius * radius
            analytical_det = -sigma * sigma
            analytical_g = radius * radius + spin * spin + 2.0 * spin * spin / radius
            max_det_error = max(max_det_error, abs(det_model - analytical_det))
            max_g_error = max(max_g_error, abs(g_model - analytical_g))

    assert np.isfinite(max_det_error)
    assert np.isfinite(max_g_error)
    assert max_det_error < 1.0e-8
    assert max_g_error < 1.0e-8


@pytest.mark.skipif(
    not os.environ.get("BLACKRAY_LIBRARY") and not (Path(__file__).parents[1] / "build").exists(),
    reason="build the shared library or set BLACKRAY_LIBRARY",
)
def test_trace_grid_hit_statistics_stabilize_across_screen_resolutions():
    spin = 0.99
    radius = 100.0
    inclination = 1.0

    x8 = np.linspace(-8.0, 8.0, 8)
    y8 = np.linspace(-8.0, 8.0, 8)
    xx8, yy8 = np.meshgrid(x8, y8, indexing="ij")

    x16 = np.linspace(-8.0, 8.0, 16)
    y16 = np.linspace(-8.0, 8.0, 16)
    xx16, yy16 = np.meshgrid(x16, y16, indexing="ij")

    records8, status8 = trace_grid(
        xx8.ravel(), yy8.ravel(), observer_radius=radius,
        inclination=inclination, spin=spin, r_isco=find_isco(spin)["radius"], r_out=30.0,
    )
    records16, status16 = trace_grid(
        xx16.ravel(), yy16.ravel(), observer_radius=radius,
        inclination=inclination, spin=spin, r_isco=find_isco(spin)["radius"], r_out=30.0,
    )

    coarse_hits = records8[status8 == 0]
    fine_hits = records16[status16 == 0]

    assert coarse_hits.shape[0] > 0
    assert fine_hits.shape[0] > 0

    coarse_rate = coarse_hits.shape[0] / records8.shape[0]
    fine_rate = fine_hits.shape[0] / records16.shape[0]

    assert abs(fine_rate - coarse_rate) < 0.20


@pytest.mark.skipif(
    not os.environ.get("BLACKRAY_LIBRARY") and not (Path(__file__).parents[1] / "build").exists(),
    reason="build the shared library or set BLACKRAY_LIBRARY",
)
@pytest.mark.parametrize("spin", [0.2, 0.5, 0.8, 0.99])
def test_isco_matches_analytic_kerr_formula(spin):
    z1 = 1.0 + (1.0 - spin * spin) ** (1.0 / 3.0) * (
        (1.0 + spin) ** (1.0 / 3.0) + (1.0 - spin) ** (1.0 / 3.0)
    )
    z2 = np.sqrt(3.0 * spin * spin + z1 * z1)
    expected = 3.0 + z2 - np.sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2))

    actual = find_isco(spin=spin)["radius"]

    assert np.isfinite(actual)
    assert abs(actual - expected) < 1.0e-5
