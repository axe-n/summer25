"""Null-geodesic conservation test for a deformed Johannsen metric."""

import os
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "python"))
from pyblackray import max_null_norm


@pytest.mark.skipif(
    not os.environ.get("BLACKRAY_LIBRARY") and not (Path(__file__).parents[1] / "build").exists(),
    reason="build the shared library or set BLACKRAY_LIBRARY",
)
def test_deformed_photon_null_norm_over_100_rays():
    screen = np.linspace(-10.0, 10.0, 10)
    x_obs, y_obs = np.meshgrid(screen, screen, indexing="ij")
    maximum = max_null_norm(
        x_obs.ravel(), y_obs.ravel(), observer_radius=100.0,
        inclination=1.0, spin=0.7, alpha13=0.5, epsilon3=0.2,
        step=1.0e-2, steps=100,
    )
    assert np.isfinite(maximum)
    assert maximum < 1.0e-7