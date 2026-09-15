"""Metric-domain screening and safe termination checks."""

import os
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "python"))
from pyblackray import metric_diagnostics, trace_grid


@pytest.mark.skipif(
    not os.environ.get("BLACKRAY_LIBRARY") and not (Path(__file__).parents[1] / "build").exists(),
    reason="build the shared library or set BLACKRAY_LIBRARY",
)
def test_unphysical_metric_regions_are_detectable():
    found_ctc = False
    found_bad_determinant = False
    for radius in np.geomspace(1.01, 20.0, 80):
        for alpha13 in np.linspace(-20.0, 20.0, 81):
            determinant, g_phiphi = metric_diagnostics(
                radius, np.pi / 2.0, spin=0.9, alpha13=alpha13,
            )
            found_ctc |= g_phiphi < 0.0
            found_bad_determinant |= determinant <= 0.0 or not np.isfinite(determinant)
            if found_ctc and found_bad_determinant:
                break
        if found_ctc and found_bad_determinant:
            break
    assert found_ctc or found_bad_determinant, "pathology scan found no invalid metric region"


def test_invalid_ray_inputs_terminate_safely():
    if not os.environ.get("BLACKRAY_LIBRARY") and not (Path(__file__).parents[1] / "build").exists():
        pytest.skip("build the shared library or set BLACKRAY_LIBRARY")
    screen = np.array([0.0])
    records, status = trace_grid(
        screen, screen, observer_radius=3.0, inclination=1.0,
        spin=0.9, alpha13=-20.0, r_isco=1.5, r_out=10.0,
    )
    assert records.shape == (1, 3)
    assert status[0] in (1, 2, 3)
    assert np.all(np.isfinite(records))