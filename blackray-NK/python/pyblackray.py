"""ctypes wrapper for the BlackRay C ABI and lightweight spectrum helpers."""

from __future__ import annotations

import ctypes
import os
from pathlib import Path

import numpy as np


def _load_library() -> ctypes.CDLL:
    root = Path(__file__).resolve().parents[1]
    configured = os.environ.get("BLACKRAY_LIBRARY")
    if configured:
        return ctypes.CDLL(configured)
    candidates = [
        root / "build" / "libblackray_core.so",
        root / "build" / "libblackray_core.dylib",
        root / "build" / "blackray_core.dll",
    ]
    for path in candidates:
        if path.exists():
            return ctypes.CDLL(str(path))
    raise FileNotFoundError(
        "libblackray_core not found; build with CMake or set BLACKRAY_LIBRARY"
    )


def _configure(library: ctypes.CDLL) -> None:
    function = library.blackray_trace_grid
    function.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int,
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int),
    ]
    function.restype = ctypes.c_int
    null_norm = library.blackray_max_null_norm
    null_norm.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int,
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_int, ctypes.POINTER(ctypes.c_double),
    ]
    null_norm.restype = ctypes.c_int
    diagnostics = library.blackray_metric_diagnostics
    diagnostics.argtypes = [
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ]
    diagnostics.restype = ctypes.c_int


def trace_grid(x_obs, y_obs, *, observer_radius=1000.0, inclination=1.0,
               spin=0.5, alpha13=0.0, alpha22=0.0, alpha52=0.0,
               epsilon3=0.0, r_isco=0.0, r_out=50.0, library=None):
    """Trace screen coordinates and return an array of ray-hit records."""
    x = np.ascontiguousarray(x_obs, dtype=np.float64).ravel()
    y = np.ascontiguousarray(y_obs, dtype=np.float64).ravel()
    if x.size != y.size:
        raise ValueError("x_obs and y_obs must have the same size")
    r_disk = np.empty(x.size, dtype=np.float64)
    g_factor = np.empty(x.size, dtype=np.float64)
    cos_theta_e = np.empty(x.size, dtype=np.float64)
    status = np.empty(x.size, dtype=np.int32)
    lib = library or _load_library()
    _configure(lib)
    code = lib.blackray_trace_grid(
        x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), int(x.size),
        observer_radius, inclination, spin, alpha13, alpha22, alpha52,
        epsilon3, r_isco, r_out,
        r_disk.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        g_factor.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        cos_theta_e.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        status.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
    )
    if code != 0:
        raise RuntimeError(f"BlackRay C API failed with code {code}")
    return np.column_stack((r_disk, g_factor, cos_theta_e)), status


def max_null_norm(x_obs, y_obs, *, observer_radius=1000.0, inclination=1.0,
                  spin=0.5, alpha13=0.0, alpha22=0.0, alpha52=0.0,
                  epsilon3=0.0, step=0.01, steps=100, library=None):
    x = np.ascontiguousarray(x_obs, dtype=np.float64).ravel()
    y = np.ascontiguousarray(y_obs, dtype=np.float64).ravel()
    if x.size != y.size:
        raise ValueError("x_obs and y_obs must have the same size")
    maximum = ctypes.c_double()
    lib = library or _load_library()
    _configure(lib)
    code = lib.blackray_max_null_norm(
        x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), int(x.size),
        observer_radius, inclination, spin, alpha13, alpha22, alpha52,
        epsilon3, step, steps, ctypes.byref(maximum),
    )
    if code != 0:
        raise RuntimeError(f"BlackRay null-norm API failed with code {code}")
    return maximum.value


def find_isco(spin, *, alpha13=0.0, alpha22=0.0, alpha52=0.0,
              epsilon3=0.0, library=None):
    """Return the analytic Kerr ISCO radius for the zero-deformation benchmark."""
    spin = float(spin)
    if not np.isfinite(spin) or abs(spin) >= 1.0:
        raise ValueError("spin must satisfy |spin| < 1")
    if not np.isclose(alpha13, 0.0, atol=1.0e-12) or not np.isclose(
        alpha22, 0.0, atol=1.0e-12
    ) or not np.isclose(alpha52, 0.0, atol=1.0e-12) or not np.isclose(
        epsilon3, 0.0, atol=1.0e-12
    ):
        raise ValueError(
            "find_isco currently supports the zero-deformation Kerr benchmark only"
        )
    z1 = 1.0 + (1.0 - spin * spin) ** (1.0 / 3.0) * (
        (1.0 + spin) ** (1.0 / 3.0) + (1.0 - spin) ** (1.0 / 3.0)
    )
    z2 = np.sqrt(3.0 * spin * spin + z1 * z1)
    return float(
        3.0 + z2 - np.sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2))
    )


def metric_diagnostics(r, theta, *, spin=0.5, alpha13=0.0, alpha22=0.0,
                       alpha52=0.0, epsilon3=0.0, library=None):
    determinant = ctypes.c_double()
    g_phiphi = ctypes.c_double()
    lib = library or _load_library()
    _configure(lib)
    code = lib.blackray_metric_diagnostics(
        r, theta, spin, alpha13, alpha22, alpha52, epsilon3,
        ctypes.byref(determinant), ctypes.byref(g_phiphi),
    )
    if code != 0:
        raise RuntimeError(f"BlackRay metric API failed with code {code}")
    return determinant.value, g_phiphi.value