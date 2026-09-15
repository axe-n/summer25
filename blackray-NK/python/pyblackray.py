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

    isco = library.blackray_find_isco
    isco.argtypes = [
        ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int),
    ]
    isco.restype = ctypes.c_int

    conv = library.blackray_convolve_spectrum
    conv.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int,
        ctypes.c_double, ctypes.c_double, ctypes.c_int,
        ctypes.c_int, ctypes.c_double, ctypes.c_double,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    conv.restype = ctypes.c_int


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


def _analytic_kerr_isco(spin: float) -> float:
    """Analytic Kerr ISCO radius (prograde) for the zero-deformation benchmark."""
    z1 = 1.0 + (1.0 - spin * spin) ** (1.0 / 3.0) * (
        (1.0 + spin) ** (1.0 / 3.0) + (1.0 - spin) ** (1.0 / 3.0)
    )
    z2 = np.sqrt(3.0 * spin * spin + z1 * z1)
    return float(3.0 + z2 - np.sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2)))


def find_isco(spin, *, alpha13=0.0, alpha22=0.0, alpha52=0.0,
              epsilon3=0.0, library=None, inner_radius=0.0,
              outer_radius=50.0, radial_scan_step=0.02,
              root_tolerance=1.0e-10, maximum_iterations=100):
    """Find the ISCO radius for the given Johannsen deformation parameters.

    When all deformation parameters are zero, this calls the C++ solver which
    will return the same value as the analytic Kerr formula
    (_analytic_kerr_isco). For non-zero deformation parameters, the C++ solver
    performs a radial bisection on the marginal-stability condition.

    Returns a dict with keys: radius, omega, u_t, energy, angular_momentum,
    radial_curvature, vertical_curvature, vertical_instability_warning.
    """
    spin = float(spin)
    if not np.isfinite(spin) or abs(spin) >= 1.0:
        raise ValueError("spin must satisfy |spin| < 1")
    if not np.isfinite(alpha13) or not np.isfinite(alpha22) or \
       not np.isfinite(alpha52) or not np.isfinite(epsilon3):
        raise ValueError("deformation parameters must be finite")

    lib = library or _load_library()
    _configure(lib)

    r_isco = ctypes.c_double()
    omega = ctypes.c_double()
    u_t = ctypes.c_double()
    energy = ctypes.c_double()
    angular_momentum = ctypes.c_double()
    radial_curvature = ctypes.c_double()
    vertical_curvature = ctypes.c_double()
    vertical_instability = ctypes.c_int()

    code = lib.blackray_find_isco(
        spin, alpha13, alpha22, alpha52, epsilon3,
        inner_radius, outer_radius, radial_scan_step,
        root_tolerance, maximum_iterations,
        ctypes.byref(r_isco), ctypes.byref(omega), ctypes.byref(u_t),
        ctypes.byref(energy), ctypes.byref(angular_momentum),
        ctypes.byref(radial_curvature), ctypes.byref(vertical_curvature),
        ctypes.byref(vertical_instability),
    )
    if code != 0:
        raise RuntimeError(f"BlackRay ISCO API failed with code {code}")

    return {
        "radius": r_isco.value,
        "omega": omega.value,
        "u_t": u_t.value,
        "energy": energy.value,
        "angular_momentum": angular_momentum.value,
        "radial_curvature": radial_curvature.value,
        "vertical_curvature": vertical_curvature.value,
        "vertical_instability_warning": bool(vertical_instability.value),
    }


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


def load_xillver_spectrum(path, *, max_angles=None):
    """Load an XILLVER/RELXILL grid file into arrays for convolve_spectrum.

    The file format is whitespace-separated columns:
        E_lo  E_hi  flux(angle=0)  flux(angle=1)  ...  flux(angle=N-1)

    Returns (energies, cos_inclinations, fluxes) where:
      - energies: bin-edge energies (length n_energies + 1)
      - cos_inclinations: cos(theta) for each angle column (length n_angles)
      - fluxes: row-major flattened flux[angle * n_energies + energy]
    """
    import numpy as np
    data = np.loadtxt(path, ndmin=2)
    n_angles = data.shape[1] - 2  # first two cols are E_lo, E_hi
    if max_angles is not None:
        n_angles = min(n_angles, max_angles)
    energies = np.concatenate([data[0, 0:2], data[-1, 1:2]])
    # Reconstruct full energy bin edges
    e_lo = data[:, 0]
    e_hi = data[:, 1]
    full_energies = np.empty(len(data) + 1)
    full_energies[:-1] = e_lo
    full_energies[-1] = e_hi[-1]
    cos_inclinations = np.array([
        np.cos(np.deg2rad(18.194874 + i * (87.13402 - 18.194874) / max(n_angles - 1, 1)))
        for i in range(n_angles)
    ])
    fluxes = data[:, 2:2 + n_angles].ravel(order="F")  # column-major per angle
    return full_energies, cos_inclinations, fluxes


def convolve_spectrum(rays_r_disk, rays_g_factor, rays_cos_theta_e,
                      spec_energies, spec_cos_inclinations, spec_fluxes,
                      *, emissivity_index=-3.0, fixed_cos_inclination=0.5,
                      angle_mode="ray_traced", output_bins=2000,
                      output_energy_min=0.00035, output_energy_max=2000.0,
                      library=None):
    """Convolve ray-hit observables with a local reflection spectrum.

    Parameters
    ----------
    rays_r_disk, rays_g_factor, rays_cos_theta_e : array-like
        Per-ray disk radius, redshift factor, and emission-angle cosine.
    spec_energies : array-like
        Energy bin edges of the local spectrum table (length n_energies + 1).
    spec_cos_inclinations : array-like
        cos(theta) values for each angle slice (length n_angles).
    spec_fluxes : array-like
        Flux values, row-major: [angle0_energy0..angle0_energyN-1,
                                 angle1_energy0..angle1_energyN-1, ...].
    angle_mode : {"ray_traced", "fixed"}
        Whether to use the ray-traced emission angle or a fixed inclination.
    """
    e_arr = np.ascontiguousarray(spec_energies, dtype=np.float64)
    c_arr = np.ascontiguousarray(spec_cos_inclinations, dtype=np.float64)
    f_arr = np.ascontiguousarray(spec_fluxes, dtype=np.float64)
    r_arr = np.ascontiguousarray(rays_r_disk, dtype=np.float64).ravel()
    g_arr = np.ascontiguousarray(rays_g_factor, dtype=np.float64).ravel()
    ct_arr = np.ascontiguousarray(rays_cos_theta_e, dtype=np.float64).ravel()

    n_energies = len(e_arr) - 1
    n_angles = len(c_arr)
    if len(f_arr) != n_energies * n_angles:
        raise ValueError(
            f"spec_fluxes has {len(f_arr)} entries but expected "
            f"{n_energies} energies x {n_angles} angles = {n_energies * n_angles}"
        )

    energy_lo = np.empty(output_bins, dtype=np.float64)
    energy_hi = np.empty(output_bins, dtype=np.float64)
    flux = np.empty(output_bins, dtype=np.float64)

    lib = library or _load_library()
    _configure(lib)

    mode = 1 if angle_mode == "fixed" else 0
    code = lib.blackray_convolve_spectrum(
        r_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        g_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ct_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        int(r_arr.size),
        e_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        c_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        f_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        int(n_energies), int(n_angles),
        emissivity_index, fixed_cos_inclination, mode,
        output_bins, output_energy_min, output_energy_max,
        energy_lo.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        energy_hi.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        flux.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    )
    if code != 0:
        raise RuntimeError(f"BlackRay convolve API failed with code {code}")
    return energy_lo, energy_hi, flux