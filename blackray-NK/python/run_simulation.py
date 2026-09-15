"""Run a ray-tracing grid and optionally convolve it with a local spectrum."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from pyblackray import convolve_spectrum, trace_grid


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mass", type=float, default=1.0)
    parser.add_argument("--spin", type=float, default=0.5)
    parser.add_argument("--alpha13", type=float, default=0.0)
    parser.add_argument("--alpha22", type=float, default=0.0)
    parser.add_argument("--alpha52", type=float, default=0.0)
    parser.add_argument("--epsilon3", type=float, default=0.0)
    parser.add_argument("--rstep", type=int, default=128)
    parser.add_argument("--pstep", type=int, default=128)
    parser.add_argument("--inclination", type=float, default=1.0)
    parser.add_argument("--observer-radius", type=float, default=1000.0)
    parser.add_argument("--screen-radius", type=float, default=20.0)
    parser.add_argument("--r-isco", type=float, default=0.0)
    parser.add_argument("--r-out", type=float, default=50.0)
    parser.add_argument("--output", type=Path, default=Path("rays.dat"))
    parser.add_argument("--spectrum-output", type=Path, default=Path("spectrum.dat"))
    parser.add_argument("--line-energy", type=float, default=6.4)
    parser.add_argument("--spectrum-bins", type=int, default=400)
    parser.add_argument("--spectrum-table", type=Path, default=None,
                       help="XILLVER/RELXILL grid file to convolve with the ray table")
    parser.add_argument("--emissivity-index", type=float, default=-3.0)
    parser.add_argument("--angle-mode", choices=["ray_traced", "fixed"],
                       default="ray_traced")
    parser.add_argument("--fixed-cos-inclination", type=float, default=0.5)
    parser.add_argument("--plot", type=Path, default=None)
    args = parser.parse_args()
    if args.mass <= 0 or args.rstep < 1 or args.pstep < 1:
        parser.error("mass, rstep, and pstep must be positive")

    x = np.linspace(-args.screen_radius, args.screen_radius, args.rstep)
    y = np.linspace(-args.screen_radius, args.screen_radius, args.pstep)
    xx, yy = np.meshgrid(x, y, indexing="ij")
    records, status = trace_grid(
        xx.ravel(), yy.ravel(), observer_radius=args.observer_radius,
        inclination=args.inclination, spin=args.spin, alpha13=args.alpha13,
        alpha22=args.alpha22, alpha52=args.alpha52, epsilon3=args.epsilon3,
        r_isco=args.r_isco, r_out=args.r_out,
    )
    table = np.column_stack((xx.ravel(), yy.ravel(), records, status))
    np.savetxt(args.output, table,
               header="x_obs y_obs r_disk g_factor cos_theta_e status")
    valid = (status == 0) & (records[:, 0] > 0.0) & (records[:, 1] > 0.0)
    valid_rays = records[valid]

    if args.spectrum_table is not None and args.spectrum_table.exists():
        # Convolve with a real local spectrum table (e.g., XILLVER grid)
        from pyblackray import load_xillver_spectrum
        spec_energies, spec_cos_incl, spec_fluxes = load_xillver_spectrum(
            str(args.spectrum_table)
        )
        energy_lo, energy_hi, flux = convolve_spectrum(
            valid_rays[:, 0], valid_rays[:, 1], valid_rays[:, 2],
            spec_energies, spec_cos_incl, spec_fluxes,
            emissivity_index=args.emissivity_index,
            fixed_cos_inclination=args.fixed_cos_inclination,
            angle_mode=args.angle_mode,
            output_bins=args.spectrum_bins,
        )
        np.savetxt(args.spectrum_output,
                   np.column_stack((energy_lo, energy_hi, flux)),
                   header="energy_lo energy_hi flux")
    else:
        # Fallback: simple broadened line histogram (no local spectrum table)
        energies = args.line_energy * valid_rays[:, 1]
        weights = valid_rays[:, 1] ** 3
        edges = np.linspace(args.line_energy * 0.1, args.line_energy * 1.5,
                            args.spectrum_bins + 1)
        histogram, _ = np.histogram(energies, bins=edges, weights=weights)
        centers = 0.5 * (edges[1:] + edges[:-1])
        np.savetxt(args.spectrum_output, np.column_stack((centers, histogram)),
                   header="energy_keV flux")
    if args.plot is not None:
        try:
            import matplotlib.pyplot as plt
        except ImportError as error:
            raise RuntimeError(
                "--plot requires matplotlib; install requirements.txt"
            ) from error
        image = records[:, 1].reshape(xx.shape)
        plt.imshow(image, origin="lower", extent=[x.min(), x.max(), y.min(), y.max()])
        plt.colorbar(label="g factor")
        plt.xlabel("x_obs")
        plt.ylabel("y_obs")
        plt.tight_layout()
        plt.savefig(args.plot, dpi=160)
        plt.close()


if __name__ == "__main__":
    main()