# blackray-deformed-kerr: Open-Source Ray-Tracing and Reflection Spectroscopy in Non-Kerr Spacetimes

<!-- [![License: GPL-3.0](https://img.shields.io/badge/License-GPL--3.0-blue.svg)](LICENSE) -->
[![C++](https://img.shields.io/badge/C%2B%2B-17-00599C.svg)](https://isocpp.org/)
[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB.svg)](https://www.python.org/)

## Overview

`blackray-deformed-kerr` is an open-source toolkit for backward ray tracing,
relativistic disk spectroscopy, and local reflection-spectrum convolution in
deformed Kerr spacetimes. Its C++17 core provides metric, Christoffel,
geodesic, ISCO, tetrad-camera, convolver, and OpenMP screen-grid modules.
Python wrappers and notebooks support analysis and publication figures.

reference copies came from [ABHModels/blackray](https://github.com/ABHModels/blackray)
and [ABHModels/raytransfer](https://github.com/ABHModels/raytransfer).

## Mathematical Background

The implementation uses Boyer-Lindquist coordinates and geometric units
$G=c=M=1$. The Johannsen deformation functions are

$$A_1=1+\frac{\alpha_{13}}{r^3},\quad A_2=1+\frac{\alpha_{22}}{r^2},\quad A_5=1+\frac{\alpha_{52}}{r^2},\quad f=\frac{\epsilon_3}{r}.$$

Setting all four parameters to zero recovers Kerr. The code evaluates metric
components, inverse metrics, Christoffel symbols, circular orbits, ISCO
stability, photon geodesics, redshifts, emission angles, and reflection
convolution.

## Installation

Requirements are a C++17 compiler, CMake 3.16+, OpenMP, Python 3.10+, and the
packages in `requirements.txt`.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
cmake -S . -B build -DBLACKRAY_NATIVE=ON
cmake --build build --parallel
```

Use `-DBLACKRAY_NATIVE=OFF` for a portable build. `make` is an equivalent
alternative and creates the same shared library and CLI targets.

If the repository is located on a FAT/VFAT or other mount that does not
preserve Unix executable bits, build into a native Linux filesystem instead:

```bash
cmake -S . -B /tmp/blackray-build -DBLACKRAY_NATIVE=ON
cmake --build /tmp/blackray-build --parallel
/tmp/blackray-build/blackray_cli --help
export BLACKRAY_LIBRARY="/tmp/blackray-build/libblackray_core.so"
```

<!--The project directory under `/media/...` is currently on VFAT with executable
bits masked, so `/tmp/blackray-build` or a directory under `$HOME` is required
for running compiled binaries. -->

## Quickstart

Run these commands from the repository root. The `build/` directory is created
by the CMake configure step; the CLI and Python wrapper cannot run before the
library has been compiled.

### 1. Create the Python environment

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Matplotlib is needed only for `--plot`, while `nbformat` is included for users
who want to inspect or manipulate notebooks programmatically. The notebook
generator itself uses only the Python standard library.

### 2. Configure and build the C++ engine

```bash
cmake -S . -B build -DBLACKRAY_NATIVE=ON
cmake --build build --parallel
```

For a portable build, use `-DBLACKRAY_NATIVE=OFF`. The alternative Makefile
build is:

```bash
make clean
make
```

Confirm that the expected artifacts exist:

```bash
test -x build/blackray_cli
test -f build/libblackray_core.so
```

If either check fails, rerun the configure and build commands. A previous
temporary validation build under `/tmp` does not create the repository-local
`build/` directory.

### 3. Run the C++ command-line ray tracer

```bash
./build/blackray_cli \
	--spin 0.99 \
	--inclination 1.0 \
	--rstep 128 \
	--pstep 128 \
	--r-out 50 \
	--output rays.csv
```

This fires a 128 by 128 observer-screen grid. The command writes
`rays.csv` with the columns `x_obs`, `y_obs`, `r_disk`, `g_factor`,
`cos_theta_e`, and `status`. Status values are `0` for disk hits, `1` for
horizon termination, `2` for escape, and `3` for integration failure.

### 4. Run the Python simulation wrapper

Point the ctypes wrapper at the shared library built in Step 2:

```bash
export BLACKRAY_LIBRARY="$PWD/build/libblackray_core.so"
export PYTHONPATH="$PWD/python"
python python/run_simulation.py \
	--mass 1.0 \
	--spin 0.99 \
	--alpha13 0.5 \
	--epsilon3 0.2 \
	--inclination 1.0 \
	--rstep 128 \
	--pstep 128 \
	--r-out 50 \
	--output rays.dat \
	--spectrum-output spectrum.dat \
	--plot spectrum.png
```

The script writes the ray diagnostics to `rays.dat`, a broadened iron-line
histogram to `spectrum.dat`, and an observer-plane `g`-factor image to
`spectrum.png`. The `--plot` option requires Matplotlib. The `--mass` option
is accepted for workflow compatibility; the current geometric-unit engine
uses $M=1$ and dimensionless radii.

The same parameters can be supplied through the CLI's `key=value` config
file support:

```text
# simulation.cfg
spin=0.99
inclination=1.0
rstep=128
pstep=128
r-out=50
output=rays.csv
```

```bash
./build/blackray_cli --config simulation.cfg
```
<!--
### 5. Generate and open the 14-cell notebook

The generator writes a valid JSON `.ipynb` directly, so `nbformat` is not
required to create it:

```bash
python notebooks/generate_demo_notebook.py
jupyter notebook notebooks/demo_deformed_kerr_analysis.ipynb
```

The notebook contains 14 alternating markdown and Python cells:

1. Title and analysis overview.
2. Imports, library path setup, and figure-directory setup.
3. Kerr/deformed-Kerr case description.
4. The two ray-tracing simulations.
5. Figure 1 explanation.
6. Figure 1 iron-line profile.
7. Figure 2 explanation.
8. Figure 2 radius/redshift joint histogram.
9. Figure 3 explanation.
10. Figure 3 emission-angle histogram.
11. Figure 4 explanation.
12. Figure 4 disk-radius histogram.
13. Figure 5 explanation.
14. Figure 5 observer image-plane impact map.

All five PNG files are saved at 300 DPI under `notebooks/figures/`.

### 6. Run the validation suite

```bash
BLACKRAY_LIBRARY="$PWD/build/libblackray_core.so" \
PYTHONPATH="$PWD/python" \
pytest -q tests
```

The suite checks the Kerr limit, null-norm conservation for 100 deformed-metric
rays, and safe detection/termination of pathological metric regions.
-->
## Repository Layout

```text
include/       Public C++ headers and C ABI
src/           Physics modules, convolver, and CLI
python/        ctypes wrapper and simulation CLI
tests/         Kerr, null-norm, and pathology tests
notebooks/     Notebook generator and generated research demo
data/          XILLVER/RELXILL and other local model tables
docs/          Mathematical and usage documentation
```
<!--
## Legacy Comparison

The original `blackray` contains optimized ray tracing and disk intersection
routines. `RayTransfer` adds transfer-function and reflection-table workflows.
This repository separates those concerns into reusable modules, adds adaptive
integration and null checks, supports Johannsen parameters, and exposes Python
and CLI interfaces while preserving the scientific workflow. Cite the
ABHModels repositories and relevant papers when using their algorithms or data.

## Citation

Until a DOI and `CITATION.cff` are published, cite this repository URL, release
or commit hash, and access date. Also cite the Johannsen metric paper, the
XILLVER/RELXILL model papers and table release, the ABHModels repositories, and
any paper defining the transfer function or emissivity prescription used.
Verify exact bibliographic metadata for the versions used.

## License

The intended project license is GPL-3.0. Add the complete `LICENSE` text before
publication and confirm that retained legacy code and model tables are
redistributed under compatible terms.

## Preparing the First Commit

Run these commands from the repository root only after reviewing provenance and
licenses:

```bash
git status
git diff -- README.md requirements.txt .gitignore CMakeLists.txt Makefile
rm -rf _legacy_blackray _legacy_raytransfer
git add .
git commit -m "Build deformed Kerr ray tracing toolkit"
git branch -M main
gh auth status
gh repo create <your-github-username>/blackray-deformed-kerr --public --source=. --remote=origin \
	--description "Ray tracing and reflection spectroscopy in deformed Kerr spacetimes"
git push -u origin main
```

Replace `<your-github-username>` with your GitHub account. `gh repo create`
uses the account selected by `gh auth status`.-->
