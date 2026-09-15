# BlackRay Deformed Kerr — Practical Run Guide

This guide is meant to be a simple, friend-friendly starting point for running the project on a laptop such as the Dell Latitude 7400, or on a similar machine, using the already-built files when available.

## 1) First check: do you already have compiled files?

From the project root, run:

```bash
cd /media/user/DUALBOOT/TVM/numGR/projects/blackray-m.bambi/blackray-deformed-kerr
ls -l build/blackray_cli build/libblackray_core.so
```

If those files exist, you already have the compiled C++ core and CLI. In that case, you do not need to rebuild again unless you changed the C++ code.

If they do not exist, build first:

```bash
cmake -S . -B build -DBLACKRAY_NATIVE=ON
cmake --build build --parallel
```

If the filesystem is blocking executable bits, use a native Linux location instead, for example:

```bash
cmake -S . -B /tmp/blackray-build -DBLACKRAY_NATIVE=ON
cmake --build /tmp/blackray-build --parallel
```

Then set:

```bash
export BLACKRAY_LIBRARY="/tmp/blackray-build/libblackray_core.so"
```

---

## 2) Set the runtime environment

For the Python wrapper and notebook run, set these in the terminal before launching any Python or Jupyter command:

```bash
cd /media/user/DUALBOOT/TVM/numGR/projects/blackray-m.bambi/blackray-deformed-kerr
export BLACKRAY_LIBRARY="$PWD/build/libblackray_core.so"
export PYTHONPATH="$PWD/python"
export OMP_NUM_THREADS=4
```

### Where to set `OMP_NUM_THREADS`

Set `OMP_NUM_THREADS` in the shell before running the command, not inside the compiled C++ code.

Recommended places:

- in your terminal session before `python`, `jupyter`, or `./build/blackray_cli`
- in a job script on an HPC system
- in the shell that launches the notebook

Best practice:

```bash
export OMP_NUM_THREADS=4
python notebooks/generate_demo_notebook.py
```

For a heavier run on a stronger machine, you can use:

```bash
export OMP_NUM_THREADS=8
```

On a laptop, a moderate value such as `4` is usually safer than using all cores blindly.

---

## 3) Use the already-built files if they already exist

If the repository already contains the compiled artifacts and you are currently running Python cells such as `run_aux_case(...)` and `trace_grid(...)`, then you can skip rebuilding and just do:

```bash
cd /media/user/DUALBOOT/TVM/numGR/projects/blackray-m.bambi/blackray-deformed-kerr
export BLACKRAY_LIBRARY="$PWD/build/libblackray_core.so"
export PYTHONPATH="$PWD/python"
export OMP_NUM_THREADS=4
python -c "from pyblackray import trace_grid; print('ready')"
```

If that prints `ready`, the existing build is usable.

---

## 4) What the notebook is doing right now

The notebook cells already use the Python wrapper `trace_grid(...)`, which in turn talks to the compiled shared library `libblackray_core.so`.

So when you run cells such as:

- `run_case(...)`
- `run_aux_case(...)`
- `trace_grid(...)`

you are already using the existing compiled C++ code.

The slow part is not the build. The slow part is the repeated numerical work performed by `trace_grid(...)` for many screen points and many parameter combinations.

---

## 5) What the main parameters mean

These are the most important parameters you will see in the notebook and command line:

- `rstep`: number of radial screen samples
- `pstep`: number of azimuthal screen samples
- `screen-radius`: half-width of the observer screen
- `observer-radius`: observer location
- `inclination`: observer inclination
- `spin`: dimensionless spin
- `alpha13`, `alpha22`, `alpha52`, `epsilon3`: deformation parameters
- `r_out`: outer disk radius

Important note:

- In the rewritten code, `rstep` and `pstep` are screen-sample counts, not the older geometric increment-style settings.

---

## 6) Best way to run the notebook on a laptop

### Option A: run only the cells you need

This is the safest choice for the Dell Latitude 7400.

Run the notebook in chunks:

1. first the imports and case setup,
2. then the base Kerr/deformed case,
3. then only one or two figures,
4. then the supplementary sweeps if you really need them.

### Option B: reduce the grid size first

If the notebook is too slow, reduce `rstep` and `pstep` before running bigger figures.

For example, in the notebook cells you can change:

```python
rstep, pstep = 16, 16
```

to something smaller such as:

```python
rstep, pstep = 8, 8
```

This reduces the number of traced rays immediately.

---

## 7) Run the CLI directly

If you want a quick direct run without the notebook:

```bash
cd /media/user/DUALBOOT/TVM/numGR/projects/blackray-m.bambi/blackray-deformed-kerr
export BLACKRAY_LIBRARY="$PWD/build/libblackray_core.so"
export PYTHONPATH="$PWD/python"
export OMP_NUM_THREADS=4
./build/blackray_cli \
  --spin 0.99 \
  --inclination 1.0 \
  --rstep 128 \
  --pstep 128 \
  --r-out 50 \
  --output rays.csv
```

This will run the compiled CLI directly.

---

## 8) Run the Python wrapper directly

```bash
cd /media/user/DUALBOOT/TVM/numGR/projects/blackray-m.bambi/blackray-m.bambi/blackray-deformed-kerr
export BLACKRAY_LIBRARY="$PWD/build/libblackray_core.so"
export PYTHONPATH="$PWD/python"
export OMP_NUM_THREADS=4
python python/run_simulation.py \
  --spin 0.99 \
  --alpha13 0.5 \
  --epsilon3 0.2 \
  --inclination 1.0 \
  --rstep 64 \
  --pstep 64 \
  --r-out 50 \
  --output rays.dat \
  --spectrum-output spectrum.dat \
  --plot spectrum.png
```

This is useful for a smaller batch run.

---

## 9) Launch the notebook

From the project root:

```bash
cd /media/user/DUALBOOT/TVM/numGR/projects/blackray-m.bambi/blackray-deformed-kerr
export BLACKRAY_LIBRARY="$PWD/build/libblackray_core.so"
export PYTHONPATH="$PWD/python"
export OMP_NUM_THREADS=4
jupyter notebook notebooks/demo_deformed_kerr_analysis.ipynb
```

Then run the notebook cells one by one.

---

## 10) When to restart the notebook kernel

Restart the notebook kernel if:

- you changed the compiled library or rebuilt the code,
- you changed `rstep` / `pstep` and want a clean state,
- memory use looks high,
- a figure cell is taking unusually long or appears stuck,
- or you want a clean run after a large sweep.

A restart is not required for every small run, but it is a good idea before a very large sweep.

---

## 11) Practical recommendation for your current situation

Because you already have some Python figure work running with `run_aux_case(...)` and `trace_grid(...)`, the most practical next move is:

1. keep the existing compiled files,
2. set `OMP_NUM_THREADS=4`,
3. run the notebook in chunks,
4. reduce `rstep` / `pstep` for exploratory plots,
5. only use the big sweep cells when you are ready for a longer run.

---

## 12) Quick start command block

If you want the shortest possible starting command set, use this:

```bash
cd /media/user/DUALBOOT/TVM/numGR/projects/blackray-m.bambi/blackray-deformed-kerr
export BLACKRAY_LIBRARY="$PWD/build/libblackray_core.so"
export PYTHONPATH="$PWD/python"
export OMP_NUM_THREADS=4
jupyter notebook notebooks/demo_deformed_kerr_analysis.ipynb
```

If you want only a quick CLI test:

```bash
cd /media/user/DUALBOOT/TVM/numGR/projects/blackray-m.bambi/blackray-deformed-kerr
export BLACKRAY_LIBRARY="$PWD/build/libblackray_core.so"
export PYTHONPATH="$PWD/python"
export OMP_NUM_THREADS=4
./build/blackray_cli --spin 0.99 --inclination 1.0 --rstep 32 --pstep 32 --r-out 50 --output rays.csv
```

---

## 13) A friend-friendly note

If your friend pulls this repository from the commit and uses the same environment, they can start from the same steps above:

- check for `build/blackray_cli` and `build/libblackray_core.so`
- if missing, build once
- then set `BLACKRAY_LIBRARY`, `PYTHONPATH`, and `OMP_NUM_THREADS`
- then launch the notebook or use the CLI

This should work well for a similar Asus TUF laptop as long as the system has the required Python and C++ toolchain installed.
