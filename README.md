# Two-phase SPH (water / oil)

3D weakly-compressible SPH solver for two fluid phases, written in C++ with OpenMP.
All particles use one kernel (Monaghan cubic spline, support radius `h`) for density,
pressure, viscosity and surface tension. See [kernel.h](kernel.h) and [sph.cpp](sph.cpp).

## Requirements

- `g++` with OpenMP support
- `make`
- Python 3 and `ffmpeg` (for the PyVista animations)

## Build and run

```bash
make            # build, run the simulation, then render animation.mp4 with PyVista
make run        # build and run the simulation only
make fresh      # delete old *.dat first, then run (use this for a clean restart)
make animate    # render the existing *.dat files to an mp4 with render_pyvista.py
make clean      # delete the binary and *.dat (this destroys finished runs)
```

`make animate` takes overrides:

```bash
make animate DIR=part_1 OUT=part1.mp4
make animate DIR=part_2 OUT=part2.mp4 ARGS="--opacity 0.5 --rotate"
```

Or build by hand:

```bash
g++ -O2 -fopenmp sph.cpp -o sph -lm
./sph
```

Control the number of threads with `OMP_NUM_THREADS`, e.g. `OMP_NUM_THREADS=4 ./sph`.

The run is long: `TS = 1000000` steps of `dt = 0.00005`. `make` runs the whole
simulation before it renders anything, so use `make run` and `make animate`
separately if you want to look at frames as they appear.

## What the run produces

- Every 1000 steps the code writes `NNN.dat` (`001.dat`, `002.dat`, ...) **in the
  current directory**. With the default settings that is 1000 frames.
- Each `.dat` file has one particle per line, comma separated:

  ```
  x, y, z, vx, vy, vz, phase
  ```

  `phase` is 1 for water and 2 for oil.
- The same step prints one line to stdout:
  `count  KE=...  rho/rho0=...  water_z=...  oil_z=...`.
  `rho/rho0` going to 1.0 with `KE` decaying means the fluid has settled.

Frame numbering restarts at `001` on every run, so an aborted run leaves the tail of
the previous run behind. Use `make fresh`, or move finished frames into a folder
(`part_1/`, `part_2/`, ...) before starting the next run.

## Initial conditions

`main()` in [sph.cpp](sph.cpp) calls `initialize(..., "file")`, which reads
`000.data`. That file must contain exactly `N1 + N2` rows (5500 + 5500 = 11000) in
the same `x, y, z, vx, vy, vz, phase` format as the output frames. The reader turns
commas into whitespace, so the spacing after the commas does not matter, and blank
lines are skipped.

Other options for the string passed to `initialize`: `"cube"`, `"sphere"`, `"layered"`.

To make a new `000.data`:

1. Set `drop_shape` to `"none"` (the file reader expects no drop particles).
2. Run a settling pass, for example with `"cube"` or `"layered"` in `main()`.
3. Copy a settled frame to `000.data`.

Set `drop_shape` back to `"sphere"` or `"cube"` to release a drop above the pool
(size set by `drop_radius_cells`, position by `dropX/dropY/dropZ`).

## Parameters you will most likely change (top of [sph.cpp](sph.cpp))

| Parameter | Meaning |
|---|---|
| `TS`, `dt` | number of steps and time step (CFL: `dt <= 0.4*h/cvel`) |
| `N1`, `N2` | number of water and oil particles |
| `drop_shape`, `drop_radius_cells`, `dropX/Y/Z` | the drop released above the pool |
| `cvelwater`, `cveloil` | numerical speed of sound (stiffness of the equation of state) |
| `av_alpha`, `av_beta`, `av_eps` | Monaghan artificial viscosity |
| `friction` | global damping (0 = off) |
| `muwater`, `muoil` | physical viscosity |
| `thresholdwater`, `thresholdoil` | surface-tension detection thresholds |

Recompile after any change.

## Visualisation

### Animations with PyVista (correct depth ordering)

[render_pyvista.py](render_pyvista.py) renders a folder of frames straight to an mp4
using a real z-buffer. It uses the virtual environment in `.venv/`. `make animate`
calls it for you.

```bash
.venv/bin/python render_pyvista.py part_1
.venv/bin/python render_pyvista.py part_1 --out part1.mp4
.venv/bin/python render_pyvista.py part_2 --opacity 0.5 --point-size 10
.venv/bin/python render_pyvista.py part_3 --rotate --azimuth 45 --elevation 20
.venv/bin/python render_pyvista.py part_1 --limit 50 --out test.mp4   # quick test
```

Run `.venv/bin/python render_pyvista.py -h` for all options. The default output is
`<folder>/animation.mp4`.

If `.venv/` is missing, recreate it:

```bash
python3 -m venv .venv
.venv/bin/pip install pyvista imageio imageio-ffmpeg
```

## Folder layout

| Path | Contents |
|---|---|
| [sph.cpp](sph.cpp), [kernel.h](kernel.h) | the solver and the SPH kernels |
| [Makefile](Makefile) | build |
| [render_pyvista.py](render_pyvista.py), `.venv/` | PyVista animation script and its environment |
| `000.data` | initial particle state (11000 rows) |

