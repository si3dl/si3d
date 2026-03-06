# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

**pSi3D** is a semi-implicit 3D hydrodynamic model for lake/reservoir simulation, written in Fortran 90 and parallelized with OpenMP. It solves the Navier-Stokes equations using a leapfrog-trapezoidal finite-difference scheme on a structured grid.

## Build Commands

Requires Intel Fortran Compiler (`ifort`) and OpenMP (from Intel oneAPI HPC Toolkit).

```bash
# Standard build (OpenMP, no GOTM turbulence library)
make omp

# Or via the shell script (sets paths, then calls make omp)
bash ./Compiler.sh
```

The output binary is `psi3d` in the working directory.

**Makefile targets:**
- `omp` — standard build with `ifort` + OpenMP (most commonly used)
- `all` — MPI + OpenMP build with `mpif90`
- `sec` — `ifort` + OpenMP + O2 optimization
- `mpi` — MPI-only build

The `Compiler.sh` script sets `GOTMDIR`, `MODDIR`, `INCDIR`, `LIBDIR` for the optional GOTM turbulence library. With `IFGOTM=false` (default), only si3d is compiled. The `si3d_mixing.f90` module has `USE turbulence` / `USE kpp` imports — these come from GOTM and must be available even if GOTM is not rebuilt.

## Running the Model

```bash
# Copy binary to run directory and execute
cp ./psi3d ./sample_runs/windDrivenFlow/
cd ./sample_runs/windDrivenFlow/
./psi3d
```

The model reads input from the current working directory. All input files must be present at runtime.

## Module Architecture

The compilation order matters (modules must be compiled before their consumers):

```
nspcg.f                    (sparse solver, external)
PlumeModels.f              (plume models, external)
si3d_types.f90             (global variables, constants, array declarations)
si3d_stwave.f90            (wave modeling module)
si3d_sed.f90               (sediment transport)
si3d_Hg.f90                (mercury biogeochemistry; USEs si3d_sed)
si3d_wq.f90                (water quality: DO, nutrients; USEs si3d_types)
si3d_ecomod.f90            (ecological model orchestration; USEs wq, sed, stwave, Hg)
si3d_utils.f90             (I/O, date utilities, space allocation; USEs ecomod)
si3d_boundaryconditions.f90 (open & surface BCs; USEs utils, ecomod)
si3d_mixing.f90            (vertical turbulence: M-Y 2.5, KPP, constant; USEs GOTM turbulence)
si3d_procedures.f90        (core hydrodynamic subroutines; USEs all above)
si3d.f90                   (main program)
```

**Key dependency**: `si3d_types` defines all shared global variables (using `SAVE`). Every other module `USE`s it directly, meaning the model is stateful through module-level globals rather than explicit argument passing.

## Input Files

All input files are read from the model run directory at runtime:

| File | Description |
|---|---|
| `si3d_inp.txt` | Main parameter file (domain size, time step, physics switches, output config) |
| `si3d_bathy` | Bathymetry grid |
| `si3d_layer.txt` | Vertical layer definition |
| `si3d_init.txt` | Initial conditions (temperature/salinity profiles) |
| `si3d_surfbc.txt` | Surface boundary conditions (wind, heat flux time series) |
| `si3d_wq_inp.txt` | Water quality parameters (when `ecomod > 0`) |
| `openbc0N.txt` | Open boundary condition files, one per boundary (`N` = boundary number) |

See `sample_files/si3d_inp_Info.txt` for annotated descriptions of every parameter in `si3d_inp.txt`.

## Key Parameters in `si3d_inp.txt`

- `iturb`: Turbulence model (`0`=constant, `1`=Mellor-Yamada 2-eq)
- `isal`: Solve scalar transport (`1`=yes)
- `ecomod`: Ecological module switch (`0`=conservative tracers, `1`=WQ active)
- `ntr`: Number of tracers
- `ipwq`: Frequency (time steps) to call WQ module
- `ifsbc`: Surface BC mode (`0`=constant, `1`=preprocessed, `2/3`=runtime)
- `idbg`: Debug verbosity (`0`=none, `1`=checkpoint prints)
- `nth`: Number of OpenMP threads

## Water Quality / Ecological Modules

The WQ stack (active when `ecomod > 0`):
- `si3d_wq.f90` — dissolved oxygen, reaeration, sediment oxygen demand
- `si3d_sed.f90` — suspended sediment settling/resuspension
- `si3d_Hg.f90` — mercury speciation and cycling
- `si3d_stwave.f90` — wave energy (STWAVE-derived)
- `si3d_ecomod.f90` — orchestrates all sub-modules; called from `si3d_procedures` via `sourcesink` array

Source/sink terms are accumulated in the `sourcesink` array (units: flux [M/L²/T]) and applied during tracer transport.

## Active Branch: `si3d_WQ_v2`

Current work is on the water quality module. Modified files relative to `main`: `si3d_procedures.f90`, `Compiler.sh`, `psi3d` binary.
