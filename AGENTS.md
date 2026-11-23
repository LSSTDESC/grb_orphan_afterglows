# GRB Orphan Afterglows Project - Context Summary

## Agents instructions: Conciseness & Precision Guidelines
- Keep explanations brief; focus on core points.                                                                 
- Use specific, concrete terms—avoid vague language.
- Prioritize direct answers over elaboration unless the user explicitly requests examples or deeper detail.      

## Project Purpose
Scientific project studying Gamma-ray Burst (GRB) orphan afterglows using Rubin LSST and FINK data to model and analyze optical light curves.

## Core Package Structure
- `orphans.grb_interface`: Wrapper around afterglowpy for light curve and spectrum computation
- `orphans.grb_configs`: Default parameter dictionary (GRB_BASE_PARAMS) for the interface
- `orphans.tools`: Utility functions for flux/magnitude conversion, galactic extinction, and pseudo-observation generation
- `orphans.plotting_lc`: Light curve visualization helpers
- `orphans.correlations`: Statistical correlation computations between simulation parameters and observables
- `orphans.pickling`: Serialization/deserialization helpers for simulation results
- `orphans.tools_rubin_sim`: Integration glue for LSST rubin_sim package

## Key Scripts
- `scripts/generate_grb_pop.py`: Orchestrates configuration generation, light-curve computation, and pseudo-observation creation
- `scripts/script_orphans.sh`: SLURM batch script for running simulations
- `scripts/launch.py`: Script for launching multiple batch jobs

## Data Structure
- Simulation data: `data/simulations/`
- Pseudo-observations: `data/pseudo_obs/`
- Dust maps: `data/SFD_dust_4096_ngp.fits`, `data/SFD_dust_4096_sgp.fits`
- Other data files: `data/` directory

## Environment Requirements
- Conda environment with Python 3.11
- Rubin LSST simulation data
- Dust maps data
- Special packages: afterglowpy, rubin_sim

## Required Environment Variables
- `SIMU`: Simulation data directory
- `OBS`: Observation data directory  
- `RUBIN_SIM_DATA`: Rubin simulation data directory
- `DUSTMAPS`: Dust maps directory

## Installation
1. Create conda environment: `conda env create -f environment.yml && conda activate orphans`
2. Install dependencies: `pip install -r requirements.txt`
3. Install package: `pip install -e .`
4. Install optional dependencies: rubin_sim and jetsimpy

## Development Commands
- Install package: `pip install -e .`
- Run tests: `eval "$(/home/bregeon/Soft/conda/bin/conda shell.bash hook)"; conda activate orphans; python -m pytest --cov=orphans --cov-report=html tests`
  - Note: the `conda activate` must be on the same line as pytest, as it does not persist across separate shells.
- Run linting: `./pylint.sh`
- Generate docs: `cd docs && make html`

## Virtual Environment (uv)
A uv-managed Python 3.11 environment is used for reproducibility (`uv.lock`, 199 packages).
- Create/activate: `uv venv -p 3.11 .venv` (add `.venv` to `.gitignore`).
- Install package + deps: `uv pip install -e .`
- Optional extras: `uv pip install --extra rubin` (rubin_sim).
  - `jetsimpy` is a **required** dependency installed from git at build time (pybind11 C++ extension, needs gcc); pinned in `uv.lock` (see `[project.dependencies]`).
- Dev tooling: `uv pip install --group dev` (pytest, pytest-cov, pylint, sphinx, jupyterlab).
- Run tests: `uv run --python .venv/bin/python python -m pytest --cov=orphans --cov-report=html tests`
- Refresh lockfile: `uv lock`
- h5py: resolves via a prebuilt PyPI wheel (3.16.0, cp311 manylinux) — no HDF5 source build needed. If a future h5py has no wheel for the target interpreter, build it against conda HDF5 headers+libs (the original `environment.yml` conda setup handles this), so keep that as a fallback.
- Note: `pytest`/`--cov` are dev-only, not core deps; use `uv run` so the resolved interpreter is picked up.

## Git Workflow
- Always `git branch --show-current` before committing/pushing.
- `git push` with no branch argument pushes the **current** branch (currently `main`). Specify the branch explicitly to avoid landing on `main`.

## Data Files
- `sed_Fnu_grb_0.txt`: Example GRB spectrum energy density
- `Fermi_LAT_GRB_table_21042022.fits`: Fermi LAT GRB table
- `Swift_GRB_Table_04012022.txt`: Swift GRB table
