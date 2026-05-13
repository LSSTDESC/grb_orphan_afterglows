# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Common development commands
- **Install the package (editable)**: `pip install -e .`
- **Create the conda environment**: `conda env create -f environment.yml && conda activate orphans`
- **Install optional dependencies** (required for full functionality):
  - Rubin Sim: see the "Optional dependencies" section below.
  - jetsimpy: `git clone https://github.com/haowang-astro/jetsimpy && cd jetsimpy && pip install -e .`
- **Run the test suite with coverage**: `python -m pytest --cov=./orphans --cov-report=html tests`
- **Run a single test**: `python -m pytest tests/<test_file>.py::TestClass::test_method`
- **Run a subset of tests** (by pattern): `python -m pytest -k "<keyword>" tests`
- **Lint the code**: `./pylint.sh`
- **Generate documentation**: `cd docs && make html`
- **Run the main simulation script** (requires environment variables `SIMU`, `OBS`, `RUBIN_SIM_DATA`, `DUSTMAPS`): `python scripts/generate_grb_pop.py`
- **Run a notebook**: launch with `jupyter lab` and open any `.ipynb` in `notebooks/`

## High‑level architecture
- The *orphans* package is the core library, exposing the following sub‑modules:
  - `orphans.grb_interface` – thin wrapper around **afterglowpy** to compute light curves and spectra.
  - `orphans.grb_configs` – default parameter dictionary (`GRB_BASE_PARAMS`) used by the interface.
  - `orphans.tools` – utility functions for flux/magnitude conversion, galactic extinction, and pseudo‑observation generation.
  - `orphans.plotting_lc` – helpers for visualising light curves.
  - `orphans.correlations` – functions that compute statistical correlations between simulation parameters and observable features, returning pandas DataFrames and optional seaborn heatmaps.
  - `orphans.pickling` – I/O helpers to serialize/deserialize simulation results.
  - `orphans.tools_rubin_sim` – thin glue to the LSST **rubin_sim** package for sky coordinate handling.
- **Scripts** in `scripts/` provide entry‑points for specific workflows:
  - `generate_grb_pop.py` – orchestrates configuration generation, light‑curve computation, and pseudo‑observation creation; it expects the environment variables listed above.
  - `agpy_test.py` – quick sanity‑check visualisation of afterglow light curves.
- **Notebooks** under `notebooks/` demonstrate analysis pipelines (e.g., pseudo‑observation fitting, correlation heatmaps) and are useful for exploratory work.
- CI definition (`.gitlab-ci.yml`) builds the conda environment, installs the package, runs tests, and builds the Sphinx docs.

## Optional dependencies
- **rubin_sim** – required for sky coordinate handling. Install via:
  ```
  git clone https://github.com/lsst/rubin_sim.git
  cd rubin_sim
  conda install -c conda-forge --file=requirements.txt
  pip install -e .
  ```
  Then set the data directory:
  ```
  export RUBIN_SIM_DATA_DIR=$MYWORKDIR/rubin_sim/rubin_sim_data
  ```
- **jetsimpy** – used for jet simulations. Install with the commands shown in the README.

## Important non‑code notes
- The repository relies on external data directories (`RUBIN_SIM_DATA`, `DUSTMAPS`) and must have the `SIMU` and `OBS` paths set before running the simulation scripts.
- The `setup.py` uses `setuptools_scm` to version the package from git tags.
- Linting uses `pylint-3` targeting both the library and the `scripts/` directory.
- CI can be run locally by executing the same steps defined in `.gitlab-ci.yml` after setting up the conda environment.
