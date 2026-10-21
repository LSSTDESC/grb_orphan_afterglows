# GRB Orphan Afterglows Project

Study Gamma-ray Bursts orphan afterglows with Rubin LSST and FINK

Read documentation on ReadTheDocs: [latest](https://grb-orphan-afterglows.readthedocs.io/en/latest/)

## Project Overview

This is a scientific project studying Gamma-ray Burst (GRB) orphan afterglows using Rubin LSST and FINK data. The goal is to model and analyze the optical light curves of orphan afterglows.

## Key Components

### Core Package Structure
- `orphans.grb_interface` – thin wrapper around **afterglowpy** to compute light curves and spectra.
- `orphans.grb_configs` – default parameter dictionary (`GRB_BASE_PARAMS`) used by the interface.
- `orphans.tools` – utility functions for flux/magnitude conversion, galactic extinction, and pseudo-observation generation.
- `orphans.plotting_lc` – helpers for visualising light curves.
- `orphans.correlations` – functions that compute statistical correlations between simulation parameters and observable features, returning pandas DataFrames and optional seaborn heatmaps.
- `orphans.pickling` – I/O helpers to serialize/deserialize simulation results.
- `orphans.tools_rubin_sim` – thin glue to the LSST **rubin_sim** package for sky coordinate handling.

### Scripts
- `scripts/generate_grb_pop.py` – orchestrates configuration generation, light-curve computation, and pseudo-observation creation
- `scripts/script_orphans.sh` – SLURM batch script for running simulations
- `scripts/launch.py` – script for launching multiple batch jobs

### Data Structure
- Simulation data in `data/simulations/`
- Pseudo-observations in `data/pseudo_obs/`
- Dust maps in `data/` (SFD_dust_4096_ngp.fits, SFD_dust_4096_sgp.fits)
- Other data files in `data/` directory

## Setup Requirements

This project requires:
- Conda environment with Python 3.11
- Rubin LSST simulation data
- Dust maps data
- Special packages like `afterglowpy` and `rubin_sim`

### Environment Setup

1. Create the conda environment:
```bash
conda env create -f environment.yml
conda activate orphans
```

2. Install Python dependencies:
```bash
pip install -r requirements.txt
```

3. Install the package in development mode:
```bash
pip install -e .
```

4. Install optional dependencies (required for full functionality):
   - Rubin Sim: see the "Optional dependencies" section below.
   - jetsimpy: `git clone https://github.com/haowang-astro/jetsimpy && cd jetsimpy && pip install -e .`

### Required Environment Variables

The simulation scripts require these environment variables to be set:

```bash
export SIMU=/home/bregeon/Work/orphans/grb_orphan_afterglows/simulation/data
export OBS=/home/bregeon/Work/orphans/grb_orphan_afterglows/observation/data
export RUBIN_SIM_DATA=/home/bregeon/Work/orphans/rubin_sim/rubin_sim_data
export DUSTMAPS=/home/bregeon/Work/orphans/dustmaps
```

## Optional Dependencies

### rubin_sim
Required for sky coordinate handling. Install via:
```bash
git clone https://github.com/lsst/rubin_sim.git
cd rubin_sim
conda install -c conda-forge --file=requirements.txt
pip install -e .
```
Then set the data directory:
```bash
export RUBIN_SIM_DATA_DIR=$MYWORKDIR/rubin_sim/rubin_sim_data
```

### jetsimpy
Used for jet simulations. Install with instructions shown in the README.

## Installation from Source

This is the recommended installation procedure in a conda environment.
It shall take care of the `afterglowpy` dependency cleanly.

```bash
> git clone git@github.com:LSSTDESC/grb_orphan_afterglows.git
or
> https://github.com/LSSTDESC/grb_orphan_afterglows.git
then
> cd orphans
> conda env create -f environment.yml
> conda activate orphans
> pip install -r requirements.txt
> conda install markupsafe==2.0.1
> pip install -e .
```

then if needed, one shall also install the `rubin_sim` package as follows:
```bash
# from the orphans directory
> cd ..
> git clone https://github.com/lsst/rubin_sim.git
> cd rubin_sim
> conda install -c conda-forge --file=requirements.txt
> pip install -e .
> cd ../orphans
```
in which case, one shall follow instructions to download the scheduler data, and then set:
```bash
> export RUBIN_SIM_DATA_DIR=$MYWORKDIR/rubin_sim/rubin_sim_data"
```

## Development Commands

### Common Development Commands
- **Install the package (editable)**: `pip install -e .`
- **Create the conda environment**: `conda env create -f environment.yml && conda activate orphans`
- **Run the test suite with coverage**: `python -m pytest --cov=./orphans --cov-report=html tests`
- **Run a single test**: `python -m pytest tests/<test_file>.py::TestClass::test_method`
- **Run a subset of tests** (by pattern): `python -m pytest -k "<keyword>" tests`
- **Lint the code**: `./pylint.sh`
- **Generate documentation**: `cd docs && make html`
- **Run the main simulation script** (requires environment variables `SIMU`, `OBS`, `RUBIN_SIM_DATA`, `DUSTMAPS`): `python scripts/generate_grb_pop.py`
- **Run a notebook**: launch with `jupyter lab` and open any `.ipynb` in `notebooks/`

## Running Tests

To run the unit tests:
```bash
source /home/bregeon/Soft/conda/etc/profile.d/conda.sh && conda activate orphans && python -m pytest tests/ -v
```

## Data Folder Content

- `sed_Fnu_grb_0.txt` : simple text file with an example of a spectrum energy density of a GRB (lambda, Fnu)
                      produced with `afterglowpy`
- `Fermi_LAT_GRB_table_21042022.fits` : FITS file with a table of Fermi LAT GRBs, analyzed as for the 2nd Fermi--LAT
                                      GRB catalog ([2019ApJ...878...52A](https://ui.adsabs.harvard.edu/abs/2019ApJ...878...52A/abstract))
  - https://heasarc.gsfc.nasa.gov/W3Browse/fermi/fermilgrb.html
  - https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3query.pl#fermilgrb
- `Swift_GRB_Table_04012022.txt` : text file with a table of SWIFT GRBs, values are tab separated
  - https://swift.gsfc.nasa.gov/archive/grb_table/

## People
- Marina Masson <marina.masson@lpsc.in2p3.fr>
- [Johan Bregeon](https://github.com/bregeon) (CNRS IN2P3 LPSC) <bregeon@in2p3.fr>

## License, Contributing etc

The code in this repository is available for re-use under the GPL v3 license.

## Important Non-code Notes

- The repository relies on external data directories (`RUBIN_SIM_DATA`, `DUSTMAPS`) and must have the `SIMU` and `OBS` paths set before running the simulation scripts.
- The `setup.py` uses `setuptools_scm` to version the package from git tags.
- Linting uses `pylint-3` targeting both the library and the `scripts/` directory.
- CI can be run locally by executing the same steps defined in `.gitlab-ci.yml` after setting up the conda environment.

## Testing Strategy

Run tests using:
```bash
python -m pytest tests/
```

The failing tests are in the correlations module and need to be debugged and fixed to ensure proper functionality.