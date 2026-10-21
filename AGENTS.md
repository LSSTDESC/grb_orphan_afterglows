# GRB Orphan Afterglows Project - Context Summary

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
- Run tests: `python -m pytest tests/`
- Run linting: `./pylint.sh`
- Generate docs: `cd docs && make html`

## Data Files
- `sed_Fnu_grb_0.txt`: Example GRB spectrum energy density
- `Fermi_LAT_GRB_table_21042022.fits`: Fermi LAT GRB table
- `Swift_GRB_Table_04012022.txt`: Swift GRB table