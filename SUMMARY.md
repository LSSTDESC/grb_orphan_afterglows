# GRB Orphan Afterglows Project - Summary

## Project Overview
This is a scientific project studying Gamma-ray Burst (GRB) orphan afterglows using Rubin LSST and FINK data. The goal is to model and analyze the optical light curves of orphan afterglows.

## Key Components

### 1. Main Simulation Scripts
- `scripts/generate_grb_pop.py` - Main script for generating GRB population and simulations
- `scripts/script_orphans.sh` - SLURM batch script for running simulations
- `scripts/launch.py` - Script for launching multiple batch jobs

### 2. Data Structure
- Simulation data in `data/simulations/`
- Pseudo-observations in `data/pseudo_obs/`
- Dust maps in `data/` (SFD_dust_4096_ngp.fits, SFD_dust_4096_sgp.fits)
- Other data files in `data/` directory

### 3. Core Dependencies
- afterglowpy - GRB afterglow modeling
- rubin_sim - Rubin LSST simulation package
- dustmaps - Dust extinction calculations
- Various scientific Python libraries (numpy, scipy, pandas, etc.)

## Setup Requirements
The project requires:
1. Conda environment with Python 3.11
2. Environment variables (SIMU, OBS, RUBIN_SIM_DATA, DUSTMAPS) 
3. Installation of rubin_sim and jetsimpy packages
4. Special handling for dustmaps and other scientific libraries

## Current Status
- Unit tests are mostly passing (20/24 tests)
- 4 tests in the correlations module are currently failing and need investigation
- The simulation workflow is designed to run in a cluster environment with SLURM

## Testing Strategy
Run tests using:
```bash
python -m pytest tests/
```

The failing tests are in the correlations module and need to be debugged and fixed to ensure proper functionality.