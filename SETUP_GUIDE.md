# Setup Guide for GRB Orphan Afterglows Project

## Prerequisites

This project requires:
- Conda environment with Python 3.11
- Rubin LSST simulation data
- Dust maps data
- Special packages like `afterglowpy` and `rubin_sim`

## Environment Setup

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

## Required Environment Variables

The simulation scripts require these environment variables to be set:

```bash
export SIMU=/home/bregeon/Work/orphans/grb_orphan_afterglows/simulation/data
export OBS=/home/bregeon/Work/orphans/grb_orphan_afterglows/observation/data
export RUBIN_SIM_DATA=/home/bregeon/Work/orphans/rubin_sim/rubin_sim_data
export DUSTMAPS=/home/bregeon/Work/orphans/dustmaps
```

## Running Tests

To run the unit tests:
```bash
source /home/bregeon/Soft/conda/etc/profile.d/conda.sh && conda activate orphans && python -m pytest tests/ -v
```

## Documentation

To build the documentation:
```bash
cd docs && make html
```

## Key Scripts

- `scripts/generate_grb_pop.py` - Main simulation script
- `scripts/script_orphans.sh` - SLURM batch script for running simulations
- `scripts/launch.py` - Script for launching multiple jobs
