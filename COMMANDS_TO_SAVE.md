# Important Commands for GRB Orphan Afterglows Project

## Environment Setup
```bash
conda env create -f environment.yml && conda activate orphans
pip install -r requirements.txt
pip install -e .
```

## Running Tests
```bash
python -m pytest tests/
```

## Building Documentation
```bash
cd docs && make html
```

## Running Simulations
```bash
# Set environment variables
export SIMU=/home/bregeon/Work/orphans/grb_orphan_afterglows/simulation/data
export OBS=/home/bregeon/Work/orphans/grb_orphan_afterglows/observation/data
export RUBIN_SIM_DATA=/home/bregeon/Work/orphans/rubin_sim/rubin_sim_data
export DUSTMAPS=/home/bregeon/Work/orphans/dustmaps

# Run simulation script
python scripts/generate_grb_pop.py
```

## Linting
```bash
./pylint.sh
```