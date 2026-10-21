# AGENTS.md

## Setup
Create the conda environment and install dependencies:
```bash
conda env create -f environment.yml && conda activate orphans
pip install -r requirements.txt
pip install -e .
```
**Will an agent miss this without help?** Yes – the full setup process is not obvious.

## Environment Variables
The main simulation script requires these environment variables:
```bash
export SIMU=/home/bregeon/Work/orphans/grb_orphan_afterglows/simulation/data
export OBS=/home/bregeon/Work/orphans/grb_orphan_afterglows/observation/data
export RUBIN_SIM_DATA=/home/bregeon/Work/orphans/rubin_sim/rubin_sim_data
export DUSTMAPS=/home/bregeon/Work/orphans/dustmaps
```
**Will an agent miss this without help?** Yes – environment variables are not obvious from the code.

## Running Tests
```bash
python -m pytest tests/
```
**Will an agent miss this without help?** No – standard pytest usage.

## Unit Test Results
The following tests are currently failing:
- test_duration_between_first_and_peak
- test_heatmap_structure
- test_minimal_magnitude
- test_peak_time

These failures are in the correlations module and need to be investigated and fixed.

## Linting
```bash
./pylint.sh
```
**Will an agent miss this without help?** No – standard linting command.

## Documentation
```bash
cd docs && make html
```
**Will an agent miss this without help?** Yes – not obvious from README.