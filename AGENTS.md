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
export SIMU=/path/to/simulation/data
export OBS=/path/to/observation/data
export RUBIN_SIM_DATA=/path/to/rubin_sim/data
export DUSTMAPS=/path/to/dustmaps/data
```
**Will an agent miss this without help?** Yes – environment variables are not obvious from the code.

## Running Tests
```bash
python -m pytest tests/
```
**Will an agent miss this without help?** No – standard pytest usage.

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