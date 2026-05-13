# -*- coding: utf-8 -*-
"""Generate GRB population simulations and pseudo‑observations.

This script creates GRB configurations, computes light‑curves, and produces
pseudo‑observations using the utilities in :mod:`orphans.pickling`.
Environment variables required:
- ``SIMU``: directory for simulation outputs
- ``OBS``: directory for observation outputs
- ``SLURM_JOB_ID``: identifier for the current job (used in file names)
- ``RUBIN_SIM_DATA`` and ``DUSTMAPS``: paths needed for pseudo‑observation
  generation.
"""

import os
import warnings
import numpy as np
import orphans.pickling as pkg

warnings.filterwarnings("ignore")

# Required environment variables
SIMU = os.environ["SIMU"]
OBS = os.environ["OBS"]
RUBIN_SIM_DATA = os.environ["RUBIN_SIM_DATA"]
DUSTMAPS = os.environ["DUSTMAPS"]

# Job identifier (may be ``None`` if not set)
JOB_ID = os.getenv("SLURM_JOB_ID")

# Simulation parameters
N = 1000  # number of simulated GRBs
JET_TYPE = "PL"  # type of structured jet

print(f"Generating {N} configurations...")
CONFIG_NAME = f"{SIMU}/short_configs_{JOB_ID}"
# Generate configurations
pkg.generate_configs(N, popType="realistic", grbType="short", filename=CONFIG_NAME)

print("Calculating results...")
time_grid = np.geomspace(1.0e2, 1.0e9, 300)
config_path = CONFIG_NAME
SIM_NAME = f"{SIMU}/short_simulations_{JET_TYPE}_{JOB_ID}"
# Compute light curves and results
pkg.calculate_results(N, time_grid, filename_in=config_path, filename_out=SIM_NAME)

# Load results into a DataFrame
results = pkg.open_results(N, filename=SIM_NAME)

# Select observable off‑axis afterglows (observed > 7 days)
observable_oa = results[(results["axis"] == "off") & (results["t_obs"] > 7.0)]

print(f"Generating pseudo‑observations for {len(observable_oa)} orphans")
pseudo_obs_name = f"{OBS}/long_pseudo_obs_{JET_TYPE}_{JOB_ID}"
# Generate pseudo‑observations
pkg.generate_pseudo_obs(
    N,
    path_data=RUBIN_SIM_DATA,
    path_dustmaps=DUSTMAPS,
    filename_in=SIM_NAME,
    filename_out=pseudo_obs_name,
)

print("Done!")
