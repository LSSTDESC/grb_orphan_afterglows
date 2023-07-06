import afterglowpy as grb
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
import math
import sys
import os
from itertools import chain
import math
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.cosmology import Planck18 as cosmo

import warnings
warnings.filterwarnings('ignore')

#sys.path.append('/pbs/home/m/mmasson/lsst/orphans/orphans')

#import pickling as pkg
import orphans.pickling as pkg
#from tools import ObsTime, mag_to_flux, flux_to_mag, pseudo_obs_with_points, galactic_extinction, get_wl_and_nu_band
from orphans.tools import ObsTime, mag_to_flux, flux_to_mag, pseudo_obs_with_points, galactic_extinction, get_wl_and_nu_band


# Where to save the output files
simu = os.environ['$SIMU']
obs = os.environ['$OBS']

# Job ID
job_id = os.get_env('SLURM_JOB_ID')

# Where are saved the rubin_sim data and dust factor table
rubin_sim_data = os.environ['$RUBIN_SIM_DATA']
path_dustmaps = os.environ['$DUSTMAPS']


# Parameters chosen for the GRB jet
N = 1000           # number of simulated GRBs
jetType = 'PL'     # type of structured jet


# Generate all the configurations
print(f'Generating {N} configurations...')
name = f'{simu}/short_configs_{job_id}'
#name = f'/home/masson/orphans/data/simulations/long_configs'
pkg.generate_configs(N, popType='realistic', grbType='short', filename=name)


# Compute the light curve and calculate some results from it (observability, minimal magnitude...)
print('Calculating results...')
t = np.geomspace(1.0e2, 1.0e9, 300)   # time for which the afterglow light curve is calculated

#name_config = f'/home/masson/orphans/data/simulations/long_configs_{N}'
name_config = f'{simu}/short_configs_{job_id}'
#name_simu = f'/home/masson/orphans/data/simulations/long_simulations_{jetType}_{N}'
name_simu = f'{simu}/short_simulations_{jetType}_{job_id}'

pkg.calculate_results(N, t, filename_in=name_config, filename_out=name_simu)


# Open the results in a Pandas Dataframe
results = pkg.open_results(N, filename=name_simu)

# Number of afterglows observed off-axis for at least 7 days
observable_oa = results[(results['axis'] == 'off') & (results['t_obs'] > 7.)]


# Generate pseudo-observations
print(f'Generating pseudo-observations for {len(observable_oa)} orphans')
#name_pseudo_obs = f'../data/pseudo_obs/long_pseudo_obs_{jetType}_{N}'
name_pseudo_obs = f'{obs}/long_pseudo_obs_{jetType}_{job_id}'
pkg.generate_pseudo_obs(N, path_data=rubin_sim_data, path_dustmaps=path_dustmaps,
                        filename_in=name_simu, filename_out=name_pseudo_obs)
print('Done!')