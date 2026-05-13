import os
import sys
import numpy as np
import jetsimpy
import pickle
import pandas as pd
import sqlite3
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D

from astropy.cosmology import Planck18 as cosmo
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord

sys.path.append('/home/masson/rubin_sim/')

from rubin_sim.photUtils.Bandpass import Bandpass
from rubin_sim.photUtils.Sed import Sed
from rubin_sim.data import get_baseline

from modules.tools import obs_duration_th, time_coord, flux_to_mag
from modules.functions_pseudo_obs import compute_mags, df_obs, real_obs, grb_obs_duration

import warnings
warnings.filterwarnings('ignore')


file_open = open('data/configs_jetsimpy.pkl', 'rb')
configs_oa = pickle.load(file_open)
file_open.close()


baseline_db = get_baseline()
print(baseline_db)

conn = sqlite3.connect(baseline_db)

# In the near future, 'summaryallprops' will be replaced with 'observations'
df = pd.read_sql('select * from observations;', conn)

conn.close()

path_rubin_sim_data = '/home/masson/rubin_sim_data'
fdir = os.path.join(path_rubin_sim_data, 'throughputs', 'baseline')

# Read the throughput curves
filterlist = ['u', 'g', 'r', 'i', 'z', 'y']
filtercolors = {'u': 'b', 'g': 'c', 'r': 'g', 'i': 'orange', 'z': 'r', 'y': 'm'}

lsst = {}
for f in filterlist:
    lsst[f] = Bandpass()
    lsst[f].readThroughput(os.path.join(fdir, f'total_{f}.dat'))

all_LC = []

#for i in range(len(configs_oa)):
for i in range(2000, 3000):

    if i%50 == 0:
        print(f'{i} orphans done')

    # dictionary containing all the GRB information
    LC = {'config': {},  # Dictionary P
          'grb_time': 0,  # GRB observation date
          'grb_coord': 0,  # GRB ra/dec coordinates
          'time': [],  # Time of each detection
          'mags': [],  # Magnitude of the detection
          'filt': [],  # Filter used at the moment of the detection
          'mags_lim': [],  # Limiting magnitude of the detection at the observation time
          'mags_err': []}  # Error on the magnitude

    config = configs_oa[i]['config']

    # Time in mjd and ra/dec coordinates of the GRB
    grb_time, grb_coord = time_coord()

    # Observe t_before days before and t_after days after
    t_before = TimeDelta(20, format='jd')
    t_after = TimeDelta(365, format='jd')
    obs_start = grb_time - t_before
    obs_end = grb_time + t_after

    # Get time span
    df_time = df[(df['observationStartMJD'] > obs_start.mjd) & (df['observationStartMJD'] < obs_end.mjd)]

    # Angular separation with SkyCoord.separation
    # Rubin FOV is 47 square degree for a 3.5-degree diameter, hence 1.7 deg separation radius.
    df_time['Separation'] = SkyCoord(df_time['fieldRA'], df_time['fieldDec'], unit="deg").separation(grb_coord).degree
    df_sky = df_time[df_time['Separation'] < 1.7]

    time_bins = grb_obs_duration(df_sky, grb_time)

    obs_list = df_obs(config, df_sky, time_bins, lsst)

    y_mags = []  # magnitude with extinction

    # If there is no observation, let's go to the next configuration
    if len(obs_list) == 0:
        all_LC.append(0)
        continue

    else:
        obs_df = pd.concat(obs_list)

        obs_df['observationId'] = df_sky['observationId']

        x_times, y_mags_without_ext, z_colors, mags_lim, mags_err = real_obs(obs_df, df_sky, time_bins, grb_time, lsst)

        LC['mags'] = np.array(y_mags_without_ext)
        LC['config'] = config
        LC['grb_time'] = grb_time.isot
        LC['grb_coord'] = grb_coord.to_string('hmsdms')
        LC['time'] = np.array(x_times)
        LC['filt'] = np.array(z_colors)
        LC['mags_lim'] = np.array(mags_lim)
        LC['mags_err'] = np.array(mags_err)[:, 0, 0]

        if len(LC['mags'][LC['mags'] < LC['mags_lim']]) > 5:
            all_LC.append(LC.copy())


file = open(f'/home/masson/jetsimpy/orphan_pop/data/pseudo_obs_jetsimpy_5pts_3.pkl', 'wb')
pickle.dump(all_LC, file)
file.close()