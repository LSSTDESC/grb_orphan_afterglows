import pandas as pd
import sqlite3
import sys
from rubin_sim.data import get_baseline
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
import pickle
from tqdm.notebook import tqdm

from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
from astropy import units as u

# Import the primary photometry classes from rubin_sim.photUtils
import os
import rubin_sim.photUtils.Bandpass as Bandpass
import rubin_sim.photUtils.Sed as Sed

import afterglowpy as grb
from grb_interface import make_grb_spectrum, dump_wl_Fnu_spectrum

import warnings
warnings.filterwarnings('ignore')

sys.path.append('../modules')
import resource


# Function that take randomly the time and ra/dec of a GRB during LSST
def time_coord():
    print(f'time_coord {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}')
    grb_datetime = np.random.uniform(59945, 63598) #time in mjd : t_start = 01/01/2023, t_end = 01/01/2033
    grb_time = Time(grb_datetime, format='mjd', scale='utc')
    
    grb_ra = np.random.uniform(-180.00,180.00)*u.degree 
    grb_dec = np.random.uniform(-90.00,0.00)*u.degree
    grb_coord = SkyCoord(grb_ra, grb_dec, frame='icrs') #ra/dec coordinates in Southern Hemisphere
    
    return grb_time, grb_coord


# Function that compute magnitudes
def compute_mags(i, wls, fnus, obs_t, f):
    #print(f'compute_mags {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}')
    new_grb_sed = Sed()
    new_grb_sed.wavelen = np.array(wls)
    new_grb_sed.fnu = np.array(fnus)
    # convert fnu to flambda
    new_grb_sed.fnuToflambda()
    # Calculate expected AB magnitudes. 
    new_grb_mags = {}
    new_grb_mags[str(f)] = new_grb_sed.calcMag(lsst[f])
    # time is one column
    new_grb_mags['obs_time'] = obs_t
    # Make a dataframe just to get a nice output cell.
    return pd.DataFrame(new_grb_mags, index=[i])
    
    
# Function that compute the observation time in the GRB time frame, i.e. from GRB T_0
def obs_time(df_s, grb_time):
    print(f'obs_time {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}')
    obs_times_grb_frame = df_s['observationStartMJD'] - grb_time.mjd
    time_bins = obs_times_grb_frame[obs_times_grb_frame>0]
    return time_bins
    #time_bins.to_dict()
    
    
# Function that compute magnitudes in each filter at the observation time
def df_obs(Z, df_s, time_bins):
    print(f'df_obs {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}')
    Fnu_Jy = dict()
    for obs_id, obs_t in time_bins.to_dict().items():
        wl_full_band, freq_full_band, t, Fnu_Jy[obs_id] = make_grb_spectrum(E0=Z['E0'], z=Z['z'], n0=Z['n0'], thetaObs=Z['thetaObs'], thetaCore=Z['thetaCore'], thetaWing=Z['thetaWing'], t=obs_t * grb.day2sec)
    
    obs_list = list()
    for obs_id, fnu_val in Fnu_Jy.items():
        obs_t = time_bins[obs_id]
        filt = df_s[df_s['observationId']==obs_id]['filter']
        df = compute_mags(obs_id, wl_full_band, fnu_val, obs_t, filt.values[0])
        obs_list.append(df)
        
    return obs_list
    
    
# Function that keep only "real" observations for the right filter
def real_obs(obs_df, df_s, time_bins, grb_time):
    print(f'real_obs {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}')
    x_times = list()
    y_mags = list()
    z_colors = list()
    mags_lim = list()
    for obs_id, obs in obs_df.iterrows():
        filt = df_s[df_s['observationId']==obs_id]['filter']
        lim = df_s[df_s['observationId']==obs_id]['fiveSigmaDepth']
        obs_t = time_bins[obs_id]
        
        x_times.append(obs_t + grb_time.mjd)
        y_mags.append(obs[filt].values[0])
        z_colors.append(filtercolors[filt.values[0]])
        mags_lim.append(lim.values[0])
        
    return x_times, y_mags, z_colors, mags_lim
    

if __name__ == "__main__":
    # or replace with any local file you want to use
    baseline_db = get_baseline()
    print(baseline_db)

    conn = sqlite3.connect(baseline_db)

    # In the near future, 'summaryallprops' will be replaced with 'observations'
    df = pd.read_sql('select * from observations;', conn)

    conn.close()

    fdir = '/home/masson/rubin_sim_data'
    fdir = os.path.join(fdir, 'throughputs', 'baseline')

    # Read the throughput curves
    filterlist = ['u', 'g', 'r', 'i', 'z', 'y']
    filtercolors = {'u':'b', 'g':'c', 'r':'g', 'i':'orange', 'z':'r', 'y':'m'}

    lsst = {}
    for f in filterlist:
        lsst[f] = Bandpass()
        lsst[f].readThroughput(os.path.join(fdir, f'total_{f}.dat'))
        
        
    file_open = open('../data/all_dico_results_PL_015_100000.pkl', 'rb')
    configs_open = pickle.load(file_open)
    file_open.close()

    configs = pd.DataFrame(configs_open)
    configs = configs[(configs['axis'] == 'off') & (configs['t_obs'] > 7.)]
    configs = configs[configs['mag_min']<22]

    file = open('../config/lc_configs_PL_100000.pkl', 'wb')



    all_LC = []


    for config in tqdm(configs['config']):
        
        Z = config
        print(Z['z'])
        
        # Time in mjd and ra/dec coordinates of the GRB
        grb_time, grb_coord = time_coord()
        # print(grb_time.isot, grb_coord.to_string('hmsdms'))
        
        #grb_time = 61402.36838867448
        #grb_time = Time(grb_datetime, format='mjd', scale='utc')
        
        #grb_ra = "03h38m30.03464941s"
        #grb_dec = "-54d36m49.59379716s"
        #grb_coord = SkyCoord(grb_ra, grb_dec, frame='icrs')
        
        # Observe t_before days before and t_after days after
        t_before = TimeDelta(20, format='jd')
        t_after = TimeDelta(365, format='jd')
        obs_start = grb_time - t_before
        obs_end = grb_time + t_after
        
        # Get time span
        print(f'Before time span {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}')
        df_time = df[(df['observationStartMJD']>obs_start.mjd) & (df['observationStartMJD']<obs_end.mjd)]
        print(f'After time span {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}')
        # Angular separation with SkyCoord.separation
        # Rubin FOV is 47 square degree for a 3.5-degree diameter, hence 1.7 deg separation radius.

        df_time['Separation'] = SkyCoord(df_time['fieldRA'], df_time['fieldDec'], unit="deg").separation(grb_coord).degree
        df_sky = df_time[df_time['Separation']<1.7]
        print(f'After separation {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}')
        print(f'df_sky size {len(df_sky)}')
        #print(df_sky)
        
        time_bins = obs_time(df_sky, grb_time)
        print(f'after ob time {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}')
        #print(time_bins.to_dict())
        
        obs_list = df_obs(Z, df_sky, time_bins)
        
        # If there is no observation, let's go to the next configuration
        if len(obs_list)==0:
            all_LC.append(0)
            print('no obs')
            continue
            
        else:
            obs_df = pd.concat(obs_list)
            
            obs_df['observationId'] = df_sky['observationId']
            print(f'obs id {obs_df["observationId"].head(1)}')

            x_times, y_mags, z_colors, mags_lim = real_obs(obs_df, df_sky, time_bins, grb_time)
            
            LC = {'config' :       {},        # Dictionary Z
                  'grb_time' :     0,         # GRB observation date
                  'grb_coord' :    0,         # GRB ra/dec coordinates
                  'time' :         [],        # Time of each detection
                  'mags' :         [],        # Magnitude of the detection
                  'filt' :         [],        # Filter used at the moment of the detection
                  'mags_lim' :     []}        # Limiting magnitude of the detection at the observation time

            LC['config'] = Z
            LC['grb_time'] = grb_time.isot
            LC['grb_coord'] = grb_coord.to_string('hmsdms')
            LC['time'] = x_times
            LC['mags'] = y_mags
            LC['filt'] = z_colors
            LC['mags_lim'] = mags_lim
            #l = LC.copy()
            #print(f'dict size: {l.__sizeof__()}')
            #all_LC.append(LC.copy())
            all_LC.append(LC)
            del(obs_df)
            del(df_sky)
            del(df_time)


    pickle.dump(all_LC, file)
    file.close()
