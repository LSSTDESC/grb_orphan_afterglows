"""Tools of rubin_sim module

This module provides a number of generic tools to use with the PicklePseudoObs function from the pickling module.

"""


import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
import pickle
import afterglowpy as grb

import os
import rubin_sim.phot_utils.bandpass as Bandpass
import rubin_sim.phot_utils.sed as Sed
from rubin_sim.data import get_baseline

from grb_interface import make_grb_spectrum, dump_wl_Fnu_spectrum




def ComputeMags(i, wls, fnus, obs_t, f, lsst):

    """ Compute magnitudes
    """
    
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




def GRBObsTime(df_sky, grb_time):
    """ Compute the observation time in the GRB time frame, i.e. from GRB T_0
    """

    obs_times_grb_frame = df_sky['observationStartMJD'] - grb_time.mjd
    time_bins = obs_times_grb_frame[obs_times_grb_frame > 0]
    return time_bins
    # time_bins.to_dict()
    
    
    
    
def df_obs(Z, df_sky, time_bins, lsst):

    """ Compute magnitudes in each filter at the observation time
    """
    
    Fnu_Jy = dict()
    for obs_id, obs_t in time_bins.to_dict().items():
        wl_full_band, freq_full_band, t, Fnu_Jy[obs_id] = make_grb_spectrum(E0=Z['E0'], z=Z['z'], n0=Z['n0'], thetaObs=Z['thetaObs'], thetaCore=Z['thetaCore'], thetaWing=Z['thetaWing'], t=obs_t * grb.day2sec)
    
    obs_list = list()
    for obs_id, fnu_val in Fnu_Jy.items():
        obs_t = time_bins[obs_id]
        filt = df_sky[df_sky['observationId']==obs_id]['filter']
        df = compute_mags(obs_id, wl_full_band, fnu_val, obs_t, filt.values[0], lsst)
        obs_list.append(df)
        
    return obs_list




def real_obs(obs_df, df_sky, time_bins, grb_time):

    """ Keep only "real" observations for the right filter
    """
    
    filtercolors = {'u':'b', 'g':'c', 'r':'g', 'i':'orange', 'z':'r', 'y':'m'}
    x_times = list()
    y_mags = list()
    z_colors = list()
    mags_lim = list()
    for obs_id, obs in obs_df.iterrows():
        filt = df_sky[df_sky['observationId']==obs_id]['filter']
        lim = df_sky[df_sky['observationId']==obs_id]['fiveSigmaDepth']
        obs_t = time_bins[obs_id]
        
        x_times.append(obs_t + grb_time.mjd)
        y_mags.append(obs[filt].values[0])
        z_colors.append(filtercolors[filt.values[0]])
        mags_lim.append(lim.values[0])
        
    return x_times, y_mags, z_colors, mags_lim

