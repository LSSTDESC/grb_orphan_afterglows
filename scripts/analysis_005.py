import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import afterglowpy as grb
import math 
import pandas as pd
from astropy.cosmology import Planck18 as cosmo
from scipy.integrate import quad
from tqdm.notebook import tqdm


def ObsTime(t, mag):
    index = []
    for i in range(len(mag)-1):
        if mag[i] < 24.5:
            index.append(i)
    if not index:
        return 0
    else:
        dt = (t[max(index)] - t[min(index)])*grb.sec2day
        return dt
        
        
file = open('../data/all_dico_results_PL_005_100000.pkl', 'wb')

file_open = open('../data/all_dico_005_100000.pkl', 'rb')
configs_open = pickle.load(file_open)

Z_results = {'jetType' :      -1,          # Jet Type
             'config' :       {},          # Dictionary Z
             't_obs' :        0,           # Observability duration in days
             'mag_min' :      0,           # Minimum magnitude of the afterglow
             'axis' :         'on',        # On-axis or Off-axis
             'observable' :   'no',        # Observability of the afterglow : no or yes for >4. or >7. in days
             'F_gamma' :      '0'}         # Flux of photons above 100 keV in ph/cm2/s


all_Z_results = []
dt = []


for i in range(len(configs_open)):
    
    t = np.geomspace(1.0e2, 1.0e9, 300)
    nu = 5.0e14
    
    Z = configs_open[i]
    Z['jetType'] = grb.jet.PowerLaw
    Z['thetaCore'] = 0.05

    mag = -2.5 * np.log10(grb.fluxDensity(t, nu, **Z)*1.0e-26) - 48.6
    dt.append(ObsTime(t, mag))
    
    # Filling of Z_results dictionary
    if Z['jetType'] == grb.jet.TopHat:
        Z_results['jetType'] = 'TopHat'
    elif Z['jetType'] == grb.jet.Gaussian:
        Z_results['jetType'] = 'Gaussian'
    elif Z['jetType'] == grb.jet.PowerLaw:
        Z_results['jetType'] = 'PowerLaw'
        
    Z_results['config'] = Z
    Z_results['t_obs'] = dt[i]
    Z_results['mag_min'] = min(mag)
    
    if Z['thetaWing'] < Z['thetaObs']:
        Z_results['axis'] = 'off'
    else:
        Z_results['axis'] = 'on'
    
    if dt[i] < 4.:
        Z_results['observable'] = 'no'
    elif dt[i] >= 4.:
        Z_results['observable'] = '> 4 days'
        if dt[i] >= 7.:
            Z_results['observable'] = '> 7 days'
            
    imax = np.where(mag == min(mag))
    t_max = t[imax]

    if len(t_max) == 1:
        # Range of LAT
        nua = 1.0e8 * 1.6e-19 / 6.62e-34    # 10^5 keV
        nub = 5.0e10 * 1.6e-19 / 6.62e-34    # 5 x 10^7 keV


        nu = np.geomspace(nua, nub, num=100)
        Fnu_max = grb.fluxDensity(t_max, nu, **Z)   # Flux en mJy

        delta_nu = nub - nua

        Nnu = []     # Flux calculated with wavelength
        for i in range(len(nu)):
            Nnu.append(Fnu_max[i] * 10**-29 / (nu[i] * 6.62e-34))

        Z_results['F_gamma'] = sum(Nnu)  * delta_nu / 10**4 # Flux in photon/cm2/s
        
    else:
        Z_results['F_gamma'] = 0
    
    all_Z_results.append(Z_results.copy())
    

pickle.dump(all_Z_results, file)
     
file_open.close()
file.close()
