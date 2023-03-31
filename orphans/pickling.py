"""Pickle module

"""


import pickle
import numpy as np
import afterglowpy as grb
import math 
import pandas as pd
import sqlite3
from astropy.cosmology import Planck18 as cosmo
from tqdm.notebook import tqdm
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
from scipy.stats import rv_continuous
from scipy.integrate import simps
from astropy import units as u

import os
import rubin_sim.photUtils.Bandpass as Bandpass
import rubin_sim.photUtils.Sed as Sed
from rubin_sim.data import get_baseline

from tools import ObsTime, time_coord, galactic_extinction
from tools_rubin_sim import compute_mags, df_obs, real_obs, GRBObsTime




def generate_configs(N, popType='realistic', specType=0, b=4, p=2.2, epsilon_e=0.1, epsilon_B=0.01, xi_N=1.0, filename=None):

    """ Generate and save configurations

    :param N: number of configurations
    :param popType: type of the wanted GRB population ('realistic' for a realistic population or 'boosted' for an energy boosted population). Default value is 'realistic'
    :param specType: type of emission spectrum (0 for global cooling time + no inverse compton and 1 for global cooling time + inverse compton)
    :param filename: name of the Pickle file containing all the configurations ('configs_thetaC_N.pkl' by default)
    :other parameters: default values of the parameters of the model

    :return: 1 Pickle file of configurations for each value of thetaCore (thetaCore = 0.05 and 0.15 so 2 Pickle files)
    """


    Z = {'jetType': grb.jet.PowerLaw,  # Jet Type
         'specType': specType,  # Emission Spectrum
         'b': b,  # Power Law index
         'thetaObs': 0.2,  # Viewing angle in radians
         'E0': 1.0e51,  # Isotropic-equivalent energy in erg
         'thetaWing': 0.15,  # Truncation angle in radians
         'thetaCore': 0.05,  # Half-opening angle in radians
         'n0': 0.1,  # Circumburst density in cm^{-3}
         'p': p,  # Electron energy distribution index
         'epsilon_e': epsilon_e,  # epsilon_e
         'epsilon_B': epsilon_B,  # epsilon_B
         'xi_N': xi_N,  # Fraction of electrons accelerated
         'd_L': 1.0e28,  # Luminosity distance in cm
         'z': 0.5}  # Redshift


    if filename == None:
        file1 = open('../data/simulations/configs_005_' + str(N) + '.pkl', 'wb')
        file2 = open('../data/simulations/configs_015_' + str(N) + '.pkl', 'wb')
    else:
        file1 = open(f'{filename}_005_{N}.pkl', 'wb')
        file2 = open(f'{filename}_015_{N}.pkl', 'wb')


    z_range = np.linspace(0.001, 0.1, 100)

    def z(x):
        return (1 + 3.1 * x) / (1 + (x / 2.5) ** 3)

    def normalisation(x):
        return simps(z(x), x)

    class z_distrib_gen(rv_continuous):
        def _pdf(self, x, const):
            return (1.0 / const) * z(x)

    z_distrib = z_distrib_gen(name="redshift_distribution", a=0.001)

    norm_constant = normalisation(z_range)

    z_samples = z_distrib.rvs(const=norm_constant, size=N)

    all_dico_005 = []
    all_dico_015 = []

    for i in tqdm(range(N)):

        #Z['thetaObs'] = np.random.uniform(0, np.pi/2)
        Z['thetaObs'] = np.arccos(np.random.uniform(0, 1))

        if popType == 'realistic':
            E = np.random.normal(51, 1)
        elif popType == 'boosted':
            E = np.random.uniform(53, 55)

        Z['E0'] = 1.0 * 10 ** E
        Z['n0'] = 1.0 * 10 ** (-np.random.uniform(-1, 2))
        Z['z'] = z_samples[i]
        Z['d_L'] = cosmo.luminosity_distance(Z['z']).value * 3.08e24

        Z['thetaCore'] = 0.05
        Z['thetaWing'] = math.radians(np.exp(np.random.normal(1.94, 0.5)))

        while Z['thetaWing'] < Z['thetaCore']:
            Z['thetaWing'] = math.radians(np.exp(np.random.normal(1.94, 0.5)))
        all_dico_005.append(Z.copy())

        Z['thetaCore'] = 0.15

        while Z['thetaWing'] < Z['thetaCore']:
            Z['thetaWing'] = math.radians(np.exp(np.random.normal(1.94, 0.5)))
        all_dico_015.append(Z.copy())

    pickle.dump(all_dico_005, file1)
    pickle.dump(all_dico_015, file2)

    file1.close()
    file2.close()


def calculate_results(N, t, thetaC, freq=5.0e14, jetType='PL', filename_in=None, filename_out=None):

    """ Calculate and save some results from the configurations

    :param N: number of configurations
    :param t: time for which you want to calculate the light curve
    :param thetaC: value of thetaCore you want to study the afterglow ('005' for 0.05 radians or '015' for 0.15 radians)
    :param nu: value of the frequency of the observer in Hertz. Default is 5.0e14 Hz corresponding to the r filter
    :param jetType: type of simulated jet ('PL' for Power-Law, 'G' for Gaussian and 'TH' for Top-Hat). Default value is 'PL'

    :return: 1 Pickle file of results from the configurations for the chosen value of thetaCore
    """

    if filename_out == None:
        file = open('../data/simulations/simulations_' + jetType + '_' + str(thetaC) + '_' + str(N) + '.pkl', 'wb')
    else:
        file = open(f'{filename_out}.pkl', 'wb')

    if filename_in == None:
        file_open = open('../data/simulations/configs_' + str(thetaC) + '_' + str(N) + '.pkl', 'rb')
    else:
        file_open = open(f'{filename_in}.pkl', 'rb')

    configs_open = pickle.load(file_open)

    Z_results = {'jetType': -1,  # Jet Type
                 'config': {},  # Dictionary Z
                 'time': [],  # Time
                 'lc': [],  # Flux
                 't_obs': 0,  # Observability duration in days
                 'mag_min': 0,  # Minimum magnitude of the afterglow
                 'axis': 'on',  # On-axis or Off-axis
                 'observable': 'no',  # Observability of the afterglow : no or yes for >4. or >7. in days
                 'F_gamma': '0'}  # Flux of photons above 30 MeV in ph/cm2/s

    all_Z_results = []
    dt = []

    for i in tqdm(range(len(configs_open))):

        nu = np.empty(t.shape)
        nu[:] = freq

        Z = configs_open[i]

        if jetType == 'PL':
            Z['jetType'] = grb.jet.PowerLaw
            Z_results['jetType'] = 'PowerLaw'
        elif jetType == 'G':
            Z['jetType'] = grb.jet.Gaussian
            Z_results['jetType'] = 'Gaussian'
        elif jetType == 'TH':
            Z['jetType'] = grb.jet.TopHat
            Z_results['jetType'] = 'TopHat'

        fnu = grb.fluxDensity(t, nu, **Z)
        mag = -2.5 * np.log10(fnu * 1.0e-26) - 48.6
        dt.append(ObsTime(t, mag))

        Z_results['config'] = Z
        Z_results['time'] = t * grb.sec2day
        Z_results['lc'] = fnu
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
        t_max = t[imax[0]]

        if len(t_max) == 1:

            # Range of LAT
            nua = 1.0e8 * 1.6e-19 / 6.62e-34  # 10^5 keV
            nub = 5.0e10 * 1.6e-19 / 6.62e-34  # 5 x 10^ keV

            nu_lat = np.geomspace(nua, nub, num=100)
            Fnu_max = grb.fluxDensity(t_max, nu_lat, **Z)  # Flux en mJy

            delta_nu = nub - nua

            Nnu = []  # Flux calculated with wavelength
            for i in range(len(nu_lat)):
                Nnu.append(Fnu_max[i] * 10 ** -29 / (nu_lat[i] * 6.62e-34))

            Z_results['F_gamma'] = sum(Nnu) * delta_nu / 10 ** 4  # Flux in photon/cm2/s

        else:
            Z_results['F_gamma'] = 0

        all_Z_results.append(Z_results.copy())

    pickle.dump(all_Z_results, file)

    file_open.close()
    file.close()


def open_results(N, thetaC, jetType='PL', filename=None):

    """ Shows results from the configurations

    :param N: number of configurations
    :param thetaC: value of thetaCore you want to study the afterglow ('005' for 0.05 radians or '015' for 0.15 radians)
    :param jetType: type of simulated jet ('PL' for Power-Law, 'G' for Gaussian and 'TH' for Top-Hat). Default value is 'PL'

    :return: a Panda Dataframe containing the results from the configurations for the chosen value of thetaCore
    """

    if filename == None:
        file_open = open('../data/simulations/simulations_' + jetType + '_' + str(thetaC) + '_' + str(N) + '.pkl', 'rb')
    else:
        file_open = open(f'{filename}.pkl', 'rb')

    configs_open = pickle.load(file_open)
    file_open.close()

    return pd.DataFrame(configs_open)


def generate_pseudo_obs(N, thetaC, path_data, path_dustmaps, jetType='PL', filename_in=None, filename_out=None, extinction=True):

    """ Generate and save pseudo-observations of the lights curves

    :param N: number of configurations
    :param jetType: type of simulated jet ('PL' for Power-Law, 'G' for Gaussian and 'TH' for Top-Hat). Default value is 'PL'
    :param extinction: whether you want to take into account the galactic extinction or not. Takes the values 'True' or 'False', default is 'True'

    :return: 1 Pickle file containing a dictionary for each config with the configuration, the GRB coordinates and observation date, the observed magnitude, the   limiting magnitude at the observation time and the filter used at the moment of the detection
    """

    # or replace with any local file you want to use
    baseline_db = get_baseline()
    print(baseline_db)

    conn = sqlite3.connect(baseline_db)

    # In the near future, 'summaryallprops' will be replaced with 'observations'
    df = pd.read_sql('select * from observations;', conn)

    conn.close()

    # path_data = '/home/masson/rubin_sim_data'
    fdir = os.path.join(path_data, 'throughputs', 'baseline')

    # Read the throughput curves
    filterlist = ['u', 'g', 'r', 'i', 'z', 'y']
    filtercolors = {'u': 'b', 'g': 'c', 'r': 'g', 'i': 'orange', 'z': 'r', 'y': 'm'}

    lsst = {}
    for f in filterlist:
        lsst[f] = Bandpass()
        lsst[f].readThroughput(os.path.join(fdir, f'total_{f}.dat'))

    if filename_in == None:
        file_open = open('../data/simulations/simulations_' + jetType + '_' + str(thetaC) + '_' + str(N) + '.pkl', 'rb')
    else:
        file_open = open(f'{filename_in}.pkl', 'rb')
    configs_open = pickle.load(file_open)
    file_open.close()

    configs = pd.DataFrame(configs_open)
    configs = configs[(configs['axis'] == 'off') & (configs['t_obs'] > 7.)]
    # configs = configs[configs['mag_min']<22]

    if filename_out == None:
        file = open('../data/pseudo_obs/pseudo_obs_' + jetType + '_' + str(thetaC) + '_' + str(N) + '.pkl', 'wb')
    else:
        file = open(f'{filename_out}.pkl', 'wb')


    # dictionary containing all the GRB information
    LC = {'config': {},  # Dictionary Z
          'grb_time': 0,  # GRB observation date
          'grb_coord': 0,  # GRB ra/dec coordinates
          'time': [],  # Time of each detection
          'mags': [],  # Magnitude of the detection
          'filt': [],  # Filter used at the moment of the detection
          'mags_lim': [],  # Limiting magnitude of the detection at the observation time
          'mags_err': []}   # Error on the magnitude

    all_LC = []

    for config in tqdm(configs['config']):

        Z = config
        # print(Z)

        # Time in mjd and ra/dec coordinates of the GRB
        grb_time, grb_coord = time_coord()
        # print(grb_time.isot, grb_coord.to_string('hmsdms'))

        # grb_time = 61402.36838867448
        # grb_time = Time(grb_datetime, format='mjd', scale='utc')

        # grb_ra = "03h38m30.03464941s"
        # grb_dec = "-54d36m49.59379716s"
        # grb_coord = SkyCoord(grb_ra, grb_dec, frame='icrs')

        # Observe t_before days before and t_after days after
        t_before = TimeDelta(20, format='jd')
        t_after = TimeDelta(365, format='jd')
        obs_start = grb_time - t_before
        obs_end = grb_time + t_after

        # Get time span
        df_time = df[(df['observationStartMJD'] > obs_start.mjd) & (df['observationStartMJD'] < obs_end.mjd)]

        # Angular separation with SkyCoord.separation
        # Rubin FOV is 47 square degree for a 3.5-degree diameter, hence 1.7 deg separation radius.
        df_time['Separation'] = SkyCoord(df_time['fieldRA'], df_time['fieldDec'], unit="deg").separation(
            grb_coord).degree
        df_sky = df_time[df_time['Separation'] < 1.7]

        # print(df_sky)

        time_bins = GRBObsTime(df_sky, grb_time)

        # print(time_bins.to_dict())

        obs_list = df_obs(Z, df_sky, time_bins, lsst)

        y_mags = []    # magnitude with extinction

        a_lambda_u, a_lambda_g, a_lambda_r, a_lambda_i, a_lambda_z, a_lambda_y = galactic_extinction(grb_coord,
                                                                                                    path_dustmaps)

        # If there is no observation, let's go to the next configuration
        if len(obs_list) == 0:
            all_LC.append(0)
            continue

        else:
            obs_df = pd.concat(obs_list)

            obs_df['observationId'] = df_sky['observationId']

            x_times, y_mags_without_ext, z_colors, mags_lim, mags_err = real_obs(obs_df, df_sky, time_bins, grb_time, lsst)

            if extinction == True:

                for i in range(len(y_mags_without_ext)):
                    if z_colors[i] == 'b':
                        y_mags.append(y_mags_without_ext[i] + a_lambda_u)
                    elif z_colors[i] == 'c':
                        y_mags.append(y_mags_without_ext[i] + a_lambda_g)
                    elif z_colors[i] == 'g':
                        y_mags.append(y_mags_without_ext[i] + a_lambda_r)
                    elif z_colors[i] == 'orange':
                        y_mags.append(y_mags_without_ext[i] + a_lambda_i)
                    elif z_colors[i] == 'r':
                        y_mags.append(y_mags_without_ext[i] + a_lambda_z)
                    else:
                        y_mags.append(y_mags_without_ext[i] + a_lambda_y)

                LC['mags'] = y_mags

            elif extinction == False:

                LC['mags'] = y_mags_without_ext

            LC['config'] = Z
            LC['grb_time'] = grb_time.isot
            LC['grb_coord'] = grb_coord.to_string('hmsdms')
            LC['time'] = x_times
            LC['filt'] = z_colors
            LC['mags_lim'] = mags_lim
            LC['mags_err'] = list(np.array(mags_err)[:, 0, 0])

            all_LC.append(LC.copy())

    pickle.dump(all_LC, file)
    file.close()

