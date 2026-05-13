"""Tools module

This module provides a number of generic tools.

"""


import numpy as np
from scipy.constants import lambda2nu
import afterglowpy as grb
import pandas as pd

from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
from astropy import units as u
from dustmaps.sfd import SFDQuery



def flux_to_mag(flux):
    """ Convert flux from Jansky to AB Magnitude

    1 Jy = 1e-23 erg/cm2/s/Hz
    Fnu = 3631 Jy = 3.631*1e-20 erg/cm2/s/Hz
    ABmag = 0-2.5*log10( Fnu )-48.6 = 0

    :param flux: flux in milli-Jansky
    :return: mag: as the AB Magnitude
    """

    mag = -2.5 * np.log10(flux*1.0e-26) - 48.6
    return mag




def mag_to_flux(mag):
    """ Convert flux from AB Magnitude to milli-Jansky

    1 Jy = 1e-23 erg/cm2/s/Hz
    Fnu = 3631 Jy = 3.631*1e-20 erg/cm2/s/Hz
    ABmag = 0-2.5*log10( Fnu )-48.6 = 0

    :param mag: as the AB Magnitude
    :return: flux: flux in milli-Jansky
    """
    
    flux = pow(10, (26 - (mag + 48.6) / 2.5))
    return flux




def get_wl_and_nu_band(wl_min=200, wl_max=1300):
    """ Get arrays for a given wavelength band

    Parameters
    ----------
    wl_min : `int`
        lower end of the wavelength band in nm
    wl_max : `int`
        higher end of the wavelength band in nm

    Returns
    -------
    wl_full_band : `array` of `int`
        a `numpy.array` of wavelengths
    freq_full_band : `array` of `float`
        a `numpy.array` of frequencies
    """
    # first create a wavelength range from 200 to 1300 nm
    wl_full_band = np.arange(wl_min, wl_max)
    wl_full_band_nm = wl_full_band * 1e-9
    # then convert to frequency
    freq_full_band = np.array(lambda2nu(wl_full_band_nm))
    # verify
    #print(f'wl_600_nm = {wl_full_band_nm[400]} nm = nu_5.0e14_Hz = {freq_full_band[400]:e} Hz')
    return wl_full_band, freq_full_band




def ObsTime(t, mag):
    """ Get observability duration for a orphan GRB light curve
    
    Parameters
    ----------
    t : `list`
        duration for which the light curve is calculated
    mag : `list`
        magnitude of the light curve at each time

    Returns
    -------
    dt : `float` 
        observability duration
    """
    
    index = []
    
    for i in range(len(mag)-1):
        if mag[i] < 24.5:
            index.append(i)
    if not index:
        return 0
    else:
        dt = (t[max(index)] - t[min(index)])*grb.sec2day
        return dt
        
        
        
        
def time_coord():
    """ Take randomly the time and ra/dec of a GRB during LSST

    Returns
    -------
    grb_time, grb_coord : GRB detection time during LSST and ra/dec coordinates
    """
    
    grb_datetime = np.random.uniform(59945, 63598) #time in mjd : t_start = 01/01/2023, t_end = 01/01/2033
    grb_time = Time(grb_datetime, format='mjd', scale='utc')
    
    grb_ra = np.random.uniform(-180.00,180.00)*u.degree 
    grb_dec = np.random.uniform(-90.00,0.00)*u.degree
    grb_coord = SkyCoord(grb_ra, grb_dec, frame='icrs') #ra/dec coordinates in Southern Hemisphere
    
    return grb_time, grb_coord




def galactic_extinction(coord, path_dustmaps):

    sfd = SFDQuery()
    c = SkyCoord(coord)
    ebv = sfd(c)

    #df = pd.read_csv('/home/masson/Documents/galactic_extinction/schlafly_dust_factor.csv', delimiter=',')
    df = pd.read_csv(path_dustmaps, delimiter=',', skiprows=1)
    df_lsst = df[(df.Bandpass == 'LSST_u')
                 | (df.Bandpass == 'LSST_g')
                 | (df.Bandpass == 'LSST_r')
                 | (df.Bandpass == 'LSST_i')
                 | (df.Bandpass == 'LSST_z')
                 | (df.Bandpass == 'LSST_y')]

    a_lambda_u = ebv * df_lsst.R_V31[38]
    a_lambda_g = ebv * df_lsst.R_V31[39]
    a_lambda_r = ebv * df_lsst.R_V31[40]
    a_lambda_i = ebv * df_lsst.R_V31[41]
    a_lambda_z = ebv * df_lsst.R_V31[42]
    a_lambda_y = ebv * df_lsst.R_V31[43]

    return a_lambda_u, a_lambda_g, a_lambda_r, a_lambda_i, a_lambda_z, a_lambda_y




def pseudo_obs_with_points(pseudo_obs, n_pts=1):
    """ Get the pseudo-observed light curves that have at least n_pts observable point(s)

    Parameters
    ----------
    pseudo_obs: `list` of `dict`
        pseudo-observed light curve for each configuration
    n_pts: `int`
        number of points on the pseudo-observed light curve. Default is 1

    Returns
    -------
    lc_with_pts: `list` of `dict`
        `list` of the pseudo-observed light curves that have at least n_pts observable point(s)
    """

    lc_obs = [pseudo_obs[i] for i in range(len(pseudo_obs)) if pseudo_obs[i] != 0]

    lc_with_pts = []

    for i in range(len(lc_obs)):

        mags = np.array(lc_obs[i]['mags'])
        mags_lim = np.array(lc_obs[i]['mags_lim'])

        sub = np.subtract(mags, mags_lim)

        if np.sum(np.array(sub) < 0) >= n_pts:
            lc_with_pts.append(lc_obs[i])

    return lc_with_pts
