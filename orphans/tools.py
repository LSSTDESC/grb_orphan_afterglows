"""Tools module

This module provides a number of generic tools.

"""


import numpy as np
from scipy.constants import lambda2nu


def flux_to_mag(flux):
    """ Convert flux from Jansky to AB Magnitude

    1 Jy = 1e-23 erg/cm2/s/Hz
    Fnu = 3631 Jy = 3.631*1e-20 erg/cm2/s/Hz
    ABmag = 0-2.5*log10( Fnu )-48.6 = 0

    :param flux: flux in milli-Jansky
    :return: mag: as the AB Magnitude
    """
    mag = -2.5 * np.log10(flux*1.0e-23) - 48.6
    return mag


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
