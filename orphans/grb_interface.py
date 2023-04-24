""" GRB Interface module

This module provides the main interface to afterglowpy
"""

import math
import numpy as np
from copy import deepcopy
import afterglowpy as grb


from orphans.tools import flux_to_mag, get_wl_and_nu_band
from orphans.grb_configs import GRB_BASE_PARAMS


def make_grb_light_curve(E0=1.0e53, thetaObs=0.05, thetaCore=0.1, freq=5.0e14):
    """ Compute GRB light curve

    Note that the Flux is in mJy

    :param thetaObs: Observer angle
    :param thetaCore: Jet opening angle
    :param freq: Light frequency
    :return: arrays of frequency, time and fluxes in Jy
    """
    # For convenience, place arguments into a dict.
    Z = deepcopy(GRB_BASE_PARAMS)
    Z['E0'] = E0
    Z['thetaObs'] = thetaObs
    Z['thetaCore'] = thetaCore

    # Space time points geometrically, from 10^3 s to 10^7 s
    t = np.geomspace(1.0e3, 1.0e7, 300)

    # Calculate flux in a single band (all times have same frequency)
    nu = np.empty(t.shape)
    nu[:] = freq

    # Calculate but Fnu is in mJy by default
    fnu = grb.fluxDensity(t, nu, **Z)
    # so we convert to Jy
    Fnu_Jy = fnu * 1.0e-3
    return nu, t, Fnu_Jy


def make_grb_spectrum(E0=1.0e53, z=1, n0=1., thetaObs=0.05, thetaCore=0.1, thetaWing=0.15, t=1.0 * grb.day2sec):
    """ Compute GRB SED
    1.0 * grb.day2sec is just 1 day
    """
    # For convenience, place arguments into a dict.
    Z = deepcopy(GRB_BASE_PARAMS)
    Z['E0'] = E0
    Z['z'] = z
    Z['n0'] = n0
    Z['thetaObs'] = thetaObs
    Z['thetaCore'] = thetaCore
    Z['thetaWing'] = thetaWing
    # first create a wavelength range from 200 to 1300 nm
    wl_full_band, freq_full_band = get_wl_and_nu_band()
    # Calculate but Fnu is in mJy by default
    fnu = grb.fluxDensity(t, freq_full_band, **Z)
    # so we convert to Jy
    Fnu_Jy = fnu * 1.0e-3
    return wl_full_band, freq_full_band, t, Fnu_Jy


def dump_wl_Fnu_spectrum(wavelenghts, Fnu_Jy, file_name="grb_sed.txt"):
    """ Get arrays for a given wavelength band

    Parameters
    ----------
    wavelenghts : `array` of `int`
        a `numpy.array` of wavelengths
    Fnu_Jy : `array` of `float`
        a `numpy.array` of fluxes in Jansky
    file_name : `string`
        the file path

    Returns
    -------
    0 : if file was properly written on disk
    """
    print(f"Writing {file_name}")
    with open(file_name, 'w') as f:
        f.write("# lambda(nm)   Fnu(Jy)\n")
        for wl, fnu in zip(wavelenghts, Fnu_Jy):
            f.write(f'{wl:.1f}\t{fnu:.6f}\n')
    return 0
