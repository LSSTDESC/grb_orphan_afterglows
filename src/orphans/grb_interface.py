''' GRB Interface module

This module provides the main interface to afterglowpy
'''

# Standard library imports
from copy import deepcopy

# Third‑party imports
import numpy as np
import afterglowpy as grb
from astropy.cosmology import Planck18 as cosmo

# Local imports
from orphans.tools import get_wl_and_nu_band
from orphans.grb_configs import GRB_BASE_PARAMS


def make_grb_light_curve(E0=1.0e53, thetaObs=0.05, thetaCore=0.1, freq=5.0e14):
    """ Compute GRB light curve in a single frequency band.

    This function calculates the flux density as a function of time for a fixed frequency.

    :param E0: float
        The isotropic equivalent energy of the burst (in erg).
    :param thetaObs: float
        The observer viewing angle (in radians).
    :param thetaCore: float
        The jet opening angle (in radians).
    :param freq: float
        The fixed observing frequency for this calculation (in Hz).
    :return: tuple of (nu, t, Fnu_Jy) where nu is the frequency array, t is the time array, and Fnu_Jy is the flux in Jansky.
    """
    # for convenience, place arguments into a dict.
    Z = deepcopy(GRB_BASE_PARAMS)
    Z['E0'] = E0
    Z['thetaObs'] = thetaObs
    Z['thetaCore'] = thetaCore

    # space time points geometrically, from 10^3 s to 10^7 s
    t = np.geomspace(1.0e3, 1.0e7, 300)

    # calculate flux in a single band (all times have same frequency)
    nu = np.empty(t.shape)
    nu[:] = freq

    # calculate but Fnu is in mJy by default
    fnu = grb.fluxDensity(t, nu, **Z)
    # so we convert to Jy
    Fnu_Jy = fnu * 1.0e-3
    return nu, t, Fnu_Jy


def make_grb_spectrum(
    jetType=4,
    E0=1.0e53,
    z=1,
    n0=1.,
    thetaObs=0.05,
    thetaCore=0.1,
    thetaWing=0.15,
    specType=0,
    t=1.0 * grb.day2sec,
):  # pylint: disable=too-many-arguments,too-many-positional-arguments
    """ Compute GRB Spectral Energy Distribution (SED).

    This function computes the flux density across a full wavelength band based on
    GRB physical parameters and cosmological model.

    :param jetType: int
        The jet opening angle or structure type determining the light curve shape.
    :param E0: float
        Total isotropic equivalent energy of the burst (in erg).
    :param z: float
        Redshift of the GRB.
    :param n0: float
        Ambient medium density in cm^-3.
    :param thetaObs: float
        Observer viewing angle (in radians).
    :param thetaCore: float
        Core opening angle (in radians).
    :param thetaWing: float
        Wing opening angle (in radians).
    :param specType: int
        The spectral model type to use for flux calculation. (e.g., 0 for a specific case)
    :param t: float or array-like of float
        Time points for which the spectrum is calculated, specified as one day in seconds.
    :return: tuple of (wl_full_band, freq_full_band, t, Fnu_Jy) where wl/freq are wavelength/frequency arrays in nm/Hz, and Fnu_Jy is the flux density in Jansky.
    """
    # For convenience, place arguments into a dict.
    Z = deepcopy(GRB_BASE_PARAMS)
    Z['jetType'] = jetType
    Z['specType'] = specType
    Z['E0'] = E0
    Z['z'] = z
    Z['d_L'] = cosmo.luminosity_distance(Z['z']).value * 3.08e24
    Z['n0'] = n0
    Z['thetaObs'] = thetaObs
    Z['thetaCore'] = thetaCore
    Z['thetaWing'] = thetaWing
    # first create a wavelength range from 200 to 1300 nm
    wl_full_band, freq_full_band = get_wl_and_nu_band()
    # calculate but Fnu is in mJy by default
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
    with open(file_name, 'w', encoding='utf-8') as f:
        f.write("# lambda(nm)   Fnu(Jy)\n")
        for wl, fnu in zip(wavelenghts, Fnu_Jy):
            f.write(f'{wl:.1f}\t{fnu:.6f}\n')
    return 0