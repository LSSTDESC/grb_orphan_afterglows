""" GRB Interface module

This module provides the main interface to afterglowpy
"""

import math
import numpy as np
from copy import deepcopy
import jetsimpy
from astropy.cosmology import Planck18 as cosmo

#sys.path.append('/pbs/home/m/mmasson/lsst/orphans/orphans')

from modules.tools import flux_to_mag, get_wl_and_nu_band
from modules.grb_configs import GRB_BASE_PARAMS


def make_grb_light_curve(E0=1.0e53, thetaObs=0.05, thetaCore=0.1, lf=300, freq=5.0e14):
    """ Compute GRB light curve

    Note that the Flux is in mJy

    :param thetaObs: Observer angle
    :param thetaCore: Jet opening angle
    :param freq: Light frequency
    :return: arrays of frequency, time and fluxes in Jy
    """

    # for convenience, place arguments into a dict.
    P = deepcopy(GRB_BASE_PARAMS)
    P['Eiso'] = E0
    P['theta_v'] = thetaObs
    P['theta_c'] = thetaCore
    P['lf'] = lf

    # space time points geometrically, from 10^3 s to 10^7 s
    tsecond = np.geomspace(1.0e3, 1.0e7, 100)

    # calculate flux in a single band (all times have same frequency)
    nu = np.empty(tsecond.shape)
    nu[:] = freq

    # jet without spreading
    jet = jetsimpy.Jet(
        jetsimpy.PowerLaw(P["theta_c"], P["Eiso"], lf0=P["lf"]),  # jet profile
        P["A"],  # wind number density scale
        P["n0"],  # ism number density scale
        spread=False,  # w/wo spreading effect
        grid=jetsimpy.ForwardJetRes(P["theta_c"], 129)  # resolution
    )

    # calculate the afterglow flux density (unit: mJy)
    fnu = jet.FluxDensity(
        tsecond,  # [second] observing time span
        nu,  # [Hz]     observing frequency
        P,  # parameter dictionary
    )

    # so we convert to Jy
    Fnu_Jy = fnu * 1.0e-3
    return nu, tsecond, Fnu_Jy


def make_grb_spectrum(E0=1.0e53, z=1, n0=1., thetaObs=0.05, thetaCore=0.1, lf=300, tday=1.0):
    """ Compute GRB SED
    """
    # For convenience, place arguments into a dict.
    P = deepcopy(GRB_BASE_PARAMS)
    P['Eiso'] = E0
    P['theta_v'] = thetaObs
    P['theta_c'] = thetaCore
    P['lf'] = lf
    P['n0'] = n0
    P['z'] = z
    P['d'] = cosmo.luminosity_distance(z).value

    # first create a wavelength range from 200 to 1300 nm
    wl_full_band, freq_full_band = get_wl_and_nu_band()

    # jet without spreading
    jet = jetsimpy.Jet(
        jetsimpy.PowerLaw(P["theta_c"], P["Eiso"], lf0=P["lf"]),  # jet profile
        P["A"],  # wind number density scale
        P["n0"],  # ism number density scale
        spread=False,  # w/wo spreading effect
        grid=jetsimpy.ForwardJetRes(P["theta_c"], 129)  # resolution
    )

    # calculate the afterglow flux density (unit: mJy)
    fnu = jet.FluxDensity(
        tday*3600*24,  # [second] observing time span
        freq_full_band,  # [Hz]     observing frequency
        P,  # parameter dictionary
    )

    # so we convert to Jy
    Fnu_Jy = fnu * 1.0e-3
    return wl_full_band, freq_full_band, tday, Fnu_Jy


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
