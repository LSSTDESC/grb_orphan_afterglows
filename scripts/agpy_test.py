""" Testing afterglowpy

    J. Bregeon - March 2022
"""


import math
import numpy as np
import afterglowpy as grb
import matplotlib.pyplot as plt
from copy import deepcopy

from orphans.tools import flux_to_mag
from orphans.grb_configs import GRB_BASE_PARAMS


def make_grb(thetaObs=0.05, thetaCore=0.1, freq=5.0e14):
    """ Compute GRB light curve

    Note that the Flux is in mJy

    :param thetaObs: Observer angle
    :param thetaCore: Jet opening angle
    :param freq: Light frequency
    :return: arrays of frequency, time and fluxes in Jy
    """
    # For convenience, place arguments into a dict.
    Z = deepcopy(GRB_BASE_PARAMS)
    Z['thetaObs'] = thetaObs
    Z['thetaCore'] = thetaCore

    # Space time points geometrically, from 10^3 s to 10^7 s
    t = np.geomspace(1.0e3, 1.0e7, 300)

    # Calculate flux in a single band (all times have same frequency)
    nu = np.empty(t.shape)
    nu[:] = freq

    # Calculate but Fnu is in mJy by default
    Fnu = grb.fluxDensity(t, nu, **Z)
    # so we convert to Jy
    Fnu_Jy = Fnu*1.0e-3
    return nu, t, Fnu_Jy


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    grb_lcs = list()
    for theta_obs in [0.05, 0.1, 0.15, 0.2, 0.3]:
        one = make_grb(thetaObs=theta_obs)
        grb_lcs.append(one)
        # plt.plot(np.log10(one[1]), one[2], label='Theta_Obs = %.2f°' % (math.degrees(theta_obs)))
        # plt.plot(np.log10(one[1]), np.log10(one[2]), label='Theta_Obs = %.2f°'%(math.degrees(theta_obs)))
        plt.plot(np.log10(one[1]), flux_to_mag(one[2]), label='Theta_Obs = %.2f°' % (math.degrees(theta_obs)))
        plt.gca().invert_yaxis()
        plt.legend()
    plt.xlabel('Time (log10(s))')
    # plt.ylabel('Flux (log10(Fnu))')
    # plt.yscale('log')
    plt.ylabel('AB Magnitude')
    plt.ylim(30, 10)
    plt.axhline(y=24.5, color='black', linestyle=':', label='Rubin nightly magnitude')
    plt.title('Afterglow light curve for Theta_Half_Jet=5.7° at 600 nm')
    plt.show()

