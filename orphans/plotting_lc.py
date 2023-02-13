"""Plotting afterglow light curve module

"""

import numpy as np
import afterglowpy as grb
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt

from tools import ObsTime, mag_to_flux, GalacticExtinction


def PlotSimulation(t, fnu, flux='mag'):
    """ Plot theoretical light curve of an orphan afterglow

    :param t: time for which the flux is calculated
    :param fnu: value of the flux of the afterglow in mJy
    :param flux: whether you want to plot the flux in magnitude or in mJy. Takes the values 'mag' or 'flux', default is 'mag'

    :return: plot of the afterglow light curve
    """

    plt.figure(figsize=(15, 12))
    plt.rcParams.update({'font.size': 22})

    if flux == 'mag':
        mag = -2.5 * np.log10(fnu * 1.0e-26) - 48.6
        plt.plot(t, mag, color='royalblue')
        plt.ylabel('AB Magnitude')
        plt.gca().invert_yaxis()
        plt.ylim(30, 15)
        plt.axhline(y=24.5, color='black', linewidth=0.5, label='Rubin/LSST nightly limiting magnitude')

    elif flux == 'flux':
        plt.plot(t, fnu, color='royalblue')
        plt.ylabel('Fux (mJy)')
        plt.yscale('log')
        plt.ylim(mag_to_flux(30) * 1.0e3, mag_to_flux(15) * 1.0e3)
        plt.axhline(y=mag_to_flux(24.5) * 1.0e3, color='black', linewidth=0.5,
                    label='Rubin/LSST nightly limiting magnitude')

    plt.xlabel('Time (days)')
    plt.xscale('log')

    plt.show()


def PlotPseudoObs(lc, path_dustmaps, flux='mag', extinction=True):
    """ Plot pseudo-observed light curve of an afterglow

    :param lc: dictionary containing the pseudo-observed light curve, calculated with the [...] function (in the [...] module)
    :param flux: whether you want to plot the flux in magnitude or in mJy. Takes the values 'mag' or 'flux', default is 'mag'
    :param extinction: whether you want to take into account the galactic extinction or not. Takes the values 'True' or 'False', default is 'True'

    :return: plot of the afterglow pseudo-observed light curve
    """

    plt.rcParams["figure.figsize"] = [15, 12]

    filterlist = ['u', 'g', 'r', 'i', 'z', 'y']
    filtercolors = {'u': 'b', 'g': 'c', 'r': 'g', 'i': 'orange', 'z': 'r', 'y': 'm'}

    x_times = lc['time']
    z_colors = lc['filt']

    mags = lc['mags']

    true_mags = []

    coord = lc['grb_coord']

    a_lambda_u, a_lambda_g, a_lambda_r, a_lambda_i, a_lambda_z, a_lambda_y = GalacticExtinction(coord, path_dustmaps)

    for i in range(len(mags)):
        if z_colors[i] == 'b':
            true_mags.append(mags[i] + a_lambda_u)
        elif z_colors[i] == 'c':
            true_mags.append(mags[i] + a_lambda_g)
        elif z_colors[i] == 'g':
            true_mags.append(mags[i] + a_lambda_r)
        elif z_colors[i] == 'orange':
            true_mags.append(mags[i] + a_lambda_i)
        elif z_colors[i] == 'r':
            true_mags.append(mags[i] + a_lambda_z)
        else:
            true_mags.append(mags[i] + a_lambda_y)

    if flux == 'mag':
        if extinction == True:
            y_mags = true_mags
        else:
            y_mags = lc['mags']
        mags_lim = lc['mags_lim']
        plt.gca().invert_yaxis()
        plt.ylabel('Observed magnitude', fontsize=28)
    elif flux == 'flux':
        if extinction == True:
            y_mags = list(map(mag_to_flux, true_mags))
        else:
            y_mags = list(map(mag_to_flux, lc['mags']))
        mags_lim = list(map(mag_to_flux, lc['mags_lim']))
        plt.ylabel('Observed flux (mJy)', fontsize=28)
        plt.yscale('log')

    # print(x_times, y_mags, mags_lim)

    # plot pseudo observed light curve
    for x, y, z, m in zip(x_times, y_mags, z_colors, mags_lim):
        if flux == 'mag' and y < m:
            plt.scatter(x, y, c=z, s=100)
            plt.scatter(x, m, c=z, marker='2', s=100)
        elif flux == 'flux' and y > m:
            plt.scatter(x, y, c=z, s=100)
            plt.scatter(x, m, c=z, marker='2', s=100)

    plt.title('Pseudo observed light curve', fontsize=28)
    plt.xlabel('Time in MJD', fontsize=28)

    legend_elements_2 = list()

    for filt in filterlist:
        fcolor = filtercolors[filt]
        legend_elements_2.append(Line2D([0], [0], marker='o', color=fcolor, label=filt,
                                        markerfacecolor=fcolor, markersize=10))
    plt.legend(handles=legend_elements_2)
    plt.show()
