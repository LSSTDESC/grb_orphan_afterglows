"""Tools module

This module provides a number of generic tools.

"""


import numpy as np


def flux_to_mag(flux):
    ''' Convert flux from Jansky to AB Magnitude

    1 Jy = 1e-23 erg/cm2/s/Hz
    Fnu = 3631 Jy = 3.631*1e-20 erg/cm2/s/Hz
    ABmag = 0-2.5*log10( Fnu )-48.6 = 0

    :param flux: flux in milli-Jansky
    :return: mag: as the AB Magnitude
    '''
    mag = -2.5 * np.log10(flux*1.0e-23) - 48.6
    return mag