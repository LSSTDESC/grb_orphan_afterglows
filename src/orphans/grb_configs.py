""" Some GRB parameters configurations

    J. Bregeon - March 2022
"""
import afterglowpy as grb  # pylint: disable=import-error,no-member,c-extension-no-member,E1101


GRB_BASE_PARAMS = {
    'jetType': grb.jet.TopHat,  # Top-Hat jet  # pylint: disable=no-member,E1101,I1101,c-extension-no-member
    'b': 4,  # Power-law index
    'specType': 0,  # Basic Synchrotron Spectrum
    'thetaObs': 0.2,  # Viewing angle in radians
    'E0': 1.0e53,  # Isotropic-equivalent energy in erg
    'thetaCore': 0.1,  # Half-opening angle in radians
    'n0': 1.0,  # circumburst density in cm^{-3}
    'p': 2.2,  # electron energy distribution index
    'epsilon_e': 0.1,  # epsilon_e
    'epsilon_B': 0.01,  # epsilon_B
    'xi_N': 1.0,  # Fraction of electrons accelerated
    'd_L': 1.0e28,  # Luminosity distance in cm
    'z': 0.55}  # redshift
