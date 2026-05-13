"""Jetsimpy parameters configuration

This mirrors the configuration used in the original masson/jetsimpy implementation.
"""

GRB_BASE_PARAMS = {
    'Eiso': 1e53,               # isotropic equivalent energy (erg)
    'lf': 300,                  # Lorentz factor
    'theta_c': 0.1,             # jet half-opening angle (rad)
    'n0': 1,                    # ISM number density (cm^{-3})
    'A': 0,                     # wind density amplitude
    'eps_e': 0.1,               # epsilon_e
    'eps_b': 0.01,              # epsilon_b
    'p': 2.2,                   # electron power-law index
    'theta_v': 0.2,             # viewing angle (rad)
    's': 4,                     # power-law jet slope
    'd': 3273.436515882976,    # distance (Mpc) placeholder
    'z': 0.55,                  # redshift
}
