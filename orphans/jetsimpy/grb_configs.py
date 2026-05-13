""" GRB parameters configuration for jetsimpy

    J. Bregeon - March 2022
"""


GRB_BASE_PARAMS = {'Eiso':  1e53,                # isotropic equivalent energy
    'lf': 300,                                   # lorentz factor
    'theta_c': 0.1,                              # jet half opening angle [rad]
    'n0': 1,                                     # ism number density
    'A': 0,                                      # wind number density amplitude
    'eps_e': 0.1,                                # epsilon_e
    'eps_b': 0.01,                               # epsilon_b
    'p': 2.2,                                    # electron power index
    'theta_v': 0.2,                              # viewing angle [rad]
    's': 4,                                      # power-law jet slope (required for power-law jet)
    'd': 3273.436515882976,                      # distance [Mpc]
    'z': 0.55,                                   # redshift
}