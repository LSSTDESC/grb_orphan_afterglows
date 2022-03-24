import numpy as np
import matplotlib.pyplot as plt
import afterglowpy as grb
import math 
from astropy.cosmology import Planck18 as cosmo


Z = {'jetType':     grb.jet.PowerLaw,     # Top-Hat jet
     'specType':    0,                  # Basic Synchrotron Emission Spectrum
     'thetaObs':    0.2,   # Viewing angle in radians
     'b':           4,
     'E0':          1.0e53, # Isotropic-equivalent energy in erg
     'thetaWing':   0.15,    # Truncature angle
     'thetaCore':   0.1,    # Half-opening angle in radians
     'n0':          1.0,    # circumburst density in cm^{-3}
     'p':           2.2,    # electron energy distribution index
     'epsilon_e':   0.1,    # epsilon_e
     'epsilon_B':   0.01,   # epsilon_B
     'xi_N':        1.0,    # Fraction of electrons accelerated
     'd_L':         1.0e28, # Luminosity distance in cm
     'z':           0.55}   # redshift


t = np.geomspace(1.0e2, 1.0e8, 300)
nu = 5.0e14


thetaObs = [0.16, 0.2, 0.25, 0.3, 0.4]
colors = ['mediumpurple', 'mediumaquamarine', 'orange', 'tomato', 'sienna']

redshift = [0.5, 1., 1.5, 2., 2.5, 3.]
linestyle = ['-',  (0, (5, 1)), '-.', (0, (5, 10)), ':', (0, (1, 10))]


print("Half-opening angle is %.2f°." %(math.degrees(Z['thetaCore'])))


fig, ax = plt.subplots(1, 1, figsize = (8,5))

for theta in thetaObs:
    Z['thetaObs'] = theta
    for z in redshift:
        Z['z'] = z
        Z['d_L'] = cosmo.luminosity_distance(z).value*3.08e24
        Fnu = grb.fluxDensity(t, nu, **Z)
        mag = -2.5 * np.log10(Fnu*1.0e-26) - 48.6
        plt.plot(t*grb.sec2day, mag, label = r'$\theta_{obs}$ = %.2f°' %math.degrees(theta), color=colors[thetaObs.index(theta)], linestyle=linestyle[redshift.index(z)])
        if thetaObs.index(theta)==0:
            print("z = %.1f => DL = %i Mpc." %(z,cosmo.luminosity_distance(z).value))

lines = ax.get_lines()
legend1 = plt.legend([lines[i] for i in range(len(redshift))], ['z = %.1f' %redshift[i] for i in range(len(redshift))], loc=2)
legend2 = plt.legend([lines[i] for i in [j*len(redshift) for j in range(len(thetaObs))]], [r'$\theta_{obs}$ = %.2f°' %math.degrees(thetaObs[i]) for i in range(len(thetaObs))], loc=1)
ax.add_artist(legend1)    

plt.xscale('log')
plt.gca().invert_yaxis()
plt.ylim(30,15)
plt.axhline(y=24.5, color='black', linestyle=':', label='LSST limiting magnitude')
ax.set_xlabel('t (days)')
ax.set_ylabel('Magnitude')

fig.tight_layout()