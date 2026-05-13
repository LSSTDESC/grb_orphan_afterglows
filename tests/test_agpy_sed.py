"""Script to show how to generate and plot a GRB sed

    J. Bregeon - March 2022
"""
import matplotlib.pyplot as plt
from orphans.grb_interface import make_grb_spectrum, dump_wl_Fnu_spectrum


if __name__ == '__main__':
    # generate default SED
    wl_full_band, freq_full_band, t, Fnu_Jy = make_grb_spectrum()
    # dump on disk
    dump_wl_Fnu_spectrum(wl_full_band, Fnu_Jy, file_name='my_grb_sed.txt')
    # plot
    plt.rcParams["figure.figsize"] = [9, 9]
    fig, ax = plt.subplots(2, 1)
    axs = ax.ravel()
    # SED nu, Fnu
    _nu = axs[0].plot(freq_full_band, Fnu_Jy)
    axs[0].set_title(f"SED (nu, Fnu)")
    axs[0].set_xlabel('Frequency (Hz)')
    axs[0].set_ylabel('Flux (Jy)')
    # SED lambda, Fnu
    _lambda = axs[1].plot(wl_full_band, Fnu_Jy)
    axs[1].set_title(f"SED (lambda, Fnu)")
    axs[1].set_xlabel('Wavelength (nm)')
    axs[1].set_ylabel('Flux (Jy)')
    plt.tight_layout()
    plt.show()
