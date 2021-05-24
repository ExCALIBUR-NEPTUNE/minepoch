import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import electron_mass as me, elementary_charge as qe
from scipy.constants import epsilon_0 as ep0
from scipy.optimize import curve_fit
from minepoch_py.energy import plot_field_energy, get_energies


def two_stream_analysis(fname='Data/output.dat', plot=True, check_growth=True,
                        savefig=False):
    # Cold plasma frequency
    omega = np.sqrt(8e11 * qe * qe / me / ep0)

    if plot:
        plot_field_energy(fname, xscale=omega, ylog=True)
        # Plot analytic growth
        t0 = 0.0
        t1 = 30.0
        plt.plot([t0, t1], [1e-8, 1e-8*np.exp(0.7*t1)], linestyle='--',
                 color='black', label='Theory (0.70)')
        plt.ylim(top=1e-2)
        plt.gca().yaxis.set_ticks_position('both')
        plt.xlabel(r'$\mathrm{t}\omega_{pe}$', fontsize=16)
        plt.ylabel(r'$\mathrm{Field Energy}$', fontsize=16)
        plt.tight_layout()

    if check_growth:
        # Estimate linear growth rate
        times, _, fe = get_energies(fname)
        # Estimate end of linear growth
        end = np.where(fe > np.max(fe) * 0.99)[0][0]
        # Ignore first third
        start = end // 3
        if plot:
            # Show region of interest
            plt.gca().axvspan(times[start] * omega, times[end] * omega,
                              alpha=0.25, color='red')

        times_linear = times[start:end] * omega
        energy_linear = fe[start:end]

        # Perform fit in log space
        def curve(x, a, b):
            return a * x + b

        popt, pcov = curve_fit(curve, times_linear, np.log(energy_linear))
        if plot:
            # Plot linear fit, shifted slightly vertically
            fit = 2 * np.exp(popt[1]) * np.exp(times_linear * popt[0])
            plt.plot(times_linear, fit, color='green',
                     label='Observed (%.4F)' % popt[0])

    if plot:
        plt.legend(fontsize=14, loc='upper left')
        if savefig:
            plt.savefig('two_stream.png')

    if check_growth:
        return popt[0], np.sqrt(np.diag(pcov))[0]

    return None, None
