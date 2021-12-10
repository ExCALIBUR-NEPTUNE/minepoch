import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import electron_mass as me, elementary_charge as qe
from scipy.constants import epsilon_0 as ep0
from scipy.optimize import curve_fit
from minepoch_py.energy import plot_field_energy, get_energies


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="""Two stream instability analysis for minEPOCH.

        Compares the growth of the field energy to provided analytic rate.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--electron_number_density",
                        default=8e11,
                        type=float,
                        help="Electron number density (m^{-3})",
                        nargs="?")

    parser.add_argument("--analytic_rate",
                        default=0.7,
                        type=float,
                        help="Analytic growth rate",
                        nargs="?")

    parser.add_argument("--disable_plots", help="Disable plot generation",
                        action="store_true")

    return parser.parse_args()


def two_stream_analysis(fname='Data/output.dat', plot=True, check_growth=True,
                        savefig=False, ne=8e11, analytic_rate=0.7):
    """Two stream instability analysis

    Args:
        fname: Input filename. Default: Data/output.dat
        plot: Control plotting. Default: True
        check_growth: Control checking of growth rate. Default: True
        savefig: Control saving of any plots. Default: False
        ne: Electron number density. Default: 8e11 (m^{-3})
        analytic_rate: Analytic growth rate to compare to. Default: 0.7
    """

    # Cold plasma frequency
    omega = np.sqrt(ne * qe * qe / me / ep0)

    if plot:
        plot_field_energy(fname, xscale=omega, ylog=True)
        # Plot analytic growth between t0 and t1
        t0 = 0.0
        t1 = 30.0
        plt.plot([t0, t1], [1e-8, 1e-8*np.exp(analytic_rate*t1)],
                 linestyle='--', color='black',
                 label=('Theory (%.2F)' % analytic_rate))
        plt.ylim(top=1e-2)
        plt.gca().yaxis.set_ticks_position('both')
        plt.xlabel(r'$\mathrm{t}\omega_{pe}$', fontsize=16)
        plt.ylabel(r'$\mathrm{Field Energy}$', fontsize=16)
        plt.tight_layout()

    # If required, check growth rate against analytic value
    if check_growth:
        # Estimate linear growth rate
        times, _, fe = get_energies(fname)
        # Estimate end of linear growth
        end = np.where(fe > np.max(fe) * 0.99)[0][0]
        # Measure growth rate over 10 plasma frequencies
        start = np.argmin(np.abs(times - times[end] + 10 / omega))

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

    # If growth rate calculated return value and error.
    # Otherwise return None, None
    if check_growth:
        return popt[0], np.sqrt(np.diag(pcov))[0]

    return None, None


if __name__ == "__main__":
    # Parse arguments
    options = parse_arguments()
    # Carry out analysis
    rate, _ = two_stream_analysis(ne=options.electron_number_density,
                                  analytic_rate=options.analytic_rate,
                                  plot=not options.disable_plots,
                                  savefig=True)
