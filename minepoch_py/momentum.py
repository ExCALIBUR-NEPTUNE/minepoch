#!/usr/bin/env python3
import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np


def parse_arguments():
    """Parse command line options.

    """

    parser = argparse.ArgumentParser(
        description="""Momentum conservation analysis for minEPOCH
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    group = parser.add_argument_group("Tolerances")

    group.add_argument(
        "--momentum_conservation",
        default=1e-4,
        type=float,
        help="Allowed fractional error in momentum conservation",
        nargs="?",
    )

    parser.add_argument("--disable_plots", help="Disable plot generation",
                        action="store_true")

    return parser.parse_args()


def get_momenta(fname='Data/momentum.dat'):
    """Calculate time history of field and particle momentum for an output file.

    Args:
        fname: Input filename. Default: Data/momentum.dat

    Returns:
        times, particle_momentum, field_momentum
    """

    times = []
    field_momenta = []
    particle_momenta = []
    with open(fname) as f:
        # Skip header
        _ = f.readline()
        for line in f:
            times.append(line.split()[1])
            particle_momenta.append(line.split()[2:5])
            field_momenta.append(line.split()[5:])

    times = np.array(times, dtype=float)
    particle_momenta = np.array(particle_momenta, dtype=float)
    field_momenta = np.array(field_momenta, dtype=float)

    return times, particle_momenta, field_momenta


def momentum_check(fname='Data/momentum.dat', tolerance=None, label=None,
                 plot=True, savefig=False, signedplot=True):
    """Check quality of momentum conservation from minEPOCH simulation

    Args:
        fname: Input filename. Default: Data/momentum.dat
        tolerance: Maximum acceptable fractional error. Default: None (no check)
        label: Label for data set when plotting. Default: None
        plot: Control optional plotting. Default: True
        savefig: Control saving of figure. Default: False
	signedplot: Plot +/- errors, or just absolute values. Default: True
    """

    # Calculate total energy as a function of time
    times, pp, fp = get_momenta(fname)
    total_p = np.sum((pp + fp)**2, axis=1) # Magnitude as a function of time
    

    # If plotting requested produce plot of energy conservation error as
    # a function of time
    if plot:
        if signedplot:
            data = total_p - total_p[0]
        else:
            data = np.abs(total_p - total_p[0])
        plt.plot(times, data / total_p[0], label=label)
        plt.xlabel(r'$\mathrm{t (s)}$', fontsize=16)
        plt.ylabel(r'$\frac{\Delta |P|}{|P|(t=0)}$', fontsize=16)
        plt.xlim(left=times[0])

        # Update legend if necessary
        if label is not None:
            plt.legend()

        plt.tight_layout()

        if savefig:
            plt.savefig('momentum_conservation.png')

    # If tolerance provided, check final energy conservation error
    if tolerance is not None:
        delta_p = np.abs(total_p[-1] - total_p[0]) / total_p[0]
        if delta_p > tolerance:
            print('Momentum conservation (fractional) error = %.4E' % delta_p)
            return 1
        return 0

    return None


if __name__ == "__main__":
    # Parse arguments
    options = parse_arguments()
    # Run error analysis
    sys.exit(momentum_check(tolerance=options.momentum_conservation,
                            plot=not options.disable_plots,
                            savefig=True))
