#!/usr/bin/env python3
import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np


def parse_arguments():
    """Parse command line options.

    """

    parser = argparse.ArgumentParser(
        description="""Energy conservation analysis for minEPOCH
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    group = parser.add_argument_group("Tolerances")

    group.add_argument(
        "--energy_conservation",
        default=1e-4,
        type=float,
        help="Allowed fractional error in energy conservation",
        nargs="?",
    )

    parser.add_argument("--disable_plots", help="Disable plot generation",
                        action="store_true")

    return parser.parse_args()


def get_energies(fname='Data/output.dat'):
    """Calculate time history of field and particle energies for an output file.

    Args:
        fname: Input filename. Default: Data/output.dat

    Returns:
        times, particle_energies, field_energies
    """

    times = []
    field_energies = []
    particle_energies = []
    with open(fname) as f:
        # Skip header
        _ = f.readline()
        for line in f:
            times.append(line.split()[1])
            particle_energies.append(line.split()[5])
            field_energies.append(line.split()[6])

    times = np.array(times, dtype=float)
    particle_energies = np.array(particle_energies, dtype=float)
    field_energies = np.array(field_energies, dtype=float)

    return times, particle_energies, field_energies


def plot_field_energy(fname='Data/output.dat', xscale=1.0, dbg=False, ylog=False,
                      xlog=False, **kwargs):
    """Plot the field energy as a function of time.

    Args:
        fname: Input filename. Default: Data/output.dat
        xscale: Scaling for x-axis. Default: 1 (no scaling)
        dbg: Produce debugging information. Default: False
        ylog: Plot logarithmic y-axis. Default: False
        xlog: Plot logarithmic x-axis. Default: False
    """

    # Read energies from file. Ignore particle energy
    times, _, energies = get_energies(fname)

    # If debug requested print maximum energy
    if dbg:
        print("Maximum energy = ", np.max(energies))

    plt.plot(times * xscale, energies, **kwargs)

    if ylog:
        plt.yscale('log')

    plt.xlim(0, times[-1] * xscale)
    if xlog:
        plt.xscale('log')

    plt.xlabel('Time')
    plt.ylabel('Field Energy')
    plt.tight_layout()


def energy_check(fname='Data/output.dat', tolerance=None, label=None,
                 plot=True, savefig=False):
    """Check quality of energy conservation from minEPOCH simulation

    Args:
        fname: Input filename. Default: Data/output.dat
        tolerance: Maximum acceptable fractional error. Default: None (no check)
        label: Label for data set when plotting. Default: None
        plot: Control optional plotting. Default: True
        savefig: Control saving of figure. Default: False
    """

    # Calculate total energy as a function of time
    times, pe, fe = get_energies(fname)
    total_energy = pe + fe

    # If plotting requested produce plot of energy conservation error as
    # a function of time
    if plot:
        plt.plot(times, (total_energy - total_energy[0]) / total_energy[0],
                 label=label)
        plt.xlabel(r'$\mathrm{t (s)}$', fontsize=16)
        plt.ylabel(r'$\frac{\Delta \epsilon}{\epsilon(t=0)}$', fontsize=16)
        plt.xlim(left=times[0])

        # Update legend if necessary
        if label is not None:
            plt.legend()

        plt.tight_layout()

        if savefig:
            plt.savefig('energy_conservation.png')

    # If tolerance provided, check final energy conservation error
    if tolerance is not None:
        delta_e = np.abs(total_energy[-1] - total_energy[0]) / total_energy[0]
        if delta_e > tolerance:
            print('Energy conservation (fractional) error = %.4E' % delta_e)
            return 1
        return 0

    return None


if __name__ == "__main__":
    # Parse arguments
    options = parse_arguments()
    # Run error analysis
    sys.exit(energy_check(tolerance=options.energy_conservation,
                          plot=not options.disable_plots,
                          savefig=True))
