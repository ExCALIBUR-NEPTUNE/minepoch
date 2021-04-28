#!/usr/bin/env python3
import sys
import argparse
import numpy as np
from minepoch_py.energy import energy_check
from minepoch_py.two_stream import two_stream_analysis


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="""Basic CI testing script for minepoch
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

    group.add_argument(
        "--growth_rate",
        default=0.05,
        type=float,
        help="Allowed fractional error growth rate",
        nargs="?",
    )

    return parser.parse_args()


def run_tests(plot=True, savefig=True):

    options = parse_arguments()

    energy_error = energy_check(tolerance=options.energy_conservation,
                                plot=plot, savefig=savefig)

    rate, _ = two_stream_analysis(check_growth=True, plot=plot,
                                  savefig=savefig)

    # Analytic rate
    rate_analytic = 0.7
    rate_error = np.abs(rate - rate_analytic) / rate_analytic

    # Return 1 or 0, depending on error > tolerance
    if (energy_error > 0 or rate_error > options.growth_rate):
        if energy_error > 0:
            print('Energy conservation check failed!')
        if rate_error > options.growth_rate:
            print('Two stream growth check failed!')
            print('Observed growth rate = %.4E (expected %.2F)'
                  % (rate, rate_analytic))
        return 1

    return 0


if __name__ == "__main__":
    # Run error analysis
    sys.exit(run_tests())
