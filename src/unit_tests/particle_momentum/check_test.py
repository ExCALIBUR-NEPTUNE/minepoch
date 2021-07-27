#!/usr/bin/env python3
from scipy import stats
from scipy.constants import elementary_charge as q0, electron_mass as m0, c
import numpy as np
import sys
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="""Tool for validating Maxwellian particle initialisation.
        Temperatures specified in eV, drift velocities in units of c.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    group = parser.add_argument_group("Particle Distribution")

    group.add_argument(
        "--temp_x",
        default=10.0,
        type=float,
        help="x-temperature (eV)",
        nargs="?",
    )

    group.add_argument(
        "--drift_x",
        default=0.01,
        type=float,
        help="x-dirift velocity (/c)",
        nargs="?",
    )

    group.add_argument(
        "--temp_y",
        default=1.0,
        type=float,
        help="y-temperature (eV)",
        nargs="?",
    )

    group.add_argument(
        "--drift_y",
        default=0.1,
        type=float,
        help="y-drift velocity (/c)",
        nargs="?",
    )

    group.add_argument(
        "--temp_z",
        default=200.0,
        type=float,
        help="z-temperature (eV)",
        nargs="?",
    )

    group.add_argument(
        "--drift_z",
        default=0.2,
        type=float,
        help="z-drift velocity (/c)",
        nargs="?",
    )

    parser.add_argument("--filename",
                        help="filename to check",
                        default='particle_mom.dat',
                        type=str,
                        nargs="?")

    parser.add_argument("--pvalue",
                        help="Minimum allowable pvalue",
                        default=0.05,
                        type=float,
                        nargs="?")

    return parser.parse_args()


def check_dist(data, drift, temp):

    dis = stats.norm(loc=drift, scale=temp)
    return stats.cramervonmises(data, dis.cdf)


def check_data(options):

    mom = np.loadtxt(options.filename)

    p = check_dist(mom[:, 0], options.drift_x * c * m0,
                  np.sqrt(options.temp_x * q0 * m0)).pvalue
    if (p < options.pvalue):
        print('x momentum check failed! pvalue = %.4E' % p)
        return 1

    p = check_dist(mom[:, 1], options.drift_y * c * m0,
                   np.sqrt(options.temp_y * q0 * m0)).pvalue
    if (p < options.pvalue):
        print('y momentum check failed! pvalue = %.4E' % p)
        return 1

    p = check_dist(mom[:, 2], options.drift_z * c * m0,
                   np.sqrt(options.temp_z * q0 * m0)).pvalue
    if (p < options.pvalue):
        print('z momentum check failed! pvalue = %.4E' % p)
        return 1

    print('Particle momentum checks passed!')

    return 0


if __name__ == "__main__":
    options = parse_arguments()
    # Run error analysis
    sys.exit(check_data(options))
