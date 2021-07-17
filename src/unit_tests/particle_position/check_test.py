#!/usr/bin/env python3
from scipy import stats
import numpy as np
import sys
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="""Tool for validating uniform particle loading
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    group = parser.add_argument_group("Particle Distribution")

    group.add_argument(
	"--x0",
        default=0.0,
	type=float,
	help="Centre of region in x-coordinate",
        nargs="?",
    )

    group.add_argument(
	"--dx",
	default=1.0,
	type=float,
	help="Width of region in x",
	nargs="?",
    )

    group.add_argument(
	"--y0",
        default=2.0,
	type=float,
	help="Centre of region in y-coordinate",
        nargs="?",
    )

    group.add_argument(
	"--dy",
	default=2.0,
	type=float,
	help="Width of region in y",
	nargs="?",
    )

    group.add_argument(
	"--z0",
        default=4.0,
	type=float,
	help="Centre of region in z-coordinate",
        nargs="?",
    )

    group.add_argument(
	"--dz",
	default=3.0,
	type=float,
	help="Width of region in z",
	nargs="?",
    )

    parser.add_argument("--filename",
                        help="filename to check",
                        default='particle_pos.dat',
                        type=str,
                        nargs="?")

    parser.add_argument("--pvalue",
                        help="Minimum allowable pvalue",
                        default=0.05,
                        type=float,
                        nargs="?")

    return parser.parse_args()

def check_dist(data, x0, dx):

    dis = stats.uniform(loc = x0 - 0.5*dx, scale=dx)
    return stats.cramervonmises(data, dis.cdf)

def check_data(options):

    pos = np.loadtxt(options.filename)

    p = check_dist(pos[:,0], options.x0, options.dx).pvalue
    if (p < options.pvalue):
        print('x distribution check failed! pvalue = %.4E' % p)
        return 1

    p = check_dist(pos[:,1], options.y0, options.dy).pvalue
    if (p < options.pvalue):
        print('y distribution check failed! pvalue = %.4E' % p)
        return 1

    p = check_dist(pos[:,2], options.z0, options.dz).pvalue
    if (p < options.pvalue):
        print('z distribution check failed! pvalue = %.4E' % p)
        return 1

    print('Particle position checks passed!')

    return 0

if __name__ == "__main__":
    options = parse_arguments()
    # Run error analysis
    sys.exit(check_data(options))
