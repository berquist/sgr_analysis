#!/usr/bin/env python3

from __future__ import print_function
from __future__ import division

import math
import numpy as np
import numpy.linalg as npl

# 1. read in the xyz file
# 2. the last three lines are the CO2, extract them ['C', 'O', 'O']

# 3. calculate O-C-O angle (theta)
# 4. calculate O-C + C-O distance (l12?)

# 5. calculate alpha(l12)
# 6. calculate beta(theta)

# 7. make [2, 2] matrix
# 8. diagonalize!


def distance(icoords, jcoords):
    """Calculate the distance between two 3-vectors (i,j)."""
    return npl.norm(icoords - jcoords)


def bond_angle(icoords, jcoords, kcoords):
    """Calculate the angle (in degrees) between three 3-vectors (i,j,k),
    where j is the central point.
    """
    e_ji = -(jcoords - icoords) / distance(jcoords, icoords)
    e_jk = -(jcoords - kcoords) / distance(jcoords, kcoords)
    dot = max(np.dot(e_ji, e_jk), -1.0)
    angle = math.acos(dot)
    return angle * (180/math.pi)


if __name__ == "__main__":
    from mbe.xyz_operations import read_xyz

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('xyzfilename', nargs='+')
    parser.add_argument('--alpha-slope', type=float, default=-3415.3)
    parser.add_argument('--alpha-intercept', type=float, default=9523.8)
    parser.add_argument('--beta-slope', type=float, default=1.12)
    parser.add_argument('--beta-intercept', type=float, default=-515.7)
    args = parser.parse_args()

    alpha_slope = args.alpha_slope
    alpha_intercept = args.alpha_intercept
    beta_slope = args.beta_slope
    beta_intercept = args.beta_intercept

    l_lower = []
    l_upper = []

    for xyzfilename in args.xyzfilename:
        _, _, atoms, coords = read_xyz(xyzfilename)

        coords_CO2 = np.array(coords[-3:])

        # Carbon is the first atom!
        l12 = distance(coords_CO2[0], coords_CO2[1]) + distance(coords_CO2[0], coords_CO2[2])
        angle = bond_angle(coords_CO2[1], coords_CO2[0], coords_CO2[2])
        theta = 180.0 - angle

        alpha = (alpha_slope * l12) + alpha_intercept
        beta = (beta_slope * theta) + beta_intercept

        hamiltonian = np.array(
            [[alpha, beta],
             [beta, alpha]]
        )

        eigvals = npl.eigvalsh(hamiltonian)

        print(xyzfilename, eigvals)

        l_lower.append(eigvals[0])
        l_upper.append(eigvals[1])

    lower_mean = np.mean(np.array(l_lower))
    lower_stdev = np.std(np.array(l_lower))
    upper_mean = np.mean(np.array(l_upper))
    upper_stdev = np.std(np.array(l_upper))

    print(' mean: {:10.2f} {:10.2f}'.format(lower_mean, upper_mean))
    print('stdev: {:10.2f} {:10.2f}'.format(lower_stdev, upper_stdev))
