#!/usr/bin/env python3

from __future__ import print_function
from __future__ import division

import sys
import os.path

import math
import numpy as np
import numpy.linalg as npl

from mbe.xyz_operations import read_xyz
from analysis_utils import get_outputfiles_from_path

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


def getargs():

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('xyzfilename', nargs='+')
    parser.add_argument('--do-plot', action='store_true')

    args = parser.parse_args()

    return args


def get_coords_CO2(xyzfilename):
    _, _, atoms, coords = read_xyz(xyzfilename)

    # CO2 is always the last 3 atoms in each file.
    coords_CO2 = np.array(coords[-3:])

    return coords_CO2


def calc_theta(coords_CO2):
    # Carbon is the first atom!
    angle = bond_angle(coords_CO2[1], coords_CO2[0], coords_CO2[2])
    theta = 180.0 - angle
    return theta


def calc_bond_lengths(coords_CO2):
    # Carbon is the first atom!
    d1 = distance(coords_CO2[0], coords_CO2[1])
    d2 = distance(coords_CO2[0], coords_CO2[2])
    return (d1, d2)


def calc_theta_batch(xyzfilenames):
    return [calc_theta(get_coords_CO2(xyzfilename))
            for xyzfilename in xyzfilenames]


def clyde_sort(l):
    return sorted(l, key=lambda x: int(os.path.basename(x)[7:-4]))


def print_list(l):
    for x in l:
        print(x)
    return


def write_list(l, filename):
    with open(filename, 'w') as fh:
        for x in l:
            fh.write('{}\n'.format(x))
    return


def dump(path, logfilename=sys.stdout, do_clyde_sort=False):
    print("outputfiles:", path)
    if do_clyde_sort:
        xyzfiles = clyde_sort(get_outputfiles_from_path(path, ext=".xyz"))
    else:
        xyzfiles = get_outputfiles_from_path(path, ext=".xyz")
    with open(logfilename, 'w') as fh:
        fh.write('{}\n'.format(path))
        for xyzfilename in xyzfiles:
            coords_CO2 = get_coords_CO2(xyzfilename)
            d1, d2 = calc_bond_lengths(coords_CO2)
            theta = calc_theta(coords_CO2)
            fh.write('{} {} {}\n'.format(theta, d1, d2))
    return



if __name__ == "__main__":

    paths = (
        ("/home/eric/Chemistry/calc.sgr/paper_02_CD_SC/NewXYZFiles", False),
        ("/home/eric/Dropbox/SharedwPitt/NewXYZFiles", True),
        ("/home/eric/Desktop/FullBox", True),
        ("/home/eric/Chemistry/calc.sgr/incorrect/droplets_in_box/FullBox", True),
        ("/home/eric/Chemistry/calc.sgr/incorrect/droplets_in_box/xyz", False),
        ("/home/eric/Chemistry/calc.sgr/incorrect/new/xyz", False),
    )

    for idx, path in enumerate(paths):
        dump(path[0], logfilename="path_{}".format(path[0].replace("/", "_")), do_clyde_sort=path[1])
