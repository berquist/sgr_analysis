#!/usr/bin/env python3

from __future__ import print_function
from __future__ import division

import os

from glob import glob

import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def getargs():

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('filename', nargs='+')

    args = parser.parse_args()

    return args

if __name__ == '__main__':

    args = getargs()

    # filenames = sorted(glob('PESmats/*.dat'))
    # sorted(glob('RigidPESmats/*.dat'))

    # for filename in filenames:
    #     print(filename)

    # filename = 'PESmats/PES_940.dat'

    for filename in args.filename:
        stub = os.path.splitext(filename)[0]

        mat_file_contents = np.loadtxt(filename)
        # print(mat_file_contents.shape)

        mat_l1 = mat_file_contents[:, 0]
        mat_l2 = mat_file_contents[:, 1]
        mat_en = mat_file_contents[:, 2]

        print(mat_l1)
        print(mat_l2)
        print(mat_en)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_trisurf(mat_l1, mat_l2, mat_en, cmap=plt.get_cmap('coolwarm'))
        ax.set_xlabel(r'$l_1$')
        ax.set_ylabel(r'$l_2$')
        ax.set_zlabel(r'energy (cm$^{-1}$)')
        fig.savefig('{}.pdf'.format(stub), bbox_inches='tight')

        plt.close('all')

        # fig, ax = plt.subplots()

        # xx, yy = np.meshgrid(mat_l1, mat_l2)
        # ax.pcolor(xx, yy, mat_en, cmap=plt.get_cmap('coolwarm'))
        # ax.set_xlabel(r'$l_1$')
        # ax.set_ylabel(r'$l_2$')
        # ax.set_zlabel(r'energy (cm$^{-1}$)')
        # fig.savefig('{}_pcolor.pdf'.format(stub), bbox_inches='tight')

        # plt.close('all')
