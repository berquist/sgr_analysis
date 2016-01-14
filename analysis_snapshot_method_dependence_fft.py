#!/usr/bin/env python3

from __future__ import print_function
from __future__ import division

import sys
sys.path.append('/home/eric/Chemistry/aimd')
sys.path.append('/home/eric/Chemistry/aimd/spectra')
from qchem_aimd_fft import \
    (make_vibspectrum_from_dipoles, AU_TIME_IN_SEC)

import pickle

import numpy as np

# from parse_outputs_snapshot_method_dependence_4fs import \
#     (methods, basis_sets)

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt


if __name__ == '__main__':

    with open('dipoles.pypickle', 'rb') as picklefile:
        dipoles = pickle.load(picklefile)
    with open('snapnums_dipoles.pypickle', 'rb') as picklefile:
        snapnums_d = pickle.load(picklefile)

    # time step is 50 ps
    # time_step = 50.0 * 1.0e-12 / AU_TIME_IN_SEC
    # time step is 4 fs
    time_step = 4.0 * 1.0e-15 / AU_TIME_IN_SEC
    print(time_step)

    # for method in methods:
    #     for basis_set in basis_sets:
    #         dipoles_method_basis_0_0 = np.array(dipoles[method][basis_set][0][0])
    #         frequencies, intensities = make_vibspectrum_from_dipoles(dipoles_method_basis_0_0, time_step)
    #         assert frequencies.shape == intensities.shape
    #         # print(frequencies.shape)
    #         # print(frequencies)

    #         fig, ax = plt.subplots()

    #         ax.plot(frequencies, intensities, linewidth=0.50)

    #         fig.savefig('fft_{}_{}.pdf'.format(method, basis_set), bbox_inches='tight')

    #         plt.close('all')

    dipoles_b3lyp_lp_0_0 = np.array(dipoles['flex']['b3lyp'][''][0][0])
    frequencies, intensities = make_vibspectrum_from_dipoles(dipoles_b3lyp_lp_0_0, time_step)
    assert frequencies.shape == intensities.shape
    print(frequencies.shape)

    fig, ax = plt.subplots()

    ax.plot(frequencies, intensities, linewidth=0.50)

    fig.savefig('fft_{}_{}.pdf'.format('b3lyp', 'lp'), bbox_inches='tight')

    plt.close('all')
