#!/usr/bin/env python3

from __future__ import print_function

import pickle

import numpy as np

from scipy.stats import linregress

from analysis_utils import mangle_dict_keys
from model_hamiltonian_frequencies import (distance, bond_angle)

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def getargs():

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--csv', action='store_true')

    args = parser.parse_args()

    return args


def fit_line(x, y):
    """Return slope, intercept, rsq of best fit line."""
    slope, intercept, r, p, stderr = linregress(x, y)
    return slope, intercept, r**2


def zip_frequencies_geometries(frequencies_l, snapnums_f_l, geometries_l, snapnums_ge_l):
    """"""

    intersection = set(snapnums_f_l).intersection(snapnums_ge_l)

    zip_f = list(filter(lambda x: x[0] in intersection,
                        zip(snapnums_f_l, frequencies_l)))
    zip_ge = list(filter(lambda x: x[0] in intersection,
                         zip(snapnums_ge_l, geometries_l)))

    zip_f_ge = list(zip([p[0] for p in zip_f],
                        [p[1] for p in zip_f],
                        [p[1] for p in zip_ge]))

    return zip_f_ge


if __name__ == '__main__':

    args = getargs()

    # Read in the pickle files that contain all the raw data.
    with open('frequencies.pypickle', 'rb') as picklefile:
        frequencies = mangle_dict_keys(pickle.load(picklefile))
    with open('snapnums_frequencies.pypickle', 'rb') as picklefile:
        snapnums_f = mangle_dict_keys(pickle.load(picklefile))
    with open('geometries.pypickle', 'rb') as picklefile:
        geometries = mangle_dict_keys(pickle.load(picklefile))
    with open('snapnums_geometries.pypickle', 'rb') as picklefile:
        snapnums_ge = mangle_dict_keys(pickle.load(picklefile))

    # print([(n_mm, len(frequencies[0][n_mm]))
    #        for n_mm in sorted(frequencies[0].keys())])
    # print([(n_mm, len(snapnums_f[0][n_mm]))
    #        for n_mm in sorted(snapnums_f[0].keys())])
    # print([(n_mm, len(geometries[0][n_mm]))
    #        for n_mm in sorted(geometries[0].keys())])
    # print([(n_mm, len(snapnums_ge[0][n_mm]))
    #        for n_mm in sorted(snapnums_ge[0].keys())])

    for n_qm in (0, ):
        assert frequencies[n_qm].keys() \
            == snapnums_f[n_qm].keys() \
            == geometries[n_qm].keys() \
            == snapnums_ge[n_qm].keys()

    n_qm = 0
    n_mm = 0
    zip_f_ge = zip_frequencies_geometries(frequencies[n_qm][n_mm],
                                          snapnums_f[n_qm][n_mm],
                                          geometries[n_qm][n_mm],
                                          snapnums_ge[n_qm][n_mm])

    list_frequencies = []
    list_l1 = []
    list_l2 = []
    list_o12 = []
    list_l12 = []
    list_dl = []
    list_theta = []

    results = dict()

    for (snapnum, frequency, geometry) in zip_f_ge:

        C, O1, O2 = 0, 1, 2

        d_C_O1 = distance(geometry[C], geometry[O1])
        d_C_O2 = distance(geometry[C], geometry[O2])
        d_O1_O2 = distance(geometry[O1], geometry[O2])

        bond_sum = d_C_O1 + d_C_O2
        bond_difference = abs(d_C_O1 - d_C_O2)

        angle = bond_angle(geometry[O1], geometry[C], geometry[O2])
        theta = 180.0 - angle

        l1 = d_C_O1
        l2 = d_C_O2
        o12 = d_O1_O2
        l12 = bond_sum
        dl = bond_difference

        #print('{:4d} {:4.2f}'.format(snapnum, frequency))

        list_frequencies.append(frequency)
        list_l1.append(l1)
        list_l2.append(l2)
        list_o12.append(o12)
        list_l12.append(l12)
        list_dl.append(dl)
        list_theta.append(theta)

    results['frequencies'] = np.array(list_frequencies)
    results['l1'] = np.array(list_l1)
    results['l2'] = np.array(list_l2)
    results['o12'] = np.array(list_o12)
    results['l12'] = np.array(list_l12)
    results['dl'] = np.array(list_dl)
    results['theta'] = np.array(list_theta)

    # frequency_harmonic = 2420.20
    # frequency_vci_10 = 2378.57
    # correction_coeff = frequency_harmonic - frequency_vci_10

    # results['frequencies_corrected'] = results['frequencies'] - (correction_coefficient)

    if args.csv:
        import csv
        csvfile = open('analysis_geometries.csv', 'w')
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow([
            'parameter',
            'slope',
            'intercept',
            'rsq',
        ])

    for k in (x for x in results if x != 'frequencies'):

        slope, intercept, rsq = fit_line(results[k], results['frequencies'])
        fit_text = '$f(x) = {:f}x + {:f}$\n$R^{{2}} = {:f}$'.format(slope, intercept, rsq)
        print('{:5s} {:10.4f} {:10.4f} {:.4f}'.format(k, slope, intercept, rsq))

        fig, ax = plt.subplots()

        ax.plot(results[k], results['frequencies'], marker='o', linestyle='', label=k)
        ax.text(0.50, 0.05, fit_text, transform=ax.transAxes)

        ax.set_ylabel(r'$\nu_{3}$ frequency (cm$^{-1}$)')

        ax.legend(loc='best', fancybox=True, framealpha=0.50)

        fig.savefig('correlation_{}.pdf'.format(k), bbox_inches='tight')

        if args.csv:
            csvwriter.writerow([
                k,
                slope,
                intercept,
                rsq,
            ])

    if args.csv:
        csvfile.close()
