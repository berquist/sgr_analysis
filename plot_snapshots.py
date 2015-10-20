#!/usr/bin/env python3

from __future__ import print_function
from __future__ import division

import pickle

import numpy as np

from analysis_utils import mangle_dict_keys
from analysis_utils import read_snapshot_file


def filter_snapshots(desired_snapshot_numbers, snapnums_d, results_d):

    assert snapnums_d.keys() == results_d.keys()
    for n_qm in snapnums_d:
        assert snapnums_d[n_qm].keys() == results_d[n_qm].keys()
        for n_mm in snapnums_d[n_qm]:
            snapnums = snapnums_d[n_qm][n_mm]
            if len(snapnums) > 0:
                results = results_d[n_qm][n_mm]
                indices = [i for i, x in enumerate(snapnums)
                           if x in desired_snapshot_numbers]
                # assert sorted(np.asarray(snapnums)[indices]) == sorted(desired_snapshot_numbers)
                new_snapnums = np.asarray(snapnums)[indices]
                new_results = np.asarray(results)[indices]
                snapnums_d[n_qm][n_mm] = new_snapnums
                results_d[n_qm][n_mm] = new_results

    return


def sort(n_qm, n_mm, snapnums_d, results_d):

    z = sorted([(sn, r) for (sn, r) in zip(snapnums_d[n_qm][n_mm],
                                           results_d[n_qm][n_mm])
                 if sn in snapnums])

    z_sn, z_r = [p[0] for p in z], [p[1] for p in z]

    return z_sn, z_r


if __name__ == "__main__":

    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt

    with open('frequencies.pypickle', 'rb') as picklefile:
        frequencies_CO2_d = pickle.load(picklefile)
    with open('snapnums_frequencies.pypickle', 'rb') as picklefile:
        snapnums_frequencies_d = pickle.load(picklefile)

    frequencies_CO2_d = mangle_dict_keys(frequencies_CO2_d)
    snapnums_frequencies_d = mangle_dict_keys(snapnums_frequencies_d)

    assert frequencies_CO2_d.keys() == snapnums_frequencies_d.keys()
    for n_qm in frequencies_CO2_d:
        assert frequencies_CO2_d[n_qm].keys() == snapnums_frequencies_d[n_qm].keys()

    # For dictionary access:
    # n_qm, then n_mm

    snapnums = read_snapshot_file("/home/eric/Chemistry/calc.sgr/droplets/inputs_freq/representative_snapshots_3qm")

    z3_sn, z3_f = sort(3, 0, snapnums_frequencies_d, frequencies_CO2_d)
    z2_sn, z2_f = sort(2, 0, snapnums_frequencies_d, frequencies_CO2_d)
    z1_sn, z1_f = sort(1, 0, snapnums_frequencies_d, frequencies_CO2_d)
    z0_sn, z0_f = sort(0, 0, snapnums_frequencies_d, frequencies_CO2_d)

    print("np.mean(z0_f): {:7.2f} len_f: {:4d} len_sn: {:4d}".format(np.mean(z0_f), len(z0_f), len(z0_sn)))
    print("np.mean(z1_f): {:7.2f} len_f: {:4d} len_sn: {:4d}".format(np.mean(z1_f), len(z1_f), len(z1_sn)))
    print("np.mean(z2_f): {:7.2f} len_f: {:4d} len_sn: {:4d}".format(np.mean(z2_f), len(z2_f), len(z2_sn)))
    print("np.mean(z3_f): {:7.2f} len_f: {:4d} len_sn: {:4d}".format(np.mean(z3_f), len(z3_f), len(z3_sn)))

    fig, ax = plt.subplots()

    ax.plot(z0_sn, z0_f, marker='o', label='0 QM')
    ax.plot(z1_sn, z1_f, marker='o', label='1 QM')
    ax.plot(z2_sn, z2_f, marker='o', label='2 QM')
    ax.plot(z3_sn, z3_f, marker='o', label='3 QM')

    ax.set_xlabel('snapshot #')
    ax.set_ylabel(r'$\nu_3$ frequency (cm$^{-1}$)')

    ax.legend(fancybox=True, loc='best', framealpha=0.50)

    fig.savefig('plot_snapshots.pdf', bbox_inches='tight')

    x = range(4)

    fig, ax = plt.subplots()

    for i in range(len(z3_sn)):
        snapnum = z3_sn[i]
        frequencies = [z0_f[i], z1_f[i], z2_f[i], z3_f[i]]
        ax.plot(x, frequencies, marker='o', label=snapnum)

    fig.savefig('plot_snapshots2.pdf'.format(snapnum), bbox_inches='tight')

    # filter the results so only those appears that are from the
    # snapshot numbers we've read in
    # filter_snapshots(snapnums, snapnums_frequencies_d, frequencies_CO2_d)

    for n_qm in range(4):
        x = [n_mm
             for n_mm in snapnums_frequencies_d[n_qm]]
        y = [np.mean(frequencies_CO2_d[n_qm][n_mm])
             for n_mm in snapnums_frequencies_d[n_qm]]
        pairs = sorted(zip(x, y))
        n_mm, f = [p[0] for p in pairs], [p[1] for p in pairs]
        ax.plot(n_mm, f, marker='o', label=n_qm)

    ax.set_xscale('symlog', basex=2)
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
    ax.yaxis.set_major_formatter(y_formatter)
    ax.tick_params(direction='out', top='off', right='off')
    ax.set_xlabel("# IL pairs treated as point charges")
    ax.legend(loc='best', fancybox=True, framealpha=0.50)
    fig.savefig('frequency_convergence_limited.pdf', bbox_inches='tight')
