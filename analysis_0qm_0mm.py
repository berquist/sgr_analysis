#!/usr/bin/env python

from __future__ import print_function
from __future__ import division

import os
import pickle

from collections import OrderedDict

import numpy as np
import numpy.linalg as npl

from analysis_utils import get_CO2_frequencies
from analysis_utils import get_dipoles
from analysis_utils import get_outputfiles_from_path

import matplotlib as mpl
# mpl.rc(usetex=True)
mpl.use("Agg")
import matplotlib.pyplot as plt


def sort(snapnums_dict, results_dict):
    assert snapnums_dict.keys() == results_dict.keys()
    for k in keys:
        sorting_indices = [i[0] for i in sorted(enumerate(snapnums_dict[k]),
                                                key=lambda x: x[1])]
        sorted_results = [i[1] for i in sorted(zip(sorting_indices, results_dict[k]),
                                               key=lambda x: x[0])]
        sorted_snapnums = [i[1] for i in sorted(zip(sorting_indices, snapnums_dict[k]),
                                                key=lambda x: x[0])]
        # assert sorted_snapnums == list(range(min(snapnums_dict[k]), max(snapnums_dict[k]) + 1))
        snapnums_dict[k] = sorted_snapnums
        results_dict[k] = sorted_results
    return


def get_single_snapshot_results(snapnum, snapnums_dict, results_dict):
    assert snapnums_dict.keys() == results_dict.keys()
    single_dict = []
    for k in keys:
        idx = snapnums_dict[k].index(snapnum)
        single_result = results_dict[k][idx]
        single_dict.append((k, single_result))
    return OrderedDict(single_dict)


def plot_single_snapshot_dipoles(snapnum, snapnums_d, dipoles, inp_fig=None, inp_ax=None):

    dipoles_snap = get_single_snapshot_results(snapnum, snapnums_d, dipoles)

    fig, ax = plt.subplots()
    if inp_fig:
        fig = inp_fig
    if inp_ax:
        ax = inp_ax

    plot_list = [npl.norm(dipole) for dipole in dipoles_snap.values()]
    print(plot_list)
    ax.plot(ticks,
            plot_list,
            label=snapnum,
            marker='o')

    if not inp_ax:
        ax.tick_params(direction='out', top='off', right='off')
        ax.set_xticklabels(xticklabels)
        ax.set_xlabel('method', fontsize='large')
        ax.set_ylabel("total dipole moment (Debye)", fontsize='large')
        ax.set_title("snapshot {}".format(snapnum), fontsize='large')
        # ax.legend(loc='best', fancybox=True, framealpha=0.50)
    if not inp_fig:
        print('Saving dipole_snap{}.pdf'.format(snapnum))
        fig.savefig('dipole_snap{}.pdf'.format(snapnum), bbox_inches='tight')

        plt.close(fig)

    return


def plot_single_snapshot_frequencies(snapnum, snapnums_f, frequencies, inp_fig=None, inp_ax=None):

    frequencies_snap = get_single_snapshot_results(snapnum, snapnums_f, frequencies)

    fig, ax = plt.subplots()
    if inp_fig:
        fig = inp_fig
    if inp_ax:
        ax = inp_ax

    plot_list = [frequencies_snap[k] for k in keys]
    # plot_list = [frequency for frequency in frequencies_snap.values()]
    print(plot_list)
    ax.plot(ticks,
            plot_list,
            label=snapnum,
            marker='o')

    if not inp_ax:
        ax.tick_params(direction='out', top='off', right='off')
        ax.set_xticklabels(xticklabels)
        ax.set_xlabel('method', fontsize='large')
        ax.set_ylabel(r'$\nu_{3}$ frequency (cm$^{-1}$)', fontsize='large')
        ax.set_title("snapshot {}".format(snapnum), fontsize='large')
        # ax.legend(loc='best', fancybox=True, framealpha=0.50)
    if not inp_fig:
        print('Saving frequency_snap{}.pdf'.format(snapnum))
        fig.savefig('frequency_snap{}.pdf'.format(snapnum), bbox_inches='tight')

        plt.close(fig)

    return


def filter_outputfiles(l):
    return list(filter(lambda x: '_0mm' in x, l))


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("file_operation", choices=("none", "save", "read"))
    parser.add_argument("parse_operation", choices=("none", "save", "read"))
    parser.add_argument('--debug', action='store_true')

    args = parser.parse_args()

    keys = (
        'hf',
        'b3lyp',
        'wb97x-d',
        'blyp',
        'tpss',
        # 'ri-mp2',
    )

    labels = OrderedDict([
        ('blyp', 'BLYP'),
        ('tpss', 'TPSS'),
        ('b3lyp', 'B3LYP'),
        ('wb97x-d', r'$\omega$B97X-D'),
        ('hf', 'HF'),
        ('ri-mp2', 'RI-MP2'),
    ])

    basis_sets = OrderedDict([
        ('6-31gdp', '6-31G**'),
        ('cc-pvtz', 'cc-pVTZ'),
    ])

    if args.file_operation == "save":
        print("Trying to find output files...")
        basedir = "/home/eric/Chemistry/calc.sgr/droplets/snapshot_method_dependence"
        outputfiles = dict()
        for basis_set in basis_sets:
            outputfiles[basis_set] = []
            for k in keys:
                curr_outputfiles = filter_outputfiles(get_outputfiles_from_path(os.path.join(basedir, "inputs_freq_0qm_{k}_{basis_set}".format(k=k, basis_set=basis_set))))
                outputfiles[basis_set].append((k, curr_outputfiles))
            outputfiles[basis_set] = OrderedDict(outputfiles[basis_set])
        with open('outputfiles.pypickle', 'wb') as picklefile:
            pickle.dump(outputfiles, picklefile)
    elif args.file_operation == "read":
        print("Reading list of output files from: {}".format(os.path.abspath("outputfiles.pypickle")))
        with open("outputfiles.pypickle", "rb") as picklefile:
            outputfiles = pickle.load(picklefile)
    elif args.file_operation == "none":
        pass
    else:
        raise Exception

    if args.debug:
        print(outputfiles)

    if args.parse_operation == "save":
        print("Extracting valuable information from outputs...")
        frequencies = dict()
        intensities = dict()
        snapnums_f = dict()
        dipoles = dict()
        snapnums_d = dict()
        print("Parsing frequencies/intensities...")
        for basis_set in basis_sets:
            frequencies[basis_set] = []
            intensities[basis_set] = []
            snapnums_f[basis_set] = []
            dipoles[basis_set] = []
            snapnums_d[basis_set] = []
            for k in outputfiles[basis_set]:
                k_frequencies, k_intensities, k_snapnums_f = get_CO2_frequencies(outputfiles[basis_set][k])
                frequencies[basis_set].append((k, k_frequencies))
                intensities[basis_set].append((k, k_intensities))
                snapnums_f[basis_set].append((k, k_snapnums_f))
        print("Parsing dipoles...")
        for basis_set in basis_sets:
            dipoles[basis_set] = []
            snapnums_d[basis_set] = []
            for k in outputfiles[basis_set]:
                k_dipoles, k_snapnums_d = get_dipoles(outputfiles[basis_set][k])
                dipoles[basis_set].append((k, k_dipoles))
                snapnums_d[basis_set].append((k, k_snapnums_d))
        for basis_set in basis_sets:
            frequencies[basis_set] = OrderedDict(frequencies[basis_set])
            intensities[basis_set] = OrderedDict(intensities[basis_set])
            snapnums_f[basis_set] = OrderedDict(snapnums_f[basis_set])
            dipoles[basis_set] = OrderedDict(dipoles[basis_set])
            snapnums_d[basis_set] = OrderedDict(snapnums_d[basis_set])
        with open('frequencies.pypickle', 'wb') as picklefile:
            pickle.dump(frequencies, picklefile)
        with open('intensities.pypickle', 'wb') as picklefile:
            pickle.dump(intensities, picklefile)
        with open('dipoles.pypickle', 'wb') as picklefile:
            pickle.dump(dipoles, picklefile)
        with open('snapnums_frequencies.pypickle', 'wb') as picklefile:
            pickle.dump(snapnums_f, picklefile)
        with open('snapnums_dipoles.pypickle', 'wb') as picklefile:
            pickle.dump(snapnums_d, picklefile)
    elif args.parse_operation == "read":
        print("Reading frequency data from: {}".format(os.path.abspath('frequencies.pypickle')))
        with open('frequencies.pypickle', 'rb') as picklefile:
            frequencies = pickle.load(picklefile)
        print("Reading intensity data from: {}".format(os.path.abspath('intensities.pypickle')))
        with open('intensities.pypickle', 'rb') as picklefile:
            intensities = pickle.load(picklefile)
        print("Reading dipole data from: {}".format(os.path.abspath('dipoles.pypickle')))
        with open('dipoles.pypickle', 'rb') as picklefile:
            dipoles = pickle.load(picklefile)
        print("Reading snapshot number data from: {}".format(os.path.abspath('snapnums_frequencies.pypickle')))
        with open('snapnums_frequencies.pypickle', 'rb') as picklefile:
            snapnums_f = pickle.load(picklefile)
        print("Reading snapshot number data from: {}".format(os.path.abspath('snapnums_dipoles.pypickle')))
        with open('snapnums_dipoles.pypickle', 'rb') as picklefile:
            snapnums_d = pickle.load(picklefile)
    elif args.parse_operation == "none":
        pass
    else:
        raise Exception

    ticks = range(len(keys))
    xticklabels = list(labels[k] for k in keys)

    cutoff = -1
    means_frequencies = [np.mean(frequencies[k][:cutoff])
                         for k in keys]
    means_dipole_moments = [np.mean([npl.norm(dipole) for dipole in dipoles[k][:cutoff]])
                            for k in keys]

    ##################################################

    fig, ax = plt.subplots()

    ax.plot(ticks, means_frequencies, marker='s', color='green')

    ax.set_xticks(ticks)
    ax.set_xticklabels(xticklabels)
    ax.set_xlabel('method', fontsize='large')
    ax.set_ylabel(r'$\nu_{3}$ frequency (cm$^{-1}$)', fontsize='large')

    fig.savefig('frequencies_line.pdf', bbox_inches='tight')

    ##################################################

    fig, ax = plt.subplots()

    ax.plot(ticks, means_dipole_moments, marker='s', color='green')

    ax.set_xticks(ticks)
    ax.set_xticklabels(xticklabels)
    ax.set_xlabel('method', fontsize='large')
    ax.set_ylabel('total dipole moment (Debye)', fontsize='large')

    fig.savefig('dipole_moments_line.pdf', bbox_inches='tight')

    ##################################################

    ax1_color = 'blue'
    ax2_color = 'green'

    fig, ax1 = plt.subplots()

    ax1.plot(ticks, means_frequencies, marker='s', linestyle='-', color=ax1_color)

    ax1.set_ylabel(r'$\nu_{3}$ frequency (cm$^{-1}$)', fontsize='large')

    for tl in ax1.get_yticklabels():
        tl.set_color(ax1_color)

    ax2 = ax1.twinx()

    ax2.plot(ticks, means_dipole_moments, marker='*', linestyle='--', color=ax2_color)

    ax2.set_ylabel('total dipole moment (Debye)', fontsize='large')

    ax2.set_ylim(ax2.get_ylim()[::-1])

    for tl in ax2.get_yticklabels():
        tl.set_color(ax2_color)

    ax1.set_xticks(ticks)
    ax1.set_xticklabels(xticklabels)
    ax1.set_xlabel('method', fontsize='large')

    ax1.tick_params(direction='out', top='off', right='off')
    ax2.tick_params(direction='out', top='off', left='off')

    fig.savefig('frequencies_dipole_moments_combined_line.pdf', bbox_inches='tight')

    ##################################################

    snapnums = (25, 50, 75, 100)
    # snapnums = [random.randrange(1, cutoff + 1) for _ in range(5)]

    # for snapnum in snapnums:
    #     print(snapnum, get_single_snapshot_results(snapnum, snapnums_f, frequencies))

    # for snapnum in snapnums:
    #     plot_single_snapshot_frequencies(snapnum, snapnums_f, frequencies)
    #     plot_single_snapshot_dipoles(snapnum, snapnums_d, dipoles)

    ##################################################

    fig, ax = plt.subplots()
    for snapnum in snapnums:
        plot_single_snapshot_dipoles(snapnum, snapnums_d, dipoles, inp_ax=ax, inp_fig=fig)
    ax.tick_params(direction='out', top='off', right='off')
    ax.set_xticks(ticks)
    ax.set_xticklabels(xticklabels)
    ax.set_xlabel('method', fontsize='large')
    ax.set_ylabel("total dipole moment (Debye)", fontsize='large')
    ax.legend(loc='best', fancybox=True, framealpha=0.50)
    print('Saving dipole_snapshots.pdf')
    fig.savefig('dipole_snapshots.pdf', bbox_inches='tight')

    plt.close(fig)

    fig, ax = plt.subplots()
    for snapnum in snapnums:
        plot_single_snapshot_frequencies(snapnum, snapnums_d, frequencies, inp_ax=ax, inp_fig=fig)
    ax.tick_params(direction='out', top='off', right='off')
    ax.set_xticks(ticks)
    ax.set_xticklabels(xticklabels)
    ax.set_xlabel('method', fontsize='large')
    ax.set_ylabel(r'$\nu_{3}$ frequency (cm$^{-1}$)', fontsize='large')
    ax.legend(loc='best', fancybox=True, framealpha=0.50)
    print('Saving frequency_snapshots.pdf')
    fig.savefig('frequency_snapshots.pdf', bbox_inches='tight')

    plt.close(fig)

    ###

    cutoff = 100
    # xticklabels_traj = list(snapnums_f.values())[:cutoff]
    # ticks_traj = range(1, len() + 1)

    fig, ax = plt.subplots()

    for k in keys:
        ax.plot(frequencies[k][:cutoff], label=labels[k])

    # ax.set_xticks(ticks_traj)
    # ax.set_xticklabels(xticklabels_traj)
    ax.set_xlabel('snapshot #', fontsize='large')
    ax.set_ylabel(r'$\nu_{3}$ frequency (cm$^{-1}$)', fontsize='large')
    ax.tick_params(direction='out', top='off', right='off')
    ax.legend(fancybox=True, loc='upper right', framealpha=0.50)

    fig.savefig('frequencies.pdf', bbox_inches='tight')

    ##################################################
