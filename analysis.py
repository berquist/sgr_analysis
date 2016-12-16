#!/usr/bin/env python3

"""analysis.py: Where most of the analysis for the 'droplet' snapshots
is.
"""

from __future__ import print_function
from __future__ import division

import pickle
import csv

from copy import deepcopy
from functools import partial

import numpy as np
import numpy.linalg as npl
import scipy.stats as sps

from analysis_utils import filter_snapshots
from analysis_utils import get_single_snapshot_results
from analysis_utils import mangle_dict_keys
from analysis_utils import pprint_linregress
from analysis_utils import read_snapshot_file
from analysis_utils import slice
from model_hamiltonian_frequencies import distance


def condon():
    """Testing whether or not the Condon approximation is appropriate."""

    fig, ax = plt.subplots()

    frequencies_all = []
    intensities_all = []

    csvfile = open('condon_analysis_linear_regression.csv', 'w')
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow([
        '# QM',
        '# MM',
        '# points',
        'slope',
        'intercept',
        'rsq',
    ])

    list_l12 = []
    geometries = geometries_d[0][0]
    C, O1, O2 = 0, 1, 2
    for geometry in geometries:
        d_C_O1 = distance(geometry[C], geometry[O1])
        d_C_O2 = distance(geometry[C], geometry[O2])
        d_O1_O2 = distance(geometry[O1], geometry[O2])
        bond_sum = d_C_O1 + d_C_O2
        # bond_difference = abs(d_C_O1 - d_C_O2)
        list_l12.append(bond_sum)
    list_l12 = np.array(list_l12)

    for n_qm in sorted(frequencies_CO2_d):
        print("Forming Condon approximation plot for {}".format(labels[n_qm]))
        frequencies_single_qm_all_mm = []
        intensities_single_qm_all_mm = []
        # This is only necessary to get this mess to work, so the list
        # lengths are correct. The CO2 geometry will always be the
        # same.
        geometries_single_qm_all_mm = []
        for n_mm in possible_keys:
            f = frequencies_CO2_d[n_qm][n_mm]
            i = intensities_CO2_d[n_qm][n_mm]
            s = snapnums_frequencies_d[n_qm][n_mm]
            # filter the geometry results based on the current
            # snapshots
            indices = [(snapnum - 1) for snapnum in s]
            g = list_l12[indices]
            assert len(f) == len(i) == len(g)
            frequencies_single_qm_all_mm.extend(f)
            intensities_single_qm_all_mm.extend(i)
            geometries_single_qm_all_mm.extend(g)
            frequencies_all.extend(f)
            intensities_all.extend(i)
            print('{} QM/{} MM'.format(n_qm, n_mm))
            try:
                slope, intercept, rsq = pprint_linregress(f, i)
                csvwriter.writerow([n_qm, n_mm, len(f), slope, intercept, rsq])
            except:
                pass
        assert len(frequencies_single_qm_all_mm) == len(intensities_single_qm_all_mm)
        # ax.scatter(frequencies_single_qm_all_mm,
        #            intensities_single_qm_all_mm,
        #            marker=markers[n_qm],
        #            label=labels[n_qm],
        #            color=colors[n_qm])
        ax.scatter(geometries_single_qm_all_mm,
                   intensities_single_qm_all_mm,
                   marker=markers[n_qm],
                   label=labels[n_qm],
                   color=colors[n_qm])
        print('{} QM/all MM'.format(n_qm))
        slope, intercept, rsq = pprint_linregress(frequencies_single_qm_all_mm,
                                                  intensities_single_qm_all_mm)
        csvwriter.writerow([n_qm, 'all', len(frequencies_single_qm_all_mm), slope, intercept, rsq])

    assert len(frequencies_all) == len(intensities_all)
    print('all QM/all MM')
    slope, intercept, rsq = pprint_linregress(frequencies_all, intensities_all)
    csvwriter.writerow(['all', 'all', len(frequencies_all), slope, intercept, rsq])

    ax.set_ylim((0.0, 1000.0))
    y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
    ax.yaxis.set_major_formatter(y_formatter)
    ax.tick_params(direction='out')
    ax.set_xlabel(r"$\nu_{3}$ frequency (cm$^{-1}$)")
    ax.set_ylabel(r"$\nu_{3}$ intensity (km/mol)")
    ax.legend(loc='lower right',
              fancybox=True,
              framealpha=0.50,
              numpoints=1,
              scatterpoints=1)
    if args.do_condon_plots:
        fig.savefig('condon_approximation.pdf', bbox_inches='tight')

    # now add the no CT data

    if args.include_noCT:
        for n_qm in sorted(frequencies_noCT_CO2_d):
            frequencies_single_qm_all_mm = []
            intensities_single_qm_all_mm = []
            for n_mm in possible_keys:
                f = frequencies_noCT_CO2_d[n_qm][n_mm]
                i = intensities_noCT_CO2_d[n_qm][n_mm]
                assert len(f) == len(i)
                frequencies_single_qm_all_mm.extend(f)
                intensities_single_qm_all_mm.extend(i)
                frequencies_all.extend(f)
                intensities_all.extend(i)
                print('{} QM/{} MM'.format(n_qm, n_mm))
            assert len(frequencies_single_qm_all_mm) == len(intensities_single_qm_all_mm)
            ax.scatter(frequencies_single_qm_all_mm,
                       intensities_single_qm_all_mm,
                       marker=markers_noCT[n_qm],
                       label=labels_noCT[n_qm],
                       color=colors_noCT[n_qm])
            print('{} QM/all MM'.format(n_qm))

        ax.legend(loc='lower right',
                  fancybox=True,
                  framealpha=0.50,
                  numpoints=1,
                  scatterpoints=1)

        if args.do_condon_plots:
            fig.savefig('condon_approximation_noCT.pdf', bbox_inches='tight')

    csvfile.close()

    plt.close(fig)

    return


def do_result_convergence_plots(results_d,
                                name='frequency',
                                n_qm_start=0,
                                n_qm_end=2,
                                func_to_apply=lambda x: x,
                                ylabel=r"$\nu_{3}$ frequency (cm$^{-1}$)",
                                labels=None,
                                colors=None,
                                errorbars=False):

    slice_partial = partial(slice, start=n_qm_start, end=n_qm_end + 1)

    print('Doing {} convergence plots'.format(name))

    fig, ax = plt.subplots()

    for n_qm in filter(slice_partial, sorted(results_d)):
        if labels:
            print("Doing plots for {}".format(labels[n_qm]))
        else:
            print("Doing plots of some kind.")
        ticks = []
        results_single_qm_all_mm = []
        results_single_qm_all_mm_mean = []
        results_single_qm_all_mm_stdev = []
        for n_mm in possible_keys:
            results_single_qm_single_mm = [func_to_apply(x) for x in results_d[n_qm][n_mm]]
            if len(results_single_qm_single_mm) > 0:
                results_single_qm_all_mm.extend(results_single_qm_single_mm)
                ticks.append(n_mm)
                results_single_qm_all_mm_mean.append(np.mean(results_single_qm_single_mm))
                results_single_qm_all_mm_stdev.append(np.std(results_single_qm_single_mm))
        # What's a cleaner way to do this...
        if markers:
            marker = markers[n_qm]
        else:
            marker = None
        if labels:
            label = labels[n_qm]
        else:
            label = None
        if colors:
            color = colors[n_qm]
        else:
            color = None
        # Make sure the data is offset properly.
        ticks = np.array(ticks)
        ticks += n_qm
        # Undo the re-categorization of the last point as always being 256.
        ticks[-1] -= n_qm
        if errorbars:
            ax.errorbar(ticks,
                        results_single_qm_all_mm_mean,
                        yerr=results_single_qm_all_mm_stdev,
                        marker=marker,
                        label=label,
                        color=color)
        else:
            ax.plot(ticks,
                    results_single_qm_all_mm_mean,
                    marker=marker,
                    label=label,
                    color=color)


    ax.set_xscale('symlog', basex=2)
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

    y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
    ax.yaxis.set_major_formatter(y_formatter)
    ax.tick_params(direction='out')
    ax.set_xlabel("# IL pairs in solvent box")
    ax.set_ylabel(ylabel)
    ax.legend(loc='best', fancybox=True, framealpha=0.50)
    fig.savefig('{}_convergence.pdf'.format(name), bbox_inches='tight')

    ax.set_xscale('linear')

    rlim = -5
    ax.set_xticks(possible_keys[:rlim + 1])
    ax.set_xlim((possible_keys[0], possible_keys[rlim]))
    fig.savefig('{}_convergence_{}.pdf'.format(name, possible_keys[rlim]),
                bbox_inches='tight')

    plt.close(fig)

    fig, ax = plt.subplots()

    ticks = list(filter(slice_partial, sorted(results_d)))
    for n_mm in possible_keys:
        results = [np.mean([func_to_apply(x)
                            for x in results_d[n_qm][n_mm]])
                   for n_qm in ticks]
        ax.plot(ticks,
                results,
                marker='o',
                label=n_mm)

    ax.set_xticks(ticks)
    ax.set_xticklabels(ticks)
    ax.yaxis.set_major_formatter(y_formatter)
    ax.tick_params(direction='out')
    ax.set_xlabel("# QM IL pairs")
    ax.set_ylabel('mean {}'.format(ylabel))
    ax.legend(loc='best', fancybox=True, framealpha=0.50)
    fig.savefig('{}_convergence_n_qm.pdf'.format(name), bbox_inches='tight')

    plt.close(fig)

    return


def do_result_convergence_plots_gaps(results_d,
                                     name='frequency',
                                     func_to_apply=lambda x: x,
                                     ylabel=r'$\nu_{3}$ frequency (cm$^{-1}$)',
                                     symbol='\omega'):

    fig, ax = plt.subplots()

    gaps_0_1 = []
    gaps_1_2 = []
    gaps_2_3 = []
    n_mm_ticks_0_1 = []
    n_mm_ticks_1_2 = []
    n_mm_ticks_2_3 = []

    for n_mm in possible_keys:
        try:
            gap_0_1 = np.mean(func_to_apply(results_d[1][n_mm])) - np.mean(func_to_apply(results_d[0][n_mm]))
            gaps_0_1.append(gap_0_1)
            n_mm_ticks_0_1.append(n_mm)
        except:
            pass
        try:
            gap_1_2 = np.mean(func_to_apply(results_d[2][n_mm])) - np.mean(func_to_apply(results_d[1][n_mm]))
            gaps_1_2.append(gap_1_2)
            n_mm_ticks_1_2.append(n_mm)
        except:
            pass
        try:
            gap_2_3 = np.mean(func_to_apply(results_d[3][n_mm])) - np.mean(func_to_apply(results_d[2][n_mm]))
            gaps_2_3.append(gap_2_3)
            n_mm_ticks_2_3.append(n_mm)
        except:
            pass

    ax.plot(n_mm_ticks_0_1, gaps_0_1, marker='s', color='red',
            label='$\Delta {symbol}_{{1-0\,\mathrm{{QM}}}}$'.format(**locals()))
    ax.plot(n_mm_ticks_1_2, gaps_1_2, marker='p', color='green',
            label='$\Delta {symbol}_{{2-1\,\mathrm{{QM}}}}$'.format(**locals()))
    ax.plot(n_mm_ticks_2_3, gaps_2_3, marker='*', color='blue',
            label='$\Delta {symbol}_{{3-2\,\mathrm{{QM}}}}$'.format(**locals()))

    ax.set_xscale('symlog', basex=2)
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

    ax.set_ylim(ax.get_ylim()[::-1])
    y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
    ax.yaxis.set_major_formatter(y_formatter)
    ax.tick_params(direction='out')
    ax.set_xlabel('# IL pairs treated as point charges')
    ax.set_ylabel(r'difference in {}'.format(ylabel))
    # ax.set_title('gaps')
    ax.legend(loc='best', fancybox=True, framealpha=0.50)
    filename = '{}_convergence_gaps.pdf'.format(name)
    print('Saving {}'.format(filename))
    fig.savefig(filename, bbox_inches='tight')

    ax.set_xscale('linear')

    rlim = -5
    ax.set_xticks(possible_keys[:rlim + 1])
    ax.set_xlim((possible_keys[0], possible_keys[rlim]))
    filename = '{}_convergence_{}_gaps.pdf'.format(name, possible_keys[rlim])
    print('Saving {}'.format(filename))
    fig.savefig(filename, bbox_inches='tight')

    plt.close(fig)

    return


def do_result_convergence_analysis(results_d,
                                   name='frequency',
                                   n_qm_start=0,
                                   n_qm_end=2,
                                   func_to_apply=lambda x: x):

    slice_partial = partial(slice, start=n_qm_start, end=n_qm_end + 1)

    print('Doing {} convergence analysis'.format(name))

    csvfile = open('{}_convergence.csv'.format(name), 'w')
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow([
        '# QM',
        '# MM',
        '# points',
        'mean',
        'median',
        'mode',
        'min',
        'max',
        'range',
        'stdev',
    ])

    for n_qm in filter(slice_partial, sorted(results_d)):
        print("Doing analysis for {}".format(labels[n_qm]))
        results_single_qm_all_mm = []
        results_single_qm_all_mm_mean = []
        results_single_qm_all_mm_median = []
        results_single_qm_all_mm_mode = []
        results_single_qm_all_mm_min = []
        results_single_qm_all_mm_max = []
        results_single_qm_all_mm_range = []
        results_single_qm_all_mm_stdev = []
        for n_mm in possible_keys:
            results_single_qm_single_mm = [func_to_apply(x) for x in results_d[n_qm][n_mm]]
            if len(results_single_qm_single_mm) > 0:
                # print('{} QM/{} MM'.format(n_qm, n_mm))
                results_single_qm_all_mm.extend(results_single_qm_single_mm)
                results_single_qm_all_mm_mean.append(np.mean(results_single_qm_single_mm))
                results_single_qm_all_mm_median.append(np.median(results_single_qm_single_mm))
                results_single_qm_all_mm_mode.append(sps.mode(results_single_qm_single_mm)[0][0])
                results_single_qm_all_mm_min.append(min(results_single_qm_single_mm))
                results_single_qm_all_mm_max.append(max(results_single_qm_single_mm))
                # we've jut updated the last values for each of these
                # lists, so take them rather than recalculate them
                results_single_qm_all_mm_range.append(results_single_qm_all_mm_max[-1] - results_single_qm_all_mm_min[-1])
                results_single_qm_all_mm_stdev.append(np.std(results_single_qm_single_mm))
                # Write entries for each possible QM/MM number
                # combination.
                csvwriter.writerow([
                    n_qm,
                    n_mm,
                    len(results_single_qm_single_mm),
                    # same thing here with just having appended
                    # calculated values to these lists
                    results_single_qm_all_mm_mean[-1],
                    results_single_qm_all_mm_median[-1],
                    results_single_qm_all_mm_mode[-1],
                    results_single_qm_all_mm_min[-1],
                    results_single_qm_all_mm_max[-1],
                    results_single_qm_all_mm_range[-1],
                    results_single_qm_all_mm_stdev[-1],
                ])
        # Write entries for all MM values combined for a single QM
        # value.
        # print('{} QM/all MM'.format(n_qm))
        val_min = min(results_single_qm_all_mm)
        val_max = max(results_single_qm_all_mm)
        val_range = val_max - val_min
        csvwriter.writerow([
            n_qm,
            'all',
            len(results_single_qm_all_mm),
            np.mean(results_single_qm_all_mm),
            np.median(results_single_qm_all_mm),
            sps.mode(results_single_qm_all_mm)[0][0],
            val_min,
            val_max,
            val_range,
            np.std(results_single_qm_all_mm),
        ])

    csvfile.close()

    return


def plot_single_snapshot_results(snapnum,
                                 snapnums_results_d,
                                 results_d,
                                 name='frequency',
                                 func_to_apply=lambda x: x,
                                 ylabel=r"$\nu_{3}$ frequency (cm$^{-1}$)",
                                 inp_fig=None,
                                 inp_ax=None,
                                 do_manip_fig=True,
                                 do_manip_ax=True):

    results_snap_d = get_single_snapshot_results(snapnum, snapnums_results_d, results_d)

    fig, ax = plt.subplots()
    if inp_fig:
        fig = inp_fig
    if inp_ax:
        ax = inp_ax

    for n_qm in sorted(results_snap_d):
        ticks = []
        results = []
        for n_mm in possible_keys:
            if len(results_snap_d[n_qm][n_mm]) > 0:
                if n_mm + n_qm >= 256:
                    ticks.append(256)
                else:
                    ticks.append(n_mm + n_qm)
                results.append(func_to_apply(results_snap_d[n_qm][n_mm][0]))

        ax.plot(ticks,
                results,
                marker=markers[n_qm],
                label=labels[n_qm],
                color=colors[n_qm])

    if do_manip_ax:
        ax.set_xscale('symlog', basex=2)
        ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
        y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(y_formatter)
        ax.tick_params(direction='out')
        ax.set_xlabel("total # of IL pairs included")
        ax.set_ylabel(ylabel)
        # ax.set_title("snapshot {}".format(snapnum))
        ax.legend(loc='best', fancybox=True, framealpha=0.50)
    if do_manip_fig:
        filename = '{}_convergence_snap{}.pdf'.format(name, snapnum)
        print('Saving {}'.format(filename))
        fig.savefig(filename, bbox_inches='tight')

    if do_manip_ax:
        ax.set_xscale('linear')
        rlim = -5
        ax.set_xticks(possible_keys[:rlim + 1])
        ax.set_xlim((possible_keys[0], possible_keys[rlim]))
    if do_manip_fig:
        filename = '{}_convergence_snap{}_{}.pdf'.format(name, snapnum, possible_keys[rlim])
        print('Saving {}'.format(filename))
        fig.savefig(filename, bbox_inches='tight')

        plt.close(fig)

    return


def plot_single_snapshot_results_qm_gaps(snapnum,
                                         snapnums_results_d,
                                         results_d,
                                         name='frequency',
                                         func_to_apply=lambda x: x,
                                         ylabel=r'$\nu_{3}$ frequency (cm$^{-1}$)',
                                         symbol='\omega'):

    results_snap_d = get_single_snapshot_results(snapnum,
                                                 snapnums_results_d,
                                                 results_d)

    fig, ax = plt.subplots()

    gaps_0_1 = []
    gaps_1_2 = []
    gaps_2_3 = []
    n_mm_ticks_0_1 = []
    n_mm_ticks_1_2 = []
    n_mm_ticks_2_3 = []

    for n_mm in possible_keys:
        try:
            gap_0_1 = func_to_apply(results_snap_d[1][n_mm][0]) - func_to_apply(results_snap_d[0][n_mm][0])
            gaps_0_1.append(gap_0_1)
            n_mm_ticks_0_1.append(n_mm)
        except:
            pass
        try:
            gap_1_2 = func_to_apply(results_snap_d[2][n_mm][0]) - func_to_apply(results_snap_d[1][n_mm][0])
            gaps_1_2.append(gap_1_2)
            n_mm_ticks_1_2.append(n_mm)
        except:
            pass
        try:
            gap_2_3 = func_to_apply(results_snap_d[3][n_mm][0]) - func_to_apply(results_snap_d[2][n_mm][0])
            gaps_2_3.append(gap_2_3)
            n_mm_ticks_2_3.append(n_mm)
        except:
            pass

    ax.plot(n_mm_ticks_0_1, gaps_0_1, marker='s', color='red',
            label='$\Delta {symbol}_{{1-0\,\mathrm{{QM}}}}$'.format(**locals()))
    ax.plot(n_mm_ticks_1_2, gaps_1_2, marker='p', color='green',
            label='$\Delta {symbol}_{{2-1\,\mathrm{{QM}}}}$'.format(**locals()))
    ax.plot(n_mm_ticks_2_3, gaps_2_3, marker='*', color='blue',
            label='$\Delta {symbol}_{{3-2\,\mathrm{{QM}}}}$'.format(**locals()))

    ax.set_xscale('symlog', basex=2)
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

    ax.set_ylim(ax.get_ylim()[::-1])
    y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
    ax.yaxis.set_major_formatter(y_formatter)
    ax.tick_params(direction='out')
    ax.set_xlabel('# IL pairs treated as point charges')
    ax.set_ylabel(r'difference in {}'.format(ylabel))
    # ax.set_title('snapshot {} gaps'.format(snapnum))
    ax.legend(loc='best', fancybox=True, framealpha=0.50)
    filename = '{}_convergence_snap{}_gaps.pdf'.format(name, snapnum)
    print('Saving {}'.format(filename))
    fig.savefig(filename, bbox_inches='tight')

    ax.set_xscale('linear')

    rlim = -5
    ax.set_xticks(possible_keys[:rlim + 1])
    ax.set_xlim((possible_keys[0], possible_keys[rlim]))
    filename = '{}_convergence_snap{}_{}_gaps.pdf'.format(name, snapnum, possible_keys[rlim])
    print('Saving {}'.format(filename))
    fig.savefig(filename, bbox_inches='tight')

    plt.close(fig)

    return


def getargs():

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--mpl-usetex", action="store_true")
    parser.add_argument("--do-condon-plots", action="store_true")
    parser.add_argument("--do-snapshot-plots", action="store_true")
    parser.add_argument("--include-noCT", action="store_true")

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = getargs()

    import matplotlib as mpl
    if args.mpl_usetex:
        mpl.rc(usetex=True)
    mpl.use("Agg")
    import matplotlib.pyplot as plt

    # Read in the pickle files that contain all the raw data.
    with open('frequencies.pypickle', 'rb') as picklefile:
        frequencies_CO2_d = pickle.load(picklefile)
    with open('intensities.pypickle', 'rb') as picklefile:
        intensities_CO2_d = pickle.load(picklefile)
    with open('frequencies_noCT.pypickle', 'rb') as picklefile:
        frequencies_noCT_CO2_d = pickle.load(picklefile)
    with open('intensities_noCT.pypickle', 'rb') as picklefile:
        intensities_noCT_CO2_d = pickle.load(picklefile)
    with open('dipoles.pypickle', 'rb') as picklefile:
        dipoles_d = pickle.load(picklefile)
    with open('geometries.pypickle', 'rb') as picklefile:
        geometries_d = pickle.load(picklefile)
    with open('snapnums_frequencies.pypickle', 'rb') as picklefile:
        snapnums_frequencies_d = pickle.load(picklefile)
    with open('snapnums_frequencies_noCT.pypickle', 'rb') as picklefile:
        snapnums_frequencies_noCT_d = pickle.load(picklefile)
    with open('snapnums_dipoles.pypickle', 'rb') as picklefile:
        snapnums_dipoles_d = pickle.load(picklefile)
    with open('snapnums_geometries.pypickle', 'rb') as picklefile:
        snapnums_geometries_d = pickle.load(picklefile)

    # Until I come up with a better idea, here's where I mangle some
    # of the keys (253, 254, 255, 256) into 256.

    # Make a copy beforehand, just in case...
    frequencies_CO2_d_unmangled = deepcopy(frequencies_CO2_d)
    intensities_CO2_d_unmangled = deepcopy(intensities_CO2_d)
    frequencies_noCT_CO2_d_unmangled = deepcopy(frequencies_noCT_CO2_d)
    intensities_noCT_CO2_d_unmangled = deepcopy(intensities_noCT_CO2_d)
    dipoles_d_unmangled = deepcopy(dipoles_d)
    geometries_d_unmangled = deepcopy(geometries_d)
    snapnums_frequencies_d_unmangled = deepcopy(snapnums_frequencies_d)
    snapnums_frequencies_noCT_d_unmangled = deepcopy(snapnums_frequencies_noCT_d)
    snapnums_dipoles_d_unmangled = deepcopy(snapnums_dipoles_d)
    snapnums_geometries_d_unmangled = deepcopy(snapnums_geometries_d)

    # Do the mangling.
    frequencies_CO2_d = mangle_dict_keys(frequencies_CO2_d)
    intensities_CO2_d = mangle_dict_keys(intensities_CO2_d)
    frequencies_noCT_CO2_d = mangle_dict_keys(frequencies_noCT_CO2_d)
    intensities_noCT_CO2_d = mangle_dict_keys(intensities_noCT_CO2_d)
    dipoles_d = mangle_dict_keys(dipoles_d)
    geometries_d = mangle_dict_keys(geometries_d)
    snapnums_frequencies_d = mangle_dict_keys(snapnums_frequencies_d)
    snapnums_frequencies_noCT_d = mangle_dict_keys(snapnums_frequencies_noCT_d)
    snapnums_dipoles_d = mangle_dict_keys(snapnums_dipoles_d)

    possible_keys = list(range(0, 18, 2)) + [32, 64, 128, 256]

    markers = [
        'o',
        's',
        'D',
        '*',
    ]
    markers_noCT = markers
    labels = [
        '0 QM pairs',
        '1 QM pair',
        '2 QM pairs',
        '3 QM pairs',
    ]
    labels_noCT = [
        '',
        '1 QM pair (no CT)',
        '2 QM pair (no CT)',
        '3 QM pair (no CT)',
    ]
    colors = [
        'black',
        'red',
        'green',
        'blue',
    ]
    colors_noCT = [
        '',
        'orange',
        'lime',
        'cyan',
    ]

    ###################################

    # Do some simple statistical analysis on the data sets and dump
    # them to CSV files.

    do_result_convergence_analysis(frequencies_CO2_d,
                                   name='frequency',
                                   n_qm_start=0,
                                   n_qm_end=3)
    # do_result_convergence_analysis(intensities_CO2_d,
    #                                name='intensity',
    #                                n_qm_start=0,
    #                                n_qm_end=2)
    # do_result_convergence_analysis(dipoles_d,
    #                                name='dipole',
    #                                n_qm_start=1,
    #                                n_qm_end=2,
    #                                func_to_apply=npl.norm)
    # if args.include_noCT:
    #     do_result_convergence_analysis(frequencies_noCT_CO2_d,
    #                                    name='frequency_noCT',
    #                                    n_qm_start=1,
    #                                    n_qm_end=2)
    #     do_result_convergence_analysis(intensities_noCT_CO2_d,
    #                                    name='intensity_noCT',
    #                                    n_qm_start=1,
    #                                    n_qm_end=2)

    ###################################

    # plots!

    # do_result_convergence_plots(frequencies_CO2_d,
    #                             name='frequency',
    #                             n_qm_start=0,
    #                             n_qm_end=3,
    #                             ylabel=r"$\nu_{3}$ frequency (cm$^{-1}$)",
    #                             labels=labels,
    #                             colors=colors)
    # do_result_convergence_plots_gaps(frequencies_CO2_d,
    #                                  name='frequency',
    #                                  func_to_apply=lambda x: x,
    #                                  ylabel=r'$\nu_{3}$ frequency (cm$^{-1}$)',
    #                                  symbol='\omega')
    # do_result_convergence_plots(intensities_CO2_d,
    #                             name='intensity',
    #                             n_qm_start=0,
    #                             n_qm_end=3,
    #                             ylabel=r"$\nu_{3}$ intensity (cm$^{-1}$)",
    #                             labels=labels,
    #                             colors=colors)
    # do_result_convergence_plots(dipoles_d,
    #                             name='dipole',
    #                             n_qm_start=1,
    #                             n_qm_end=3,
    #                             ylabel='total dipole moment (Debye)',
    #                             func_to_apply=npl.norm,
    #                             labels=labels,
    #                             colors=colors)
    # do_result_convergence_plots(dipoles_d,
    #                             name='dipole_0qm',
    #                             n_qm_start=0,
    #                             n_qm_end=0,
    #                             ylabel='total dipole moment (Debye)',
    #                             func_to_apply=npl.norm,
    #                             labels=labels,
    #                             colors=colors)
    # if args.include_noCT:
    #     do_result_convergence_plots(frequencies_noCT_CO2_d,
    #                                 name='frequency_noCT',
    #                                 n_qm_start=1, n_qm_end=2,
    #                                 ylabel=r"$\nu_{3}$ frequency (cm$^{-1}$)",
    #                                 labels=labels_noCT,
    #                                 colors=colors_noCT)
    #     do_result_convergence_plots(intensities_noCT_CO2_d,
    #                                 name='intensity_noCT',
    #                                 n_qm_start=1,
    #                                 n_qm_end=2,
    #                                 ylabel=r"$\nu_{3}$ intensity (cm$^{-1}$)",
    #                                 labels=labels_noCT,
    #                                 colors=colors_noCT)

    condon()

    # Read in the most "restrictive" set of snapshot numbers; this
    # will let us compare sets of equal size.
    snapnums = read_snapshot_file("/home/eric/Chemistry/calc.sgr/droplets/inputs_freq/representative_snapshots_3qm")
    filter_snapshots(snapnums, snapnums_frequencies_d, frequencies_CO2_d)

    do_result_convergence_plots(frequencies_CO2_d,
                                name='frequency_same_set',
                                n_qm_start=0,
                                n_qm_end=3,
                                ylabel=r"$\nu_{3}$ frequency (cm$^{-1}$)",
                                labels=labels,
                                colors=colors)
    do_result_convergence_plots(frequencies_CO2_d,
                                name='frequency_same_set_2QM',
                                n_qm_start=0,
                                n_qm_end=2,
                                ylabel=r"$\nu_{3}$ frequency (cm$^{-1}$)",
                                labels=labels,
                                colors=colors)
    do_result_convergence_plots_gaps(frequencies_CO2_d,
                                     name='frequency_same_set',
                                     func_to_apply=lambda x: x,
                                     ylabel=r'$\nu_{3}$ frequency (cm$^{-1}$)',
                                     symbol='\omega')

    if args.do_snapshot_plots:

        print('snapshot numbers:', snapnums)
        for snapnum in snapnums:
            plot_single_snapshot_results(snapnum,
                                         snapnums_frequencies_d,
                                         frequencies_CO2_d,
                                         name='frequency',
                                         func_to_apply=lambda x: x,
                                         ylabel=r'$\nu_{3}$ frequency (cm$^{-1}$)')
            plot_single_snapshot_results_qm_gaps(snapnum,
                                                 snapnums_frequencies_d,
                                                 frequencies_CO2_d,
                                                 name='frequency',
                                                 func_to_apply=lambda x: x,
                                                 ylabel=r'$\nu_{3}$ frequency (cm$^{-1}$)',
                                                 symbol='\omega')
            # plot_single_snapshot_results(snapnum,
            #                              snapnums_frequencies_d,
            #                              intensities_CO2_d,
            #                              name='intensity',
            #                              func_to_apply=lambda x: x,
            #                              ylabel=r'$\nu_{3}$ intensity (km/mol)')
            # plot_single_snapshot_results_qm_gaps(snapnum,
            #                                      snapnums_frequencies_d,
            #                                      intensities_CO2_d,
            #                                      name='intensity',
            #                                      func_to_apply=lambda x: x,
            #                                      ylabel=r'$\nu_{3}$ intensity (km/mol)',
            #                                      symbol='I')
            # plot_single_snapshot_results(snapnum,
            #                              snapnums_dipoles_d,
            #                              dipoles_d,
            #                              name='dipole',
            #                              func_to_apply=npl.norm,
            #                              ylabel='total dipole moment (Debye)')
            # plot_single_snapshot_results_qm_gaps(snapnum,
            #                                      snapnums_dipoles_d,
            #                                      dipoles_d,
            #                                      name='dipole',
            #                                      func_to_apply=npl.norm,
            #                                      ylabel='total dipole moment (Debye)',
            #                                      symbol='\mu')

    ###################################
