#!/usr/bin/env python3

from __future__ import print_function
from __future__ import division

import csv
import math
import pickle

import numpy as np
import scipy.stats as sps

from parse_outputs_snapshot_method_dependence import \
    (methods, basis_sets)
from analysis_snapshot_method_dependence import \
    (linestyles, colors)

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# These are the total number of snapshots we want to sample for
# each number of QM pairs available; number is decreasing due to
# increasing cost.
map_n_qm_to_num_snapshots = {
    1: 100,
    2: 25,
    3: 10,
}


def norm_pdf(x, mu, sigma):
    """Return the normal probability distribution function."""

    prefac = 1 / (sigma * math.sqrt(2 * math.pi))
    expt = np.exp(-(x - mu)**2 / (2 * sigma**2))

    return prefac * expt


def make_bin_edges(numbins, center, binwidth):
    """Given an integer number of bins, the absolute center of the bins,
    and the width of each, return an array of bin edges.
    """

    start = center - (numbins * binwidth / 2)
    stop = center + (numbins * binwidth / 2)

    bins_linspace = np.linspace(start=start, stop=stop, num=numbins + 1)

    return bins_linspace


def bin_snapshots(frequencies, snapnums, numbins, center, binwidth):
    """Given a set of frequencies and their snapshot numbers, place the
    snapshot numbers into the appropriate bins.
    """

    binned_snapshots = []

    bin_edges = make_bin_edges(numbins, center, binwidth)
    bin_indices = np.digitize(frequencies, bin_edges)

    for binidx in range(1, numbins + 1):
        mask = np.where(bin_indices == binidx)[0]
        masked_snapnums = snapnums[mask]
        # masked_frequencies = frequencies[mask]
        binned_snapshots.append(masked_snapnums)

    return binned_snapshots


def make_bin_dictionary(total_samples, binned_snapshots, samples_per_bin):

    chosen_snapshots = {ts : [] for ts in total_samples}

    # Pre-populate the first set of chosen snapshots, since this is
    # the only time we choose directly from the binned snapshots.

    for binidx, bs in enumerate(binned_snapshots, start=1):
        choices = np.array(sorted(np.random.choice(bs, size=samples_per_bin[0],
                                                   replace=False)))
        chosen_snapshots[total_samples[0]].append(choices)

    # All other sets of chosen snapshots will pull from the previous
    # set, so that every reduction in total samples/samples per bin is
    # a subset.

    for tsidx in range(1, len(total_samples)):
        for bs in chosen_snapshots[total_samples[tsidx-1]]:
            choices = np.array(sorted(np.random.choice(bs,
                                                       size=int(samples_per_bin[tsidx]),
                                                       replace=False)))
            chosen_snapshots[total_samples[tsidx]].append(choices)

    for (n_qm, num_snapshots) in map_n_qm_to_num_snapshots.items():
        print(n_qm, num_snapshots)
        for (snapbinidx, snapbin) in enumerate(chosen_snapshots[num_snapshots], start=1):
            print(snapbinidx, snapbin)

    return chosen_snapshots


def analysis_single(numbins, n_qm, n_mm):

    pairs = sorted(zip(snapnums_d[n_qm][n_mm],
                       frequencies_d[n_qm][n_mm]),
                   key=lambda x: x[1])
    print(pairs)
    # snapnums = [x[0] for x in pairs]
    frequencies = [x[1] for x in pairs]

    mmin = min(frequencies)
    mmax = max(frequencies)
    rng = mmax - mmin
    mean = np.mean(frequencies)
    median = np.median(frequencies)
    mode = sps.mode(frequencies)
    stdev_pop = np.std(frequencies, ddof=1)
    stdev_sample = np.std(frequencies, ddof=0)

    print(' min           : {:.2f}'.format(mmin))
    print(' max           : {:.2f}'.format(mmax))
    print(' range         : {:.2f}'.format(rng))
    print(' mean          : {:.3f}'.format(mean))
    print(' median        : {:.3f}'.format(median))
    print(' mode          : {}'.format(mode))
    print(' stdev (pop)   : {:.3f}'.format(stdev_pop))
    print(' stdev (sample): {:.3f}'.format(stdev_sample))

    center = mean
    binwidth = stdev_pop
    bins = make_bin_edges(numbins, center, binwidth)
    hist, bin_edges = np.histogram(frequencies, bins)
    weights = hist / sum(hist)

    print(' bins     : {}'.format(bins))
    print(' histogram: {}'.format(hist))
    print(' weights  : {}'.format(weights))

    # If you specify the bins, the edges that are returned are
    # identical to the input.
    assert bins.all() == bin_edges.all()

    fig_hist, ax_hist = plt.subplots()

    n, bins_mpl, _ = ax_hist.hist(frequencies, bins=bin_edges, normed=True)

    linspace_bins = np.linspace(bins[0], bins[-1], 300)
    linspace_all = np.linspace(mmin, mmax, 300)
    pdf_bins = norm_pdf(linspace_bins, mu=center, sigma=binwidth)
    pdf_all = norm_pdf(linspace_all, mu=center, sigma=binwidth)
    print('sum(pdf_bins): {}'.format(sum(pdf_bins)))
    print('sum(pdf_all) : {}'.format(sum(pdf_all)))
    ax_hist.plot(linspace_bins,
                 pdf_bins,
                 label='probability density function',
                 color='red',
                 linewidth=3,
                 linestyle='--')

    # Add the experimental BMIM/PF6 frequency as a stick.
    ax_hist.stem([2340.0],
                 [ax_hist.get_ylim()[1]],
                 linefmt='k:',
                 markerfmt='k:',
                 label='experiment')

    ax_hist.legend(fancybox=True, loc='upper right', framealpha=0.50)
    filename = 'hist_{}qm_{}mm_{}.pdf'.format(n_qm, n_mm, numbins)
    fig_hist.savefig(filename, bbox_inches='tight')
    print('Saving {}'.format(filename))

    plt.close(fig_hist)

    # NumPy and matplotlib are doing the same thing for the bins.
    assert bin_edges.all() == bins_mpl.all()

    return


def analysis_all_methods_all_basis_sets(numbins, n_qm=0, n_mm=0):

    if args.csv:
        csvfh = open('analysis_all_methods_all_basis_sets_{}qm_{}mm.csv'.format(n_qm, n_mm), 'w')
        csvwriter = csv.writer(csvfh)
        csvwriter.writerow([
            'method',
            'basis set',
            'min',
            'max',
            'range',
            'mean',
            'median',
            'mean-median',
            'stdev (pop)',
            'stdev (sample)',
            'bins',
            'histogram',
            'weights',
            'sum(histogram)',
        ])

    fig_snapshots_ordered, ax_snapshots_ordered = plt.subplots()
    fig_combined_hists, ax_combined_hists = plt.subplots()

    cmap = cm.get_cmap('viridis')

    num_plot_lines = len(methods)

    # Sort by increasing RI-MP2/cc-pVTZ frequency.
    idxsort = np.argsort(frequencies_d['ri-mp2']['cc-pvtz'][n_qm][n_mm])

    for basis_set in basis_sets:
        for idx, method in enumerate(methods):

            # pairs = sorted(zip(snapnums_d[method][basis_set][n_qm][n_mm],
            #                    frequencies_d[method][basis_set][n_qm][n_mm]),
            #                key=lambda x: x[1])
            # frequencies = [x[1] for x in pairs]
            frequencies = np.array(frequencies_d[method][basis_set][n_qm][n_mm])[idxsort]

            mmin = min(frequencies)
            mmax = max(frequencies)
            rng = mmax - mmin
            mean = np.mean(frequencies)
            median = np.median(frequencies)
            mean_median_diff = mean - median
            mode = sps.mode(frequencies)
            stdev_pop = np.std(frequencies, ddof=1)
            stdev_sample = np.std(frequencies, ddof=0)

            print('=' * 78)
            print('*** method   : {}'.format(methods[method]))
            print('*** basis set: {}'.format(basis_sets[basis_set]))
            print(' min           : {:.2f}'.format(mmin))
            print(' max           : {:.2f}'.format(mmax))
            print(' range         : {:.2f}'.format(rng))
            print(' mean          : {:.3f}'.format(mean))
            print(' median        : {:.3f}'.format(median))
            print(' mean-median   : {:.3f}'.format(mean_median_diff))
            print(' mode          : {}'.format(mode))
            print(' stdev (pop)   : {:.3f}'.format(stdev_pop))
            print(' stdev (sample): {:.3f}'.format(stdev_sample))

            center = mean
            binwidth = stdev_pop
            bins = make_bin_edges(numbins, center, binwidth)

            print(' bins: {}'.format(bins))

            hist, bin_edges = np.histogram(frequencies, bins)
            sum_hist = sum(hist)
            weights = hist / sum_hist

            print(' histogram     : {}'.format(hist))
            print(' weights       : {}'.format(weights))
            print(' sum(histogram): {}'.format(sum_hist))

            if args.csv:
                csvwriter.writerow([
                    methods[method],
                    basis_sets[basis_set],
                    mmin,
                    mmax,
                    rng,
                    mean,
                    median,
                    mean_median_diff,
                    stdev_pop,
                    stdev_sample,
                    bins,
                    hist,
                    weights,
                    sum_hist,
                ])

            # If you specify the bins, the edges that are returned are
            # identical to the input.
            assert bins.all() == bin_edges.all()

            fig_hist, ax_hist = plt.subplots()

            n, bins_mpl, _ = ax_hist.hist(frequencies, bins=bin_edges, normed=True)
            print('n', n)
            print('weights/n', weights/n)

            linspace_bins = np.linspace(bins[0], bins[-1], 300)
            linspace_all = np.linspace(mmin, mmax, 300)
            pdf_bins = norm_pdf(linspace_bins, mu=center, sigma=binwidth)
            pdf_all = norm_pdf(linspace_all, mu=center, sigma=binwidth)
            print('sum(pdf_bins): {}'.format(sum(pdf_bins)))
            print('sum(pdf_all) : {}'.format(sum(pdf_all)))

            ax_hist.plot(linspace_bins,
                         pdf_bins,
                         label='probability density function',
                         color='red',
                         linewidth=3,
                         linestyle='--')

            # Add the experimental BMIM/PF6 frequency as a stick.
            ax_hist.stem([2340.0],
                         [ax_hist.get_ylim()[1]],
                         linefmt='k:',
                         markerfmt='k:',
                         label='experiment')

            ax_hist.legend(fancybox=True, loc='upper right', framealpha=0.50)
            filename = 'hist_{}_{}_{}QM_{}MM_{}.pdf'.format(method,
                                                            basis_set,
                                                            n_qm,
                                                            n_mm,
                                                            numbins)
            fig_hist.savefig(filename, bbox_inches='tight')
            print('Saving {}'.format(filename))

            plt.close(fig_hist)

            # NumPy and matplotlib are doing the same thing for the bins.
            assert bin_edges.all() == bins_mpl.all()

            c = cmap(float(idx) / (num_plot_lines - 1))

            ax_snapshots_ordered.plot(frequencies,
                                      label='{}/{}'.format(methods[method], basis_sets[basis_set]),
                                      # color=colors[method],
                                      color=c,
                                      linestyle=linestyles[basis_set])

            ax_combined_hists.plot(linspace_bins,
                                   pdf_bins,
                                   label='{}/{}'.format(methods[method], basis_sets[basis_set]),
                                   # color=colors[method],
                                   color=c,
                                   linestyle=linestyles[basis_set])

    # Add the experimental BMIM/PF6 frequency as a line.
    ax_snapshots_ordered.plot(ax_snapshots_ordered.get_xlim(),
                              [2340.0, 2340.0],
                              marker='',
                              linestyle=':',
                              color='black',
                              label='experiment')

    # Having the snapshot number start at 0 looks bad and doesn't
    # make sense.
    # xticks = ax_snapshots_ordered.get_xticks()
    # xticks[0] = 1
    # ax_snapshots_ordered.set_xticks(xticks)

    # This removes all the tick labels. Ticks themselves are turned off later.
    # ax_snapshots_ordered.set_xticklabels([])

    # Only place ticks + labels at 1 and 100.
    ticks = (1, 100)
    ax_snapshots_ordered.set_xticks(ticks)
    ax_snapshots_ordered.set_xticklabels(ticks)

    ax_snapshots_ordered.legend(fancybox=True,
              loc='best',
              framealpha=0.50,
              fontsize='small')
    ax_snapshots_ordered.set_xlabel('Sorted Geometry #',
                  fontsize='large')
    ax_snapshots_ordered.set_ylabel(r'$\nu_3$ harmonic frequency (cm$^{-1}$)',
                                    fontsize='large')
    ax_snapshots_ordered.tick_params(direction='out',
                                     top='off')
    filename = 'snapshots_ordered_{}QM_{}MM.pdf'.format(n_qm, n_mm)
    fig_snapshots_ordered.savefig(filename, bbox_inches='tight')
    print('Saving {}'.format(filename))

    # for line in ax_snapshots_ordered.get_lines():
    #     print(len(line.get_ydata()), line.get_ydata())

    plt.close(fig_snapshots_ordered)

    # Add the experimental BMIM/PF6 frequency as a line.
    ax_combined_hists.plot([2340.0, 2340.0],
                           ax_combined_hists.get_ylim(),
                           marker='',
                           linestyle=':',
                           color='black',
                           label='experiment')

    ax_combined_hists.legend(fancybox=True,
                             loc='best',
                             framealpha=0.50,
                             fontsize='small')
    ax_combined_hists.set_xlabel(r'$\nu_3$ harmonic frequency (cm$^{-1}$)',
                                 fontsize='large')
    # ax_combined_hists.set_ylabel('arbitrary units',
    #                              fontsize='large')
    ax_combined_hists.set_ylabel('probability density function',
                                 fontsize='large')
    # ax_combined_hists.set_title('probability density functions')
    ax_combined_hists.tick_params(direction='out')
    filename = 'combined_pdfs_{}QM_{}MM.pdf'.format(n_qm, n_mm)
    fig_combined_hists.savefig(filename, bbox_inches='tight')
    print('Saving {}'.format(filename))

    # for line in ax_combined_hists.get_lines():
    #     print(len(line.get_ydata()), line.get_ydata())

    plt.close(fig_combined_hists)

    if args.csv:
        csvfh.close()

    return


def analysis_difference_mm(numbins):
    """Plot the difference between the 0 QM/0 MM and 0 QM/256 MM data."""

    if args.csv:
        csvfh = open('analysis_difference_mm.csv', 'w')
        csvwriter = csv.writer(csvfh)
        csvwriter.writerow([
            'method',
            'min',
            'max',
            'range',
            'mean',
            'median',
            'mean-median',
            'stdev (pop)',
            'stdev (sample)',
            'bins',
            'histogram',
            'weights',
            'sum(histogram)',
        ])

    cmap = cm.get_cmap('viridis')

    num_plot_lines = len(methods)

    fig_pdf, ax_pdf = plt.subplots()
    fig_ordered, ax_ordered = plt.subplots()

    for basis_set in basis_sets:
        for idx, method in enumerate(methods):

            pairs_0 = sorted(zip(snapnums_d[method][basis_set][0][0],
                                 frequencies_d[method][basis_set][0][0]),
                             key=lambda x: x[1])
            pairs_256 = sorted(zip(snapnums_d[method][basis_set][0][256],
                                   frequencies_d[method][basis_set][0][256]),
                               key=lambda x: x[1])
            frequencies_0 = np.array([x[1] for x in pairs_0])
            frequencies_256 = np.array([x[1] for x in pairs_256])
            frequencies_diff = frequencies_256 - frequencies_0

            mmin = min(frequencies_diff)
            mmax = max(frequencies_diff)
            rng = mmax - mmin
            mean = np.mean(frequencies_diff)
            median = np.median(frequencies_diff)
            mean_median_diff = mean - median
            mode = sps.mode(frequencies_diff)
            stdev_pop = np.std(frequencies_diff, ddof=1)
            stdev_sample = np.std(frequencies_diff, ddof=0)

            print(' min           : {:.2f}'.format(mmin))
            print(' max           : {:.2f}'.format(mmax))
            print(' range         : {:.2f}'.format(rng))
            print(' mean          : {:.3f}'.format(mean))
            print(' median        : {:.3f}'.format(median))
            print(' mean-median   : {:.3f}'.format(mean_median_diff))
            print(' mode          : {}'.format(mode))
            print(' stdev (pop)   : {:.3f}'.format(stdev_pop))
            print(' stdev (sample): {:.3f}'.format(stdev_sample))

            center_diff = np.mean(frequencies_diff)
            binwidth_diff = np.std(frequencies_diff, ddof=1)
            bins_diff = make_bin_edges(numbins, center_diff, binwidth_diff)
            hist_diff, bin_edges_diff = np.histogram(frequencies_diff, bins_diff)
            sum_hist_diff = sum(hist_diff)
            weights_diff = hist_diff / sum_hist_diff
            assert bins_diff.all() == bin_edges_diff.all()
            linspace_bins_diff = np.linspace(bins_diff[0], bins_diff[-1], 300)
            pdf_bins_diff = norm_pdf(linspace_bins_diff, mu=center_diff, sigma=binwidth_diff)

            print('{}/{}'.format(methods[method], basis_sets[basis_set]))
            print(' bins        : {}'.format(bins_diff))
            print(' histogram   : {}'.format(hist_diff))
            print(' weights     : {}'.format(weights_diff))
            print('sum(hist)    : {}'.format(sum_hist_diff))
            print('sum(pdf_bins): {}'.format(sum(pdf_bins_diff)))
            # print(frequencies_diff)

            if args.csv:
                csvwriter.writerow([
                    method,
                    mmin,
                    mmax,
                    rng,
                    mean,
                    median,
                    mean_median_diff,
                    stdev_pop,
                    stdev_sample,
                    bins_diff,
                    hist_diff,
                    weights_diff,
                    sum_hist_diff,
                ])

            c = cmap(float(idx) / (num_plot_lines - 1))

            ax_pdf.plot(linspace_bins_diff,
                        pdf_bins_diff,
                        label='{}/{}'.format(methods[method], basis_sets[basis_set]),
                        color=c,
                        linewidth=0.5,
                        linestyle=linestyles[basis_set])

            ax_ordered.plot(frequencies_diff,
                            label='{}/{}'.format(methods[method], basis_sets[basis_set]),
                            color=c,
                            linewidth=0.5,
                            linestyle=linestyles[basis_set])

    ax_pdf.set_xlabel(r'$\nu_3$ harmonic frequency (cm$^{-1}$)', fontsize='large')
    ax_pdf.set_ylabel('probability density function', fontsize='large')
    ax_ordered.set_xticklabels([])
    ax_ordered.set_ylabel(r'$\nu_3$ harmonic frequency (cm$^{-1}$)', fontsize='large')
    ax_pdf.tick_params(direction='out')
    ax_ordered.tick_params(direction='out')
    ax_pdf.legend(fancybox=True, loc='best', framealpha=0.50, fontsize='small')
    ax_ordered.legend(fancybox=True, loc='best', framealpha=0.50, fontsize='small')
    filename_pdf = 'difference_0qm_256mm_0mm_pdf.pdf'
    filename_ordered = 'difference_0qm_256mm_0mm_ordered.pdf'
    print('Saving {}'.format(filename_pdf))
    fig_pdf.savefig(filename_pdf, bbox_inches='tight')
    print('Saving {}'.format(filename_ordered))
    fig_ordered.savefig(filename_ordered, bbox_inches='tight')

    plt.close('all')

    if args.csv:
        csvfh.close()

    return


def analysis_difference_basis(numbins):
    """Plot the difference between the cc-pVTZ and 6-31G(d,p) basis sets
    for 0 QM/256 MM.
    """

    if args.csv:
        csvfh = open('analysis_difference_basis.csv', 'w')
        csvwriter = csv.writer(csvfh)
        csvwriter.writerow([
            'method',
            'min',
            'max',
            'range',
            'mean',
            'median',
            'mean-median',
            'stdev (pop)',
            'stdev (sample)',
            'bins',
            'histogram',
            'weights',
            'sum(histogram)',
        ])

    cmap = cm.get_cmap('viridis')

    num_plot_lines = len(methods)

    fig_pdf, ax_pdf = plt.subplots()
    fig_ordered, ax_ordered = plt.subplots()

    for idx, method in enumerate(methods):

        pairs_dunning = sorted(zip(snapnums_d[method]['cc-pvtz'][0][256],
                                   frequencies_d[method]['cc-pvtz'][0][256]),
                               key=lambda x: x[1])
        pairs_pople = sorted(zip(snapnums_d[method]['6-31gdp'][0][256],
                                 frequencies_d[method]['6-31gdp'][0][256]),
                             key=lambda x: x[1])
        frequencies_dunning = np.array([x[1] for x in pairs_dunning])
        frequencies_pople = np.array([x[1] for x in pairs_pople])
        frequencies_diff = frequencies_dunning - frequencies_pople

        mmin = min(frequencies_diff)
        mmax = max(frequencies_diff)
        rng = mmax - mmin
        mean = np.mean(frequencies_diff)
        median = np.median(frequencies_diff)
        mean_median_diff = mean - median
        mode = sps.mode(frequencies_diff)
        stdev_pop = np.std(frequencies_diff, ddof=1)
        stdev_sample = np.std(frequencies_diff, ddof=0)

        print(' min           : {:.2f}'.format(mmin))
        print(' max           : {:.2f}'.format(mmax))
        print(' range         : {:.2f}'.format(rng))
        print(' mean          : {:.3f}'.format(mean))
        print(' median        : {:.3f}'.format(median))
        print(' mean-median   : {:.3f}'.format(mean_median_diff))
        print(' mode          : {}'.format(mode))
        print(' stdev (pop)   : {:.3f}'.format(stdev_pop))
        print(' stdev (sample): {:.3f}'.format(stdev_sample))

        center_diff = np.mean(frequencies_diff)
        binwidth_diff = np.std(frequencies_diff, ddof=1)
        bins_diff = make_bin_edges(numbins, center_diff, binwidth_diff)
        hist_diff, bin_edges_diff = np.histogram(frequencies_diff, bins_diff)
        sum_hist_diff = sum(hist_diff)
        weights_diff = hist_diff / sum_hist_diff
        assert bins_diff.all() == bin_edges_diff.all()
        linspace_bins_diff = np.linspace(bins_diff[0], bins_diff[-1], 300)
        pdf_bins_diff = norm_pdf(linspace_bins_diff, mu=center_diff, sigma=binwidth_diff)

        print('{}'.format(methods[method]))
        print(' bins        : {}'.format(bins_diff))
        print(' histogram   : {}'.format(hist_diff))
        print(' weights     : {}'.format(weights_diff))
        print('sum(hist)    : {}'.format(sum_hist_diff))
        print('sum(pdf_bins): {}'.format(sum(pdf_bins_diff)))
        # print(frequencies_diff)

        if args.csv:
            csvwriter.writerow([
                method,
                mmin,
                mmax,
                rng,
                mean,
                median,
                mean_median_diff,
                stdev_pop,
                stdev_sample,
                bins_diff,
                hist_diff,
                weights_diff,
                sum_hist_diff,
            ])

        c = cmap(float(idx) / (num_plot_lines - 1))

        ax_pdf.plot(linspace_bins_diff,
                    pdf_bins_diff,
                    label='{}'.format(methods[method]),
                    color=c,
                    linewidth=0.5,
                    linestyle='-')

        ax_ordered.plot(frequencies_diff,
                        label='{}'.format(methods[method]),
                        color=c,
                        linewidth=0.5,
                        linestyle='-')

    ax_pdf.set_xlabel(r'$\nu_3$ harmonic frequency (cm$^{-1}$)', fontsize='large')
    ax_pdf.set_ylabel('probability density function', fontsize='large')
    ax_ordered.set_xticklabels([])
    ax_ordered.set_ylabel(r'$\nu_3$ harmonic frequency (cm$^{-1}$)', fontsize='large')
    ax_pdf.tick_params(direction='out')
    ax_ordered.tick_params(direction='out')
    ax_pdf.legend(fancybox=True, loc='best', framealpha=0.50)
    ax_ordered.legend(fancybox=True, loc='best', framealpha=0.50)
    filename_pdf = 'difference_0qm_256mm_dunning_pople_pdf.pdf'
    filename_ordered = 'difference_0qm_256mm_dunning_pople_ordered.pdf'
    print('Saving {}'.format(filename_pdf))
    fig_pdf.savefig(filename_pdf, bbox_inches='tight')
    print('Saving {}'.format(filename_ordered))
    fig_ordered.savefig(filename_ordered, bbox_inches='tight')

    plt.close('all')

    if args.csv:
        csvfh.close()

    return


def getargs():
    """Get command-line arguments."""

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('analysis_type',
                        choices=('snapshot_method_dependence',
                                 'single',
                                 'difference_mm',
                                 'difference_basis',))

    parser.add_argument('--csv', action='store_true')

    parser.add_argument('--numbins', type=int, default=5)

    parser.add_argument('--n_qm', type=int)
    parser.add_argument('--n_mm', type=int)

    args = parser.parse_args()

    return args



if __name__ == '__main__':

    args = getargs()

    # Read in the pickle files that contain all the raw data.
    with open('frequencies.pypickle', 'rb') as picklefile:
        frequencies_d = pickle.load(picklefile)
    with open('snapnums_frequencies.pypickle', 'rb') as picklefile:
        snapnums_d = pickle.load(picklefile)

    if args.analysis_type == 'snapshot_method_dependence':
        analysis_all_methods_all_basis_sets(numbins=args.numbins,
                                            n_qm=args.n_qm,
                                            n_mm=args.n_mm)
    elif args.analysis_type == 'single':
        analysis_single(numbins=args.numbins,
                        n_qm=args.n_qm,
                        n_mm=args.n_mm)
    elif args.analysis_type == 'difference_mm':
        analysis_difference_mm(numbins=args.numbins)
    elif args.analysis_type == 'difference_basis':
        analysis_difference_basis(numbins=args.numbins)
    else:
        pass
