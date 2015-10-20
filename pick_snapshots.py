#!/usr/bin/env python3

from __future__ import print_function
from __future__ import division

import sys
import math
import pickle

import numpy as np
import scipy.stats as sps

from parse_outputs_snapshot_method_dependence import methods
from parse_outputs_snapshot_method_dependence import basis_sets
from analysis_snapshot_method_dependence import linestyles
from analysis_snapshot_method_dependence import colors

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt


def getargs():
    """Get command-line arguments."""

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--csv', action='store_true')

    parser.add_argument('--sample-method', default='b3lyp')
    parser.add_argument('--sample-basis-set', default='6-31gdp')
    parser.add_argument('--numbins', type=int, default=5)

    parser.add_argument('--do-file-operations', action='store_true')

    parser.add_argument('--dir-droplets',
                        default='/home/eric/Chemistry/calc.sgr/droplets/inputs_freq',
                        help="""This is the root directory to look for all files.""")

    args = parser.parse_args()

    return args


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
        masked_frequencies = frequencies[mask]
        binned_snapshots.append(masked_snapnums)

    return binned_snapshots


def do_sanity_check():
    """Do a sanity check and ..."""

    for ts in total_samples:
        print([[map_snapshot_to_frequency[sn] for sn in bs]
               for bs in chosen_snapshots[ts]])

    return


def do_file_operations(chosen_snapshots):
    """For all the snapshots that have been randomly chosen,

    1. There is a correspondence between the # of QM pairs and number
    of snapshots we've sampled.

    2. Make directories for each # of QM pairs, with subdirectories
    for each bin.

    3. Copy over all files from the original dirs that match a glob
    for the snapshot number and the # of QM pairs; this catches all
    the possible # MM pair matches.
    """

    import shutil

    from scripts.vmd_templates import pad_left_zeros


    for (n_qm, num_snapshots) in map_n_qm_to_num_snapshots.items():
        os.chdir(dir_droplets)
        dir_main = 'inputs_freq_{}qm'.format(n_qm)
        dir_random = 'inputs_freq_{}qm_random'.format(n_qm)
        for snapbinidx, snapbin in enumerate(chosen_snapshots[num_snapshots], start=1):
            dir_bin = os.path.join(dir_random, 'bin_{}'.format(snapbinidx))
            # If this fails, don't make things more confusing by
            # copying new files over.
            try:
                os.makedirs(dir_bin, exist_ok=False)
            except OSError:
                sys.exit(1)
            for sn in snapbin:
                inputfilename_stub = 'drop_{}_{}qm*'.format(pad_left_zeros(sn, 4), n_qm)
                inputfilename_matches = glob(os.path.join(dir_main, inputfilename_stub))
                for inputfilename in inputfilename_matches:
                    shutil.copy2(inputfilename, dir_bin)


if __name__ == '__main__':

    args = getargs()

    # Read in the pickle files that contain all the raw data.
    with open('frequencies.pypickle', 'rb') as picklefile:
        frequencies_d = pickle.load(picklefile)
    with open('snapnums_frequencies.pypickle', 'rb') as picklefile:
        snapnums_d = pickle.load(picklefile)

    # Write the means and standard deviations to a CSV file if asked.
    if args.csv:
        import csv
        csvfh = open('pick_snapshots_means_stdevs.csv', 'w')
        csvwriter = csv.writer(csvfh)
        csvwriter.writerow(['method', 'basis set', 'mean', 'stdev (pop)'])

    fig, ax = plt.subplots()
    fig_combined_hists, ax_combined_hists = plt.subplots()

    for basis_set in basis_sets:
        for method in methods:

            pairs = sorted(zip(snapnums_d[method][basis_set][0][0],
                               frequencies_d[method][basis_set][0][0]),
                           key=lambda x: x[1])
            snapnums = [x[0] for x in pairs]
            frequencies = [x[1] for x in pairs]

            mmin = min(frequencies)
            mmax = max(frequencies)
            rng = mmax - mmin
            mean = np.mean(frequencies)
            median = np.median(frequencies)
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
            print(' mode          : {}'.format(mode))
            print(' stdev (pop)   : {:.3f}'.format(stdev_pop))
            print(' stdev (sample): {:.3f}'.format(stdev_sample))

            if args.csv:
                csvwriter.writerow([methods[method],
                                    basis_sets[basis_set],
                                    mean,
                                    stdev_pop])

            numbins = args.numbins
            center = mean
            binwidth = stdev_pop
            bins = make_bin_edges(numbins, center, binwidth)

            print(' bins: {}'.format(bins))

            hist, bin_edges = np.histogram(frequencies, bins)
            weights = hist / sum(hist)

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
                         label='probability distribution function',
                         color='red',
                         linewidth=3,
                         linestyle='--')

            ax_hist.legend(fancybox=True, loc='upper right', framealpha=0.50)
            fig_hist.savefig('hist_{}_{}_{}.pdf'.format(method, basis_set, numbins), bbox_inches='tight')

            plt.close(fig_hist)

            # NumPy and matplotlib are doing the same thing for the bins.
            assert bin_edges.all() == bins_mpl.all()

            ax.plot(frequencies,
                    label='{}/{}'.format(methods[method], basis_sets[basis_set]),
                    color=colors[method],
                    linestyle=linestyles[basis_set])

            ax_combined_hists.plot(linspace_bins,
                                   pdf_bins,
                                   label='{}/{}'.format(methods[method], basis_sets[basis_set]),
                                   color=colors[method],
                                   linewidth=2,
                                   linestyle=linestyles[basis_set])

        ax.legend(fancybox=True,
                  loc='best',
                  framealpha=0.50,
                  fontsize='small')
        fig.savefig('snapshots_ordered.pdf', bbox_inches='tight')

        plt.close(fig)

        # Add the experimental xxx/xxx frequency as a stick.
        # Is there a way to avoid hard-coding the "intensity" and have
        # it determined automatically?
        ax_combined_hists.stem([2340.0], [0.004],
                               linefmt='k:',
                               markerfmt='k:')

        ax_combined_hists.set_xlabel(r'$\nu_3$ frequency (cm$^{-1}$)')
        ax_combined_hists.set_ylabel('arbitrary units')
        ax_combined_hists.set_title('probability density fnctions')
        ax_combined_hists.legend(fancybox=True,
                                 loc='best',
                                 framealpha=0.50,
                                 fontsize='small')
        fig_combined_hists.savefig('combined_hists.pdf', bbox_inches='tight')

        plt.close(fig_combined_hists)

    ### Redo the histogram for a specific method/basis set and sample
    ### from each bin.
    method = args.sample_method
    basis_set = args.sample_basis_set
    print('=' * 78)
    print('****** SAMPLING')
    print('*** method   : {}'.format(methods[method]))
    print('*** basis set: {}'.format(basis_sets[basis_set]))
    pairs = sorted(zip(snapnums_d[method][basis_set][0][0],
                       frequencies_d[method][basis_set][0][0]),
                   key=lambda x: x[1])
    snapnums = np.array([x[0] for x in pairs])
    frequencies = np.array([x[1] for x in pairs])
    map_snapshot_to_frequency = {p[0] : p[1] for p in pairs}
    mmin = min(frequencies)
    mmax = max(frequencies)
    rng = mmax - mmin
    mean = np.mean(frequencies)
    median = np.median(frequencies)
    mode = sps.mode(frequencies)
    stdev_pop = np.std(frequencies, ddof=1)
    stdev_sample = np.std(frequencies, ddof=0)
    numbins = args.numbins
    center = mean
    binwidth = stdev_pop
    bin_edges = make_bin_edges(numbins, center, binwidth)
    hist, _ = np.histogram(frequencies, bin_edges)
    assert bin_edges.all() == _.all()
    hist_density, _ = np.histogram(frequencies, bin_edges, density=True)
    assert bin_edges.all() == _.all()
    print('bin edges:')
    print(bin_edges)
    print('histogram:')
    print(hist)
    weights = hist / sum(hist)
    print('weights:')
    print(weights)
    assert (hist_density * stdev_pop).all() == weights.all()
    bin_indices = np.digitize(frequencies, bin_edges)
    total_samples = np.array([100, 25, 10])
    samples_per_bin = total_samples / numbins
    binned_snapshots = bin_snapshots(frequencies, snapnums, numbins, center, binwidth)

    chosen_snapshots = {ts : [] for ts in total_samples}

    # Pre-populate the first set of chosen snapshots, since this is
    # the only time we choose directly from the binned snapshots.

    for bs in binned_snapshots:
        choices = np.array(sorted(np.random.choice(bs, size=samples_per_bin[0], replace=False)))
        chosen_snapshots[total_samples[0]].append(choices)

    # All other sets of chosen snapshots will pull from the prevous
    # set, so that every reduction in total samples/samples per bin is
    # a subset.

    for tsidx in range(1, len(total_samples)):
        for bs in chosen_snapshots[total_samples[tsidx-1]]:
            choices = np.array(sorted(np.random.choice(bs, size=int(samples_per_bin[tsidx]), replace=False)))
            chosen_snapshots[total_samples[tsidx]].append(choices)

    # As a sanity check, what are the frequencies that these snapshots
    # correspond to?

    do_sanity_check()

    # These are the total number of snapshots we want to sample for
    # each number of QM pairs available; number is decreasing due to
    # increasing cost.
    map_n_qm_to_num_snapshots = {
        1: 100,
        2: 25,
        3: 10,
    }

    # First, check to see if we've already chosen any snapshots ->
    # this corresponds to having their numbers dumped to files.

    from glob import glob
    import os

    do_file_operations_p = True
    res = glob(os.path.join(args.dir_droplets, 'representative_snapshots_*qm*'))
    if len(res) > 0:
        # If we already have results on disk, read them in and write
        # over whatever we've already chosen, but we're going to
        # assume the distribution is the same.
        chosen_snapshots = dict()
        do_file_operations_p = False
        for r in res:
            if 'bin' in r:
                stub = os.path.basename(r)
                chomp = stub.split('_')
                n_qm = int(chomp[2][0])
                binnum = int(chomp[4])
                with open(r) as binfile:
                    snapnums = sorted(map(int, binfile.readlines()))
                    chosen_snapshots[(n_qm, binnum)] = snapnums

    print('(n_qm, binnum): chosen snapshot numbers:')
    for ((n_qm, binnum), snapnums) in sorted(chosen_snapshots.items()):
        print('({}, {}): {}'.format(n_qm, binnum, snapnums))

    # Make directories and copy the snapshots over.

    if do_file_operations_p:
        print('We can do file operations.')
        # Only do this if we really, really want to.
        if args.do_file_operations:
            print('Actually going to do file operations.')
            do_file_operations(chosen_snapshots)
    else:
        print("We can't do file operations.")

    if args.csv:
        csvfh.close()
