#!/usr/bin/env python3

from __future__ import print_function
from __future__ import division

import sys
import pickle
from glob import glob
import os

import numpy as np
import scipy.stats as sps

from parse_outputs_snapshot_method_dependence import \
    (methods, basis_sets)
from analysis_histograms import \
    (make_bin_edges, bin_snapshots, make_bin_dictionary,
     map_n_qm_to_num_snapshots)


def getargs():
    """Get command-line arguments."""

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--sample-method', default='b3lyp')
    parser.add_argument('--sample-basis-set', default='6-31gdp')
    parser.add_argument('--numbins', type=int, default=5)

    parser.add_argument('--do-file-operations', action='store_true')

    parser.add_argument('--dir-droplets',
                        default='/home/eric/Chemistry/calc.sgr/paper_02_CD_SC/inputs_freq',
                        help="""This is the root directory to look for all files.""")

    args = parser.parse_args()

    return args


def do_sanity_check():
    """Do a sanity check and ...DOES THIS DO ANYTHING RIGHT NOW"""

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
        os.chdir(args.dir_droplets)
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

    return


if __name__ == '__main__':

    args = getargs()

    # Read in the pickle files that contain all the raw data.
    with open('frequencies.pypickle', 'rb') as picklefile:
        frequencies_d = pickle.load(picklefile)
    with open('snapnums_frequencies.pypickle', 'rb') as picklefile:
        snapnums_d = pickle.load(picklefile)

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

    chosen_snapshots = make_bin_dictionary(total_samples, binned_snapshots, samples_per_bin)

    # As a sanity check, what are the frequencies that these snapshots
    # correspond to?

    # do_sanity_check()

    # First, check to see if we've already chosen any snapshots ->
    # this corresponds to having their numbers dumped to files.

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
