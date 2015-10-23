#!/usr/bin/env python3

from __future__ import print_function

import os
import shutil
import glob

n_qm_vals = (1, 2, 3)

binned_snapnums = dict()
for n_qm in n_qm_vals:
    binned_snapnums[n_qm] = dict()
    binned_snapnums_filename = '/home/eric/Chemistry/calc.sgr/droplets/inputs_freq/representative_snapshots_{}qm'.format(n_qm)
    with open(binned_snapnums_filename) as fh:
        binned_snapnums_file_split = fh.read().split('#')[1:]
    for binstring in binned_snapnums_file_split:
        chomp = binstring.split()
        binnum = int(chomp[1])
        snapnums = chomp[2:]
        binned_snapnums[n_qm][binnum] = snapnums

dir_root = '/home/eric/Chemistry/calc.sgr/droplets'
for input_type in ('eda', 'eda_covp', 'freq_noCT'):
    dir_input = 'inputs_{}'.format(input_type)
    dir_input_full = os.path.join(dir_root, dir_input)
    assert os.path.exists(dir_input_full)
    for n_qm in n_qm_vals:
        dir_input_n_qm = '{}_{}qm'.format(dir_input, n_qm)
        dir_input_n_qm_full = os.path.join(dir_root, dir_input, dir_input_n_qm)
        assert os.path.exists(dir_input_n_qm_full)
        dir_random = '{}_random'.format(dir_input_n_qm)
        dir_random_full = os.path.join(dir_input_full, dir_random)
        for binnum in binned_snapnums[n_qm]:
            dir_bin = 'bin_{}'.format(binnum)
            dir_bin_full = os.path.join(dir_random_full, dir_bin)
            print(dir_bin_full)
            os.makedirs(dir_bin_full, exist_ok=False)
            for snapnum in binned_snapnums[n_qm][binnum]:
                inputfilename_glob = 'drop_{snapnum}_{n_qm}qm_*mm_{input_type}.in'.format(**locals())
                inputfilename_glob_full = os.path.join(dir_input_n_qm_full, inputfilename_glob)
                results = glob.glob(inputfilename_glob_full)
                for result in results:
                    shutil.copy2(result, dir_bin_full)
