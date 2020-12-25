"""stats.py: Print bin statistics.

This script should be run from within each 'inputs_freq_{n_qm}'
directory.
"""


import os

import numpy as np
import scipy.stats as sps

from sgr_analysis.analysis_utils import filter_outputfiles, get_CO2_frequencies, get_outputfiles_from_path


if __name__ == "__main__":

    outputfiles = get_outputfiles_from_path(os.getcwd())
    outputfiles_0mm = filter_outputfiles(outputfiles)

    # B3LYP/6-31G**
    weights_map = {
        1: 0.06060606,
        2: 0.24747475,
        3: 0.4010101,
        4: 0.22626263,
        5: 0.06464646,
    }
    print('sum of weights: {}'.format(sum(weights_map.values())))

    bins_outputfiles = dict()
    bins_f = dict()
    bins_i = dict()
    bins_sn = dict()
    for bk in weights_map:
        bins_outputfiles[bk] = [x for x in outputfiles_0mm
                                if '/bin_{}/'.format(bk) in x]
        bins_f[bk], bins_i[bk], bins_sn[bk] = get_CO2_frequencies(bins_outputfiles[bk])

    means = []

    for bk in weights_map:
        f_0mm = bins_f[bk]

        count = len(f_0mm)
        mean = np.mean(f_0mm)
        means.append(mean)
        median = np.median(f_0mm)
        mode = sps.mode(f_0mm)[0][0]
        mmin = min(f_0mm)
        mmax = max(f_0mm)
        rrange = mmax - mmin
        stdev_pop = np.std(f_0mm, ddof=1)
        stdev_sample = np.std(f_0mm, ddof=0)

        print('*** UNWEIGHTED: Bin {} ***'.format(bk))

        print(' count         : {:d}'.format(count))
        print(' mean          : {:.4f}'.format(mean))
        print(' median        : {:.4f}'.format(median))
        print(' mode          : {:.4f}'.format(mode))
        print(' min           : {:.4f}'.format(mmin))
        print(' max           : {:.4f}'.format(mmax))
        print(' range         : {:.4f}'.format(rrange))
        print(' stdev (pop)   : {:.4f}'.format(stdev_pop))
        print(' stdev (sample): {:.4f}'.format(stdev_sample))

    f_0mm = []
    for bk in weights_map:
        f_0mm.extend(bins_f[bk])

    count = len(f_0mm)
    mean = np.mean(f_0mm)
    median = np.median(f_0mm)
    mode = sps.mode(f_0mm)[0][0]
    mmin = min(f_0mm)
    mmax = max(f_0mm)
    rrange = mmax - mmin
    stdev_pop = np.std(f_0mm, ddof=1)
    stdev_sample = np.std(f_0mm, ddof=0)

    print('*** UNWEIGHTED: All bins ***')

    print(' count         : {:d}'.format(count))
    print(' mean          : {:.4f}'.format(mean))
    print(' median        : {:.4f}'.format(median))
    print(' mode          : {:.4f}'.format(mode))
    print(' min           : {:.4f}'.format(mmin))
    print(' max           : {:.4f}'.format(mmax))
    print(' range         : {:.4f}'.format(rrange))
    print(' stdev (pop)   : {:.4f}'.format(stdev_pop))
    print(' stdev (sample): {:.4f}'.format(stdev_sample))

    mean_unweighted = mean

    f_0mm = []
    weights_arr = []
    for bk, bv in weights_map.items():
        f_0mm.extend(bins_f[bk])
        weights_arr.extend([bv for _ in range(len(bins_f[bk]))])

    count = len(f_0mm)
    mean, sum_of_weights = np.average(f_0mm, weights=weights_arr, returned=True)
    median = np.median(f_0mm)
    mode = sps.mode(f_0mm)[0][0]
    mmin = min(f_0mm)
    mmax = max(f_0mm)
    rrange = mmax - mmin
    stdev_pop = np.std(f_0mm, ddof=1)
    stdev_sample = np.std(f_0mm, ddof=0)

    print('*** WEIGHTED ***')

    print(' count         : {:d}'.format(count))
    print(' mean (uw)     : {:.4f}'.format(mean_unweighted))
    print(' mean (w)      : {:.4f}'.format(mean))
    print(' diff (uw - w) : {:f}'.format(mean_unweighted - mean))
    print(' sum of weights: {:f}'.format(sum_of_weights))
    # print(' median        : {:.4f}'.format(median))
    # print(' mode          : {:.4f}'.format(mode))
    # print(' min           : {:.4f}'.format(mmin))
    # print(' max           : {:.4f}'.format(mmax))
    # print(' range         : {:.4f}'.format(rrange))
    # print(' stdev (pop)   : {:.4f}'.format(stdev_pop))
    # print(' stdev (sample): {:.4f}'.format(stdev_sample))

    x = range(1, 6)
    weights = np.array([weights_map[k] for k in x])
    weighted_means = [v * weights_map[k] for (k, v) in zip(x, means)]

    assert np.array(weighted_means).all() == (np.array(means) * weights).all()

    a1 = np.average(means, weights=weights)
    a3 = sum(weighted_means) / sum(weights)
    print(mean)
    print(a1)
    print(a3)
    # assert mean == a1 == a3

    ##########

    # import matplotlib as mpl
    # mpl.use('Agg')
    # import matplotlib.pyplot as plt

    # fig, ax = plt.subplots()

    # ax.plot(x, means, marker='o', label='unweighted')
    # ax.plot(x, weighted_means, marker='o', label='weighted')

    # fig.savefig('means.pdf', bbox_inches='tight')
