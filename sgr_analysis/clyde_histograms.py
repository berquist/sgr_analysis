#!/usr/bin/env python3

import numpy as np
import scipy.stats as sps

from analysis_histograms import \
    (make_bin_edges, norm_pdf)

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def clyde_read_file(filename):
    with open(filename) as fh:
        xl = [float(x) for x in fh.readlines()]
    return xl


def analysis_single(frequencies, numbins, filename):

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

    print(' bins: {}'.format(bins))

    fig, ax = plt.subplots()

    ax.hist(frequencies, bins=bins, normed=True)

    linspace_bins = np.linspace(bins[0], bins[-1], 300)
    pdf_bins = norm_pdf(linspace_bins, mu=center, sigma=binwidth)
    ax.plot(linspace_bins,
            pdf_bins,
            color='red',
            linewidth=3,
            linestyle='--')

    fig.savefig(filename, bbox_inches='tight')

    return


if __name__ == '__main__':

    ct0 = clyde_read_file('ct0.csv')
    ct1 = clyde_read_file('ct1.csv')
    ct2 = clyde_read_file('ct2.csv')

    assert len(ct0) == len(ct1) == len(ct2)

    snapnums = list(range(1, len(ct0) + 1))

    analysis_single(ct0, 5, 'ct0.pdf')
    analysis_single(ct1, 5, 'ct1.pdf')
    analysis_single(ct2, 5, 'ct2.pdf')
