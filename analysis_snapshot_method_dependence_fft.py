#!/usr/bin/env python3

"""analysis_snapshot_method_dependence_fft.py: Plot the actual
vibrational spectra of parsed MD snapshots separated by 4 fs (to give
a spectral width of ~8333 cm^-1) by taking the FFT of the dipole
autocorrelation function (ACF).

An important note is that although some dictionary lookups say '256',
implying the presence of 256 IL pairs as point charges, there are far
fewer than that, and the amount varies. It is simply for convenience.

"""

from __future__ import print_function
from __future__ import division

import sys
sys.path.append('/home/eric/Chemistry/aimd')
sys.path.append('/home/eric/Chemistry/aimd/spectra')
from qchem_aimd_fft import \
    (make_vibspectrum_from_dipoles, AU_TIME_IN_SEC)

import pickle

import numpy as np
import scipy.stats as sps

# from parse_outputs_snapshot_method_dependence_4fs import \
#     (methods, basis_sets)

from analysis_histograms import \
    (make_bin_edges, norm_pdf)

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt


if __name__ == '__main__':

    with open('frequencies.pypickle', 'rb') as picklefile:
        frequencies = pickle.load(picklefile)
    with open('dipoles.pypickle', 'rb') as picklefile:
        dipoles = pickle.load(picklefile)
    with open('snapnums_frequencies.pypickle', 'rb') as picklefile:
        snapnums_f = pickle.load(picklefile)
    with open('snapnums_dipoles.pypickle', 'rb') as picklefile:
        snapnums_d = pickle.load(picklefile)

    # time step is 50 ps
    # time_step = 50.0 * 1.0e-12 / AU_TIME_IN_SEC
    # time step is 4 fs
    time_step = 4.0 * 1.0e-15 / AU_TIME_IN_SEC
    print(time_step)

    # for method in methods:
    #     for basis_set in basis_sets:
    #         dipoles_method_basis_0_0 = np.array(dipoles[method][basis_set][0][0])
    #         frequencies, intensities = make_vibspectrum_from_dipoles(dipoles_method_basis_0_0, time_step)
    #         assert frequencies.shape == intensities.shape
    #         # print(frequencies.shape)
    #         # print(frequencies)

    #         fig, ax = plt.subplots()

    #         ax.plot(frequencies, intensities, linewidth=0.50)

    #         fig.savefig('fft_{}_{}.pdf'.format(method, basis_set), bbox_inches='tight')

    #         plt.close('all')

    # these are the raw frequencies; I want to plot the
    # distribution/histogram
    frequencies_flex_b3lyp_lp_0_0 = np.array(frequencies['flex']['b3lyp']['lp'][0][0])
    frequencies_flex_b3lyp_lp_0_all = np.array(frequencies['flex']['b3lyp']['lp'][0][256])
    frequencies_rigid_b3lyp_lp_0_0 = np.array(frequencies['rigid']['b3lyp']['lp'][0][0])
    frequencies_rigid_b3lyp_lp_0_all = np.array(frequencies['rigid']['b3lyp']['lp'][0][256])
    dipoles_flex_b3lyp_lp_0_0 = np.array(dipoles['flex']['b3lyp']['lp'][0][0])
    dipoles_flex_b3lyp_lp_0_all = np.array(dipoles['flex']['b3lyp']['lp'][0][256])
    # dipoles_rigid_b3lyp_lp_0_0 = np.array(dipoles['rigid']['b3lyp']['lp'][0][0])
    dipoles_rigid_b3lyp_lp_0_all = np.array(dipoles['rigid']['b3lyp']['lp'][0][256])
    print(dipoles_flex_b3lyp_lp_0_0.shape, dipoles_flex_b3lyp_lp_0_all.shape, dipoles_rigid_b3lyp_lp_0_all.shape)
    frequencies_flex_0_0, intensities_flex_0_0 = make_vibspectrum_from_dipoles(dipoles_flex_b3lyp_lp_0_0, time_step)
    frequencies_flex_0_all, intensities_flex_0_all = make_vibspectrum_from_dipoles(dipoles_flex_b3lyp_lp_0_all, time_step)
    # frequencies_rigid_0_0, intensities_rigid_0_0 = make_vibspectrum_from_dipoles(dipoles_rigid_b3lyp_lp_0_0, time_step)
    frequencies_rigid_0_all, intensities_rigid_0_all = make_vibspectrum_from_dipoles(dipoles_rigid_b3lyp_lp_0_all, time_step)
    assert frequencies_flex_0_0.shape == intensities_flex_0_0.shape
    assert frequencies_flex_0_all.shape == intensities_flex_0_all.shape
    # assert frequencies_rigid_0_0.shape == intensities_rigid_0_0.shape
    assert frequencies_rigid_0_all.shape == intensities_rigid_0_all.shape

    with open('data_flex_0_0.dat', 'w') as fh:
        for (x, y) in zip(frequencies_flex_0_0, intensities_flex_0_0):
            fh.write('{} {}\n'.format(x, y))
    with open('data_flex_0_all.dat', 'w') as fh:
        for (x, y) in zip(frequencies_flex_0_all, intensities_flex_0_all):
            fh.write('{} {}\n'.format(x, y))
    with open('data_rigid_0_all.dat', 'w') as fh:
        for (x, y) in zip(frequencies_rigid_0_all, intensities_rigid_0_all):
            fh.write('{} {}\n'.format(x, y))

    ## Plot the IR spectrum for the flexible and rigid snapshots as
    ## the FFT of the dipole autocorrelation function.

    wavenumber_cutoff_lo = 2420.0
    wavenumber_cutoff_hi = 2480.0
    # wavenumber_cutoff_lo = 0.0
    # wavenumber_cutoff_hi = 5000.0

    mask_lo = frequencies_flex_0_0 >= wavenumber_cutoff_lo
    mask_hi = frequencies_flex_0_0 <= wavenumber_cutoff_hi
    mask = np.logical_and(mask_lo, mask_hi)
    # What does the ix matrix look like?
    ix = np.where(mask)[0]

    fig, ax = plt.subplots()

    ax.plot(frequencies_flex_0_0[ix], intensities_flex_0_0[ix], label='flex (0, 0)', linewidth=0.50)
    ax.plot(frequencies_flex_0_all[ix], intensities_flex_0_all[ix], label='flex (0, all)', linewidth=0.50)
    # ax.plot(frequencies_rigid_0_0[ix], intensities_rigid_0_0[ix], label='rigid (0, 0)', linewidth=0.50)
    ax.plot(frequencies_rigid_0_all[ix], intensities_rigid_0_all[ix], label='rigid (0, all)', linewidth=0.50)

    ax.set_xlabel('frequency (cm$^{-1}$)')
    ax.set_ylabel('intensity (a.u.)')

    ax.legend(loc='best', fancybox=True, framealpha=0.50)

    fig.savefig('fft_{}_{}.pdf'.format('b3lyp', 'lp'), bbox_inches='tight')

    plt.close('all')

    ## Plot the distribution of \nu_3 harmonic frequencies for the
    ## flexible and rigid snapshots.

    numbins = 10

    center_flex = np.mean(frequencies_flex_b3lyp_lp_0_0)
    center_rigid = np.mean(frequencies_rigid_b3lyp_lp_0_0)
    # standard deviation of the population, not sample!
    binwidth_flex = np.std(frequencies_flex_b3lyp_lp_0_0, ddof=1)
    binwidth_rigid = np.std(frequencies_rigid_b3lyp_lp_0_0, ddof=1)
    bins_flex = make_bin_edges(numbins, center_flex, binwidth_flex)
    bins_rigid = make_bin_edges(numbins, center_rigid, binwidth_rigid)
    hist_flex, bin_edges_flex = np.histogram(frequencies_flex_b3lyp_lp_0_0, bins_flex)
    hist_rigid, bin_edges_rigid = np.histogram(frequencies_rigid_b3lyp_lp_0_0, bins_rigid)
    linspace_bins_flex = np.linspace(bins_flex[0], bins_flex[-1], 300)
    linspace_bins_rigid = np.linspace(bins_rigid[0], bins_rigid[-1], 300)
    pdf_bins_flex = norm_pdf(linspace_bins_flex, mu=center_flex, sigma=binwidth_flex)
    pdf_bins_rigid = norm_pdf(linspace_bins_rigid, mu=center_rigid, sigma=binwidth_rigid)

    mmin = min(frequencies_flex_b3lyp_lp_0_0)
    mmax = max(frequencies_flex_b3lyp_lp_0_0)
    rng = mmax - mmin
    mean = np.mean(frequencies_flex_b3lyp_lp_0_0)
    median = np.median(frequencies_flex_b3lyp_lp_0_0)
    mode = sps.mode(frequencies_flex_b3lyp_lp_0_0)
    print(' min           : {:.2f}'.format(mmin))
    print(' max           : {:.2f}'.format(mmax))
    print(' range         : {:.2f}'.format(rng))
    print(' mean          : {:.3f}'.format(mean))
    print(' median        : {:.3f}'.format(median))
    print(' mode          : {}'.format(mode))
    print(' stdev (pop)   : {:.3f}'.format(np.std(frequencies_flex_b3lyp_lp_0_0, ddof=1)))
    print(' stdev (sample): {:.3f}'.format(np.std(frequencies_flex_b3lyp_lp_0_0, ddof=0)))


    fig, ax = plt.subplots()

    # ax.hist(frequencies_flex_b3lyp_lp_0_0, bins=bins_flex, normed=True, label='flex (0, 0) hist')
    # ax.hist(frequencies_flex_b3lyp_lp_0_all, bins=bins_flex, normed=True, label='flex (0, all) hist')
    # # ax.hist(frequencies_rigid_b3lyp_lp_0_0, bins=bins_rigid, normed=True, label='rigid (0, 0) hist')
    # ax.hist(frequencies_rigid_b3lyp_lp_0_all, bins=bins_rigid, normed=True, label='rigid (0, all) hist')
    # ax.plot(linspace_bins_flex, pdf_bins_flex, label='flex (0, 0) PDF', linewidth=0.50, linestyle='--')
    # ax.plot(linspace_bins_flex_0_all, pdf_bins_flex_0_all, label='flex (0, all) PDF', linewidth=0.50, linestyle='--')
    # # ax.plot(linspace_bins_rigid_0_0, pdf_bins_rigid_0_0, label='rigid (0, 0) PDF', linewidth=0.50, linestyle='--')
    # ax.plot(linspace_bins_rigid_0_all, pdf_bins_rigid_0_all, label='rigid (0, all) PDF', linewidth=0.50, linestyle='--')

    ax.set_xlabel('harmonic frequency (cm$^{-1}$)')
    # ax.set_ylabel('')

    ax.legend(loc='best', fancybox=True, framealpha=0.50)

    fig.savefig('harmonic_distribution_{}_{}.pdf'.format('b3lyp', 'lp'), bbox_inches='tight')

    plt.close('all')

    fig, ax = plt.subplots()

    # for this to work, all the data should be normalized to something (?)
    ax.hist(frequencies_flex_b3lyp_lp_0_0, bins=bins_flex, normed=True, label='flex (0, 0) hist')
    ax.plot(linspace_bins_flex, pdf_bins_flex, label='flex PDF', linewidth=0.50, linestyle='--')
    ax.plot(frequencies_flex_0_0[ix], intensities_flex_0_0[ix], label='flex (0, 0) dip-dip', linewidth=0.50)

    ax.legend(loc='best', fancybox=True, framealpha=0.50)

    fig.savefig('combined_{}_{}.pdf'.format('b3lyp', 'lp'), bbox_inches='tight')

    plt.close('all')
