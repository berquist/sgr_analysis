import pickle

from collections import OrderedDict

import numpy as np
import numpy.linalg as npl

from sgr_analysis.analysis_utils import read_snapshot_file

from sgr_analysis.parse_outputs_snapshot_method_dependence import methods
from sgr_analysis.parse_outputs_snapshot_method_dependence import basis_sets

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm


colors = OrderedDict([
    ('blyp', 'blue'),
    ('tpss', 'cyan'),
    ('b3lyp', 'red'),
    ('wb97x-d', 'orange'),
    ('hf', 'black'),
    ('ri-mp2', 'grey'),
])


linestyles = OrderedDict([
    ('6-31gdp', '-'),
    ('cc-pvtz', '--'),
])


def get_single_snapshot_results(snapnum, snapnums_dict, results_dict):
    """For a given snapshot number, get the single 'result' (frequency,
    intensity, dipole, etc.) corresponding to it, maintaining the original
    dictionary structure.
    """

    single_dict = dict()

    assert snapnums_dict.keys() == results_dict.keys()
    for method in methods:
        assert snapnums_dict[method].keys() == results_dict[method].keys()
        single_dict[method] = dict()
        for basis_set in basis_sets:
            assert snapnums_dict[method][basis_set].keys() == results_dict[method][basis_set].keys()
            single_dict[method][basis_set] = dict()
            for n_qm in snapnums_dict[method][basis_set]:
                assert snapnums_dict[method][basis_set][n_qm].keys() == results_dict[method][basis_set][n_qm].keys()
                single_dict[method][basis_set][n_qm] = dict()
                for n_mm in snapnums_dict[method][basis_set][n_qm]:
                    if len(snapnums_dict[method][basis_set][n_qm][n_mm]) > 0:
                        idx = snapnums_dict[method][basis_set][n_qm][n_mm].index(snapnum)
                        single_result = results_dict[method][basis_set][n_qm][n_mm][idx]
                        single_dict[method][basis_set][n_qm][n_mm] = single_result

    return single_dict


def plot_single_snapshot_dipoles(snapnum, snapnums_d, dipoles, inp_fig=None, inp_ax=None):
    """Plot the dipole moment of a single snapshot as a function of its
    method and basis set.
    """

    dipoles_snap = get_single_snapshot_results(snapnum, snapnums_d, dipoles)

    fig, ax = plt.subplots()
    if inp_fig:
        fig = inp_fig
    if inp_ax:
        ax = inp_ax

    for basis_set in basis_sets:
        plot_list = [npl.norm(dipoles_snap[method][basis_set][0][0])
                     for method in methods]
        ax.plot(ticks_methods,
                plot_list,
                label=basis_sets[basis_set],
                marker='o')

    if not inp_ax:
        ax.tick_params(direction='out', top='off', right='off')
        ax.set_xticklabels(xticklabels_methods)
        ax.set_xlabel('method', fontsize='large')
        ax.set_ylabel("total dipole moment (Debye)", fontsize='large')
        ax.set_title("snapshot {}".format(snapnum), fontsize='large')
        ax.legend(loc='best', fancybox=True, framealpha=0.50)
    if not inp_fig:
        filename = 'dipole_snap{}.pdf'.format(snapnum)
        print('Saving {}'.format(filename))
        fig.savefig(filename, bbox_inches='tight')

        plt.close(fig)

    return


def plot_single_snapshot_frequencies(snapnum, snapnums_f, frequencies, inp_fig=None, inp_ax=None):
    """Plot the CO2 v3 frequency of a single snapshot as a function of its
    method and basis set.
    """

    frequencies_snap = get_single_snapshot_results(snapnum, snapnums_f, frequencies)

    fig, ax = plt.subplots()
    if inp_fig:
        fig = inp_fig
    if inp_ax:
        ax = inp_ax

    for basis_set in basis_sets:
        plot_list = [frequencies_snap[method][basis_set][0][0]
                     for method in methods]
        ax.plot(ticks_methods,
                plot_list,
                label=basis_sets[basis_set],
                marker='o')

    if not inp_ax:
        ax.tick_params(direction='out', top='off', right='off')
        ax.set_xticklabels(xticklabels_methods)
        ax.set_xlabel('method', fontsize='large')
        ax.set_ylabel(r'$\nu_{3}$ frequency (cm$^{-1}$)', fontsize='large')
        ax.set_title("snapshot {}".format(snapnum), fontsize='large')
        ax.legend(loc='best', fancybox=True, framealpha=0.50)
    if not inp_fig:
        filename = 'frequency_snap{}.pdf'.format(snapnum)
        print('Saving {}'.format(filename))
        fig.savefig(filename, bbox_inches='tight')

        plt.close(fig)

    return


if __name__ == '__main__':

    # Read in the pickle files that contain all the raw data.
    with open('frequencies.pypickle', 'rb') as picklefile:
        frequencies = pickle.load(picklefile)
    with open('dipoles.pypickle', 'rb') as picklefile:
        dipoles = pickle.load(picklefile)
    with open('snapnums_frequencies.pypickle', 'rb') as picklefile:
        snapnums_f = pickle.load(picklefile)
    with open('snapnums_dipoles.pypickle', 'rb') as picklefile:
        snapnums_d = pickle.load(picklefile)

    ticks_methods = range(len(methods))
    xticklabels_methods = list(methods[method] for method in methods)

    means_frequencies = [np.mean(frequencies[method]['6-31gdp'][0][0])
                         for method in methods]
    means_dipole_moments = [np.mean([npl.norm(dipole)
                                     for dipole in dipoles[method]['6-31gdp'][0][0]])
                            for method in methods]

    ##################################################

    # Plot the mean frequencies for each method using the 6-31G(d,p)
    # basis set over the 0,0 snapshots.

    fig, ax = plt.subplots()

    ax.plot(ticks_methods,
            means_frequencies,
            marker='s',
            color='green')

    ax.set_xticks(ticks_methods)
    ax.set_xticklabels(xticklabels_methods)
    ax.set_xlabel('method', fontsize='large')
    ax.set_ylabel(r'$\nu_{3}$ frequency (cm$^{-1}$)', fontsize='large')

    fig.savefig('frequencies_mean_6-31gdp.pdf', bbox_inches='tight')

    ##################################################

    # Plot the mean dipole moments for each method using the 6-31G(d,p)
    # basis set over the 0,0 snapshots.

    fig, ax = plt.subplots()

    ax.plot(ticks_methods,
            means_dipole_moments,
            marker='s',
            color='green')

    ax.set_xticks(ticks_methods)
    ax.set_xticklabels(xticklabels_methods)
    ax.set_xlabel('method', fontsize='large')
    ax.set_ylabel('total dipole moment (Debye)', fontsize='large')

    fig.savefig('dipole_moments_mean_6-31gdp.pdf', bbox_inches='tight')

    ##################################################

    # Plot the mean frequencies (increasing scale) and mean dipole
    # moments (decreasing scale) together for each method using the
    # 6-31G(d,p) basis set over the 0,0 snapshots.

    ax1_color = 'blue'
    ax2_color = 'green'

    fig, ax1 = plt.subplots()

    ax1.plot(ticks_methods,
             means_frequencies,
             marker='s',
             linestyle='-',
             color=ax1_color)

    ax1.set_ylabel(r'$\nu_{3}$ frequency (cm$^{-1}$)', fontsize='large')

    for tl in ax1.get_yticklabels():
        tl.set_color(ax1_color)

    ax2 = ax1.twinx()

    ax2.plot(ticks_methods,
             means_dipole_moments,
             marker='*',
             linestyle='--',
             color=ax2_color)

    ax2.set_ylabel('total dipole moment (Debye)', fontsize='large')

    ax2.set_ylim(ax2.get_ylim()[::-1])

    for tl in ax2.get_yticklabels():
        tl.set_color(ax2_color)

    ax1.set_xticks(ticks_methods)
    ax1.set_xticklabels(xticklabels_methods)
    ax1.set_xlabel('method', fontsize='large')

    ax1.tick_params(direction='out', top='off', right='off')
    ax2.tick_params(direction='out', top='off', left='off')

    filename = 'frequencies_dipole_moments_mean_6-31gdp.pdf'
    print('Saving {}'.format(filename))
    fig.savefig(filename, bbox_inches='tight')

    ##################################################

    snapnums = read_snapshot_file("/home/eric/Chemistry/calc.sgr/paper_02_CD_SC/inputs_freq/representative_snapshots_3qm")
    print('snapshot numbers:', snapnums)
    for snapnum in snapnums:
        plot_single_snapshot_frequencies(snapnum, snapnums_f, frequencies)
        plot_single_snapshot_dipoles(snapnum, snapnums_d, dipoles)

    ##################################################

    cutoff = 50
    nticks = 5
    step = cutoff // nticks
    ticks_traj = list(range(0, (cutoff + step), step))
    ticks_traj[0] = 1

    cmap = cm.get_cmap('viridis')

    for basis_set in basis_sets:

        fig, ax = plt.subplots()

        # one series for each method used
        num_plot_lines = len(methods)
        for idx, method in enumerate(methods):
            c = cmap(float(idx) / (num_plot_lines - 1))
            ax.plot(frequencies[method][basis_set][0][256][:cutoff],
                    label=methods[method],
                    linewidth=1,
                    color=c)

        # Add the experimental BMIM/PF6 frequency as a line.
        ax.plot(ax.get_xlim(),
                [2340.0, 2340.0],
                marker='',
                linestyle=':',
                color='black',
                label='experiment')

        ax.set_xticks(ticks_traj)
        ax.set_xlim((min(ticks_traj), max(ticks_traj)))
        ax.set_xlabel('snapshot #', fontsize='large')
        ax.set_ylabel(r'$\nu_{3}$ harmonic frequency (cm$^{-1}$)', fontsize='large')
        ax.tick_params(direction='out')
        ax.legend(fancybox=True,
                  loc='upper right',
                  framealpha=0.50,
                  fontsize='small')

        filename = 'trajectory_frequencies_{}.pdf'.format(basis_set)
        print('Saving {}'.format(filename))
        fig.savefig(filename, bbox_inches='tight')
        plt.close('all')

    ##################################################
