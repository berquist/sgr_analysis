from __future__ import division
from __future__ import print_function

import numpy as np
import pandas as pd

import matplotlib as mpl
mpl.use('Agg')
mpl.rc('text', usetex=True)
import matplotlib.pyplot as plt

from helpers import fit_line

# CP correction for geom/hessian
data = {
    'TFA':   {'no/no': 2429.31, 'yes/no': 2428.05, 'no/yes': 2457.238, 'yes/yes': 2440.879},
    'SCN_S': {'no/no': 2430.24, 'yes/no': 2427.91, 'no/yes': 2445.642, 'yes/yes': 2438.003},
    'DCA':   {'no/no': 2430.47, 'yes/no': 2426.60, 'no/yes': 2436.113, 'yes/yes': 2431.383},
    'SCN_N': {'no/no': 2431.64, 'yes/no': 2427.95, 'no/yes': 2432.691, 'yes/yes': 2437.876},
    'TfO':   {'no/no': 2431.91, 'yes/no': 2428.95, 'no/yes': 2441.921, 'yes/yes': 2443.665},
    'BF4':   {'no/no': 2434.69, 'yes/no': 2431.48, 'no/yes': 2454.767, 'yes/yes': 2450.890},
    'Tf2N':  {'no/no': 2435.80, 'yes/no': 2431.49, 'no/yes': 2451.851, 'yes/yes': 2441.309},
    'PF6':   {'no/no': 2437.74, 'yes/no': 2432.58, 'no/yes': 2479.806, 'yes/yes': 2469.473},
}

# geom/hessian
data_almo = {
    'BF4': {'scf/scf': 2434.7, 'almo/scf': 2438.56, 'almo/almo': 2441.62, 'scf/almo': 2437.69},
    'DCA': {'scf/scf': 2430.9, 'almo/scf': 2439.10, 'almo/almo': 2442.59, 'scf/almo': 2434.75},
    'PF6': {'scf/scf': 2437.5, 'almo/scf': 2438.47, 'almo/almo': 2441.03, 'scf/almo': 2440.03},
    'SCN': {'scf/scf': 2430.3, 'almo/scf': 2438.78, 'almo/almo': 2441.37, 'scf/almo': 2433.14},
    'TFA': {'scf/scf': 2429.8, 'almo/scf': 2438.39, 'almo/almo': 2441.06, 'scf/almo': 2432.93},
   'Tf2N': {'scf/scf': 2437.7, 'almo/scf': 2438.76, 'almo/almo': 2440.85, 'scf/almo': 2439.46},
    'TfO': {'scf/scf': 2433.9, 'almo/scf': 2439.37, 'almo/almo': 2442.54, 'scf/almo': 2436.75},
}


if __name__ == "__main__":
    df = pd.DataFrame(data).transpose()
    print(df)
    df_almo = pd.DataFrame(data_almo).transpose()
    print(df_almo)

    # compare to CSV for correctness in raw data
    # df_csv = pd.read_csv('cp_geometries.csv')

    formatter = mpl.ticker.ScalarFormatter(useOffset=False)
    # figsize = (5, 5)
    # fontsize_axlabel = 'medium'
    # fontsize_regression = 'medium'
    # fontsize_ticklabel = 'medium'
    figsize = (4.0, 4.0)
    fontsize_axlabel = 10
    fontsize_legend = 8
    fontsize_regression = 12
    fontsize_ticklabel = 12

    ### CP correction

    ## no/yes vs. no/no
    ## Excel: y = 3.2789x - 5526.6, R^2 = 0.44039
    ax_min, ax_max = 2425, 2485

    fig, ax = plt.subplots(figsize=figsize)
    x, y = df['no/no'], df['no/yes']
    ax.plot(x, y, marker='o', linestyle='', color='blue')
    slope, intercept, rsq = fit_line(x, y)
    print(abs(slope) - 3.2789, abs(intercept) - 5526.6, rsq - 0.44039)
    # c = pd.concat([x, y])
    c = x
    x_min, x_max = min(c), max(c)
    x_fit = np.linspace(x_min, x_max, 100)
    y_fit = (slope * x_fit) + intercept
    ax.plot(x_fit, y_fit, linestyle='-', color='blue')
    fit_string = r"""$y={:.4f}x{:+7.2f}$
$R^{{2}}={:7.5f}$""".format(slope, intercept, rsq)
    ax.text(0.10, 0.75, fit_string, transform=ax.transAxes, fontsize=fontsize_regression)
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    ax.tick_params(direction='out', labelsize=fontsize_ticklabel)
    ax.set_xlim((ax_min, ax_max))
    ax.set_ylim((ax_min, ax_max))
    ax.set_xlabel(r'uncorrected frequency @ uncorrected geometry (cm$^{-1}$)', fontsize=fontsize_axlabel)
    ax.set_ylabel(r'CP-corrected frequency @ uncorrected geometry (cm$^{-1}$)', fontsize=fontsize_axlabel)
    fig.savefig('cp_geometries_noyes_nono.pdf', bbox_inches='tight')
    plt.close('all')

    ## yes/yes vs. no/no
    ## Excel: y = 3.0672x - 5017.4, R^2 = 0.62892
    ax_min, ax_max = 2425, 2475

    fig, ax = plt.subplots(figsize=figsize)
    x, y = df['no/no'], df['yes/yes']
    ax.plot(x, y, marker='o', linestyle='', color='blue')
    slope, intercept, rsq = fit_line(x, y)
    print(abs(slope) - 3.0672, abs(intercept) - 5017.4, rsq - 0.62892)
    # c = pd.concat([x, y])
    c = x
    x_min, x_max = min(c), max(c)
    x_fit = np.linspace(x_min, x_max, 100)
    y_fit = (slope * x_fit) + intercept
    ax.plot(x_fit, y_fit, linestyle='-', color='blue')
    fit_string = r"""$y={:.4f}x{:+7.2f}$
$R^{{2}}={:7.5f}$""".format(slope, intercept, rsq)
    ax.text(0.10, 0.75, fit_string, transform=ax.transAxes, fontsize=fontsize_regression)
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    ax.tick_params(direction='out', labelsize=fontsize_ticklabel)
    ax.set_xlim((ax_min, ax_max))
    ax.set_ylim((ax_min, ax_max))
    ax.set_xlabel(r'uncorrected frequency @ uncorrected geometry (cm$^{-1}$)', fontsize=fontsize_axlabel)
    ax.set_ylabel(r'CP-corrected frequency @ CP-corrected geometry (cm$^{-1}$)', fontsize=fontsize_axlabel)
    fig.savefig('cp_geometries_yesyes_nono.pdf', bbox_inches='tight')
    plt.close('all')

    ## yes/no vs. no/no
    ## Excel: y = 0.6085x + 773.91, R^2 = 0.88757
    ax_min, ax_max = 2426, 2438

    fig, ax = plt.subplots(figsize=figsize)
    x, y = df['no/no'], df['yes/no']
    ax.plot(x, y, marker='o', linestyle='', color='blue')
    slope, intercept, rsq = fit_line(x, y)
    print(abs(slope) - 0.6085, abs(intercept) - 773.91, rsq - 0.88757)
    # c = pd.concat([x, y])
    c = x
    x_min, x_max = min(c), max(c)
    x_fit = np.linspace(x_min, x_max, 100)
    y_fit = (slope * x_fit) + intercept
    ax.plot(x_fit, y_fit, linestyle='-', color='blue')
    fit_string = r"""$y={:.4f}x{:+7.2f}$
$R^{{2}}={:7.5f}$""".format(slope, intercept, rsq)
    ax.text(0.10, 0.75, fit_string, transform=ax.transAxes, fontsize=fontsize_regression)
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    ax.tick_params(direction='out', labelsize=fontsize_ticklabel)
    ax.set_xlim((ax_min, ax_max))
    ax.set_ylim((ax_min, ax_max))
    ax.set_xlabel(r'uncorrected frequency @ uncorrected geometry (cm$^{-1}$)', fontsize=fontsize_axlabel)
    ax.set_ylabel(r'uncorrected frequency @ CP-corrected geometry (cm$^{-1}$)', fontsize=fontsize_axlabel)
    fig.savefig('cp_geometries_yesno_nono.pdf', bbox_inches='tight')
    plt.close('all')

    ## ALMO
    ax_min, ax_max = 2429, 2441

    fig, ax = plt.subplots(figsize=figsize)
    x, y = df_almo['scf/scf'], df_almo['scf/almo']
    ax.plot(x, y, marker='o', linestyle='', color='blue')
    slope, intercept, rsq = fit_line(x, y)
    print(slope, intercept, rsq)
    # c = pd.concat([x, y])
    c = x
    x_min, x_max = min(c), max(c)
    x_fit = np.linspace(x_min, x_max, 100)
    y_fit = (slope * x_fit) + intercept
    ax.plot(x_fit, y_fit, linestyle='-', color='blue')
    fit_string = r"""$y={:.4f}x{:+7.2f}$
$R^{{2}}={:7.5f}$""".format(slope, intercept, rsq)
    ax.text(0.10, 0.75, fit_string, transform=ax.transAxes, fontsize=fontsize_regression)
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    ax.tick_params(direction='out', labelsize=fontsize_ticklabel)
    ax.set_xlim((ax_min, ax_max))
    ax.set_ylim((ax_min, ax_max))
    ax.set_xlabel(r'SCF frequency @ SCF geometry (cm$^{-1}$)', fontsize=fontsize_axlabel)
    ax.set_ylabel(r'ALMO frequency @ SCF geometry (cm$^{-1}$)', fontsize=fontsize_axlabel)
    fig.savefig('compare_almo_scfalmo_scfscf.pdf', bbox_inches='tight')
    plt.close('all')

    fig, ax = plt.subplots(figsize=figsize)
    ax1_color, ax2_color = 'blue', 'red'
    x, y1, y2 = df_almo['scf/scf'], df_almo['almo/scf'], df_almo['almo/almo']
    ax.plot(x, y1, marker='o', linestyle='', color=ax1_color, label='SCF frequency @ ALMO geometry')
    ax.plot(x, y2, marker='s', linestyle='', color=ax2_color, label='ALMO frequency @ ALMO geometry')
    c = x
    x_min, x_max = min(c), max(c)
    slope1, intercept1, rsq1 = fit_line(x, y1)
    slope2, intercept2, rsq2 = fit_line(x, y2)
    print(slope1, intercept1, rsq1)
    print(slope2, intercept2, rsq2)
    x_fit = np.linspace(x_min, x_max, 100)
    y1_fit = (slope1 * x_fit) + intercept1
    y2_fit = (slope2 * x_fit) + intercept2
    ax.plot(x_fit, y1_fit, linestyle='-', color=ax1_color)
    ax.plot(x_fit, y2_fit, linestyle='--', color=ax2_color)
    fit1_string = r"""$y={:.4f}x{:+7.2f}$
$R^{{2}}={:7.5f}$""".format(slope1, intercept1, rsq1)
    fit2_string = r"""$y={:.4f}x{:+7.2f}$
$R^{{2}}={:7.5f}$""".format(slope2, intercept2, rsq2)
    ax.text(0.41, 0.50, fit1_string, transform=ax.transAxes, fontsize=fontsize_regression, color=ax1_color)
    ax.text(0.41, 0.89, fit2_string, transform=ax.transAxes, fontsize=fontsize_regression, color=ax2_color)
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    ax.tick_params(direction='out', labelsize=fontsize_ticklabel)
    ax_min, ax_max = 2428, 2444
    ax.set_xlim((ax_min, ax_max))
    ax.set_ylim((ax_min, ax_max))
    ax.set_xlabel(r'SCF frequency @ SCF geometry (cm$^{-1}$)', fontsize=fontsize_axlabel)
    ax.set_ylabel(r' frequency (cm$^{-1}$)', fontsize=fontsize_axlabel)
    ax.legend(loc='best', fancybox=True, framealpha=0.50, fontsize=fontsize_legend, numpoints=1)
    fig.savefig('compare_almo_pair.pdf', bbox_inches='tight')
    plt.close('all')
