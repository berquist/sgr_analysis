import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sgr_analysis.sapt.helpers import fit_line



if __name__ == '__main__':

    sns.set(style='white')

    # def custom_alter_pairgrid(pg):

    # Load the data.
    df_almo = pd.read_csv('data_almo.csv')
    df_sapt_dunning = pd.read_csv('./jun-cc-pvtz/ct/data_sapt.csv')
    df_sapt_pople = pd.read_csv('./631gdp/ct/data_sapt.csv')

    # What's a more automatic way of doing this?
    # 15 MCBS calculations, 15 DCBS calculations, the rest are
    # calculated statistics.
    df_sapt_dunning = df_sapt_dunning[:30]
    df_sapt_pople = df_sapt_pople[:30]

    # Filter out the extra snapshots from the ALMO results.
    assert df_sapt_dunning['snapnum'].equals(df_sapt_pople['snapnum'])
    snapnums = list(pd.Series(df_sapt_dunning['snapnum'], dtype=int))
    almo_mask = [almo_snapnum in snapnums for almo_snapnum in df_almo['snapnum']]
    df_almo = df_almo[almo_mask]

    df_almo.set_index('snapnum', inplace=True)
    df_sapt_dunning['snapnum'] = pd.Series(df_sapt_dunning['snapnum'], dtype=int)
    df_sapt_pople['snapnum'] = pd.Series(df_sapt_pople['snapnum'], dtype=int)
    df_sapt_dunning.set_index('snapnum', inplace=True)
    df_sapt_pople.set_index('snapnum', inplace=True)

    # Drop columns that clutter the plots and don't contribute anything
    # meaningful.
    # Keep snapnum?
    # df_almo.drop(['snapnum'], axis=1, inplace=True)
    # df_sapt_dunning.drop(['filename', 'snapnum', 'weight'], axis=1, inplace=True)
    # df_sapt_pople.drop(['filename', 'snapnum', 'weight'], axis=1, inplace=True)
    df_sapt_dunning.drop(['filename', 'weight'], axis=1, inplace=True)
    df_sapt_pople.drop(['filename', 'weight'], axis=1, inplace=True)

    ## Basic comparisons within each method.

    # print('g_almo')
    # g_almo = sns.PairGrid(df_almo, hue='qc_method')
    # g_almo.map_offdiag(plt.scatter)
    # g_almo.add_legend()
    # g_almo.savefig('trellis_almo.pdf')

    # print('g_sapt_dunning')
    # g_sapt_dunning = sns.PairGrid(df_sapt_dunning, hue='sapt_basis')
    # g_sapt_dunning.map_offdiag(plt.scatter)
    # g_sapt_dunning.add_legend()
    # g_sapt_dunning.savefig('trellis_sapt_dunning.pdf')

    # print('g_sapt_pople')
    # g_sapt_pople = sns.PairGrid(df_sapt_pople, hue='sapt_basis')
    # g_sapt_pople.map_offdiag(plt.scatter)
    # g_sapt_pople.add_legend()
    # g_sapt_pople.savefig('trellis_sapt_pople.pdf')

    ## Compare ALMO HF, ALMO B3LYP and SAPT dimer-centered basis with
    ## the 6-31G(d,p) basis set.

    # For ALMO, split apart B3LYP/HF and remove the column.
    df_almo_b3lyp = df_almo[df_almo['qc_method'] == 'b3lyp_6-31gss']
    df_almo_hf = df_almo[df_almo['qc_method'] == 'hf_6-31gss']
    df_almo_b3lyp.drop('qc_method', axis=1, inplace=True)
    df_almo_hf.drop('qc_method', axis=1, inplace=True)

    df_almo_b3lyp.columns = [('almo_b3lyp_' + c) for c in df_almo_b3lyp.columns]
    df_almo_hf.columns = [('almo_hf_' + c) for c in df_almo_hf.columns]

    # For SAPT, Keep the dimer-centered basis and drop any other reference
    # to DCBS/MCBS.
    df_sapt_dunning = df_sapt_dunning[df_sapt_dunning['sapt_basis'] == 'dimer']
    df_sapt_pople = df_sapt_pople[df_sapt_pople['sapt_basis'] == 'dimer']
    df_sapt_dunning.drop('sapt_basis', axis=1, inplace=True)
    df_sapt_pople.drop('sapt_basis', axis=1, inplace=True)

    df_sapt_dunning.columns = [('sapt_dunning_' + c) for c in df_sapt_dunning.columns]
    df_sapt_pople.columns = [('sapt_pople_' + c) for c in df_sapt_pople.columns]

    df_all = pd.concat([df_almo_b3lyp, df_almo_hf, df_sapt_pople], axis=1)

    # ONLY FOR TESTING, speeds things up...
    # df_all = df_all.ix[:, :6]

    plt.close('all')

    print('g_all')
    g_all = sns.PairGrid(df_all, size=(2.5 * 1.0))
    g_all.map_offdiag(plt.scatter)
    # http://stackoverflow.com/questions/34087126/plot-lower-triangle-in-a-seaborn-pairgrid
    for i, j in zip(*np.triu_indices_from(g_all.axes, 1)):
        g_ij = g_all.axes[i, j]
        g_ij.set_visible(False)
    # Disable the diagonal axes as well.
    for i in range(g_all.axes.shape[0]):
        g_all.axes[i, i].set_visible(False)
    for i, j in zip(*np.tril_indices_from(g_all.axes, -1)):
        g_ij = g_all.axes[i, j]
        ci, cj = df_all.columns[i], df_all.columns[j]
        di, dj = df_all[ci], df_all[cj]
        slope, intercept, rsq = fit_line(di, dj)
        x_min, x_max = min(di), max(di)
        x_fit = np.linspace(x_min, x_max, 100)
        g_ij.plot((slope*x_fit) + intercept, x_fit, linestyle='-', color='red')
        # Choking on LaTeX, I see...
        # rsq_str = r'$R^2 = {:.4f}$'.format(rsq)
        rsq_str = '{:.4f}'.format(rsq)
        g_ij.text(0.10, 0.75, rsq_str, transform=g_ij.transAxes, fontsize='large')
        # WHY ARE THINGS INVERTED
        xlabel = g_all.x_vars[j]
        ylabel = g_all.y_vars[i]
        g_ij.set_xlabel(xlabel)
        g_ij.set_ylabel(ylabel)
        # g_ij.xaxis.offsetText.set_visible(True)
        # g_ij.yaxis.offsetText.set_visible(True)
        for label in g_ij.get_xticklabels():
            label.set_visible(True)
        for label in g_ij.get_yticklabels():
            label.set_visible(True)
    # defaults are 0.2
    g_all.fig.subplots_adjust(hspace=0.3, wspace=0.3)
    g_all.savefig('trellis_all.pdf')

    plt.close('all')

    # fig, ax = plt.subplots()
    # i, j = 2, 0
    # ci, cj = df_all.columns[i], df_all.columns[j]
    # print(i, ci, j, cj)
    # di, dj = df_all[ci], df_all[cj]
    # ax.scatter(di, dj, marker='o', color='blue')
    # slope, intercept, r, p, stderr = sps.linregress(di, dj)
    # rsq = r ** 2
    # x_min, x_max = min(di), max(di)
    # x_fit = np.linspace(x_min, x_max, 100)
    # ax.plot(x_fit, (slope*x_fit) + intercept, linestyle='-', color='red')
    # print(slope, intercept)
    # ax.set_xlabel(di.name)
    # ax.set_ylabel(dj.name)
    # rsq_str = r'$R^2 = {:.4f}$'.format(rsq)
    # ax.text(0.10, 0.75, rsq_str, transform=ax.transAxes, fontsize='large')
    # fig.savefig('trellis_all_0_2.pdf', bbox_inches='tight')


    # linear_regression_2df(df_all, df_all)
