# rows should be method chemistry
# cols should be [frz, pol, ct, disp, tot]

import itertools
from collections import OrderedDict

import pandas as pd

from sgr_analysis.sapt.helpers import df_weighted_average

method_map = OrderedDict([
    ('hf', 'HF'),
    ('b3lyp', 'B3LYP'),
    ('b3lyp-d2', 'B3LYP-D2'),
    ('b3lyp-d3', 'B3LYP-D3'),
    ('wb97x-d', r'\omegaB97X-D'),
    ('wb97m-v', r'\omegaB97M-V'),
])

basis_map = OrderedDict([
    ('6-31gss', '6-31G(d,p)'),
    ('cc-pvtz', 'cc-pVTZ'),
])

sapt_qc_basis_map = OrderedDict([
    ('6-31gss', '6-31G(d,p)'),
    ('jun-cc-pvdz', 'jun-cc-pVDZ'),
    ('jun-cc-pvtz', 'jun-cc-pVTZ'),
])

sapt_sapt_basis_map = OrderedDict([
    ('monomer', 'MCBS'),
    ('dimer', 'DCBS'),
])

method_basis_outer_product = list(itertools.product(method_map, basis_map))

if __name__ == '__main__':

    # Load the data.
    df_almo = pd.read_csv('data_almo.csv')
    df_sapt_pople = pd.read_csv('./631gdp/ct/data_sapt.csv')
    df_sapt_dunning_small = pd.read_csv('./jun-cc-pvdz/data_sapt.csv')
    df_sapt_dunning_large = pd.read_csv('./jun-cc-pvtz/ct/data_sapt.csv')

    # What's a more automatic way of doing this?
    # 15 MCBS calculations, 15 DCBS calculations, the rest are
    # calculated statistics.
    df_sapt_pople = df_sapt_pople[:30]
    df_sapt_dunning_small = df_sapt_dunning_small[:30]
    df_sapt_dunning_large = df_sapt_dunning_large[:30]

    # Filter out the extra snapshots from the ALMO results.
    assert df_sapt_dunning_small['snapnum'].equals(df_sapt_pople['snapnum'])
    assert df_sapt_dunning_large['snapnum'].equals(df_sapt_pople['snapnum'])
    snapnums = list(pd.Series(df_sapt_dunning_small['snapnum'], dtype=int))
    almo_mask = [almo_snapnum in snapnums for almo_snapnum in df_almo['snapnum']]
    df_almo = df_almo[almo_mask]

    almo_rename_columns = {
        'ho_scf': 'higher_order_scf',
        'tot_scf': 'total_scf',
    }
    df_almo.rename(columns=almo_rename_columns, inplace=True)

    df_almo.set_index('snapnum', inplace=True)
    df_sapt_pople['snapnum'] = pd.Series(df_sapt_pople['snapnum'], dtype=int)
    df_sapt_dunning_small['snapnum'] = pd.Series(df_sapt_dunning_small['snapnum'], dtype=int)
    df_sapt_dunning_large['snapnum'] = pd.Series(df_sapt_dunning_large['snapnum'], dtype=int)
    df_sapt_pople.set_index('snapnum', inplace=True)
    df_sapt_dunning_small.set_index('snapnum', inplace=True)
    df_sapt_dunning_large.set_index('snapnum', inplace=True)

    # TODO ??????
    weights = df_sapt_dunning_large['weight'][:15]

    # Drop columns that clutter the plots and don't contribute anything
    # meaningful.
    df_sapt_pople.drop(['filename', 'weight'], axis=1, inplace=True)
    df_sapt_dunning_small.drop(['filename', 'weight'], axis=1, inplace=True)
    df_sapt_dunning_large.drop(['filename', 'weight'], axis=1, inplace=True)

    df_sapt_pople['qc_basis'] = '6-31gss'
    df_sapt_dunning_small['qc_basis'] = 'jun-cc-pvdz'
    df_sapt_dunning_large['qc_basis'] = 'jun-cc-pvtz'
    df_sapt = pd.concat([df_sapt_pople,
                         df_sapt_dunning_small,
                         df_sapt_dunning_large], axis=0)

    qc_methods = sorted(set(df_almo['qc_method']))
    weighted_averages_almo = []
    for qc_method in method_basis_outer_product:
        qc_method_joined = '_'.join(qc_method)
        if qc_method_joined in qc_methods:
            method_label = '{}/{}'.format(method_map[qc_method[0]], basis_map[qc_method[1]])
            df_method = df_almo[df_almo['qc_method'] == qc_method_joined]
            weighted_averages_qc_method = df_weighted_average(df_method, weights)
            weighted_averages_almo.append((method_label, weighted_averages_qc_method))
    weighted_averages_almo = pd.DataFrame(OrderedDict(weighted_averages_almo)).transpose()
    print(weighted_averages_almo)

    weighted_averages_sapt = []
    for qc_basis in sapt_qc_basis_map:
        for sapt_basis in sapt_sapt_basis_map:
            sapt_label = 'SAPT0/{}/{}'.format(sapt_qc_basis_map[qc_basis], sapt_sapt_basis_map[sapt_basis])
            df_sapt_basis = df_sapt[(df_sapt['qc_basis'] == qc_basis) &
                                    (df_sapt['sapt_basis'] == sapt_basis)]
            weighted_averages_sapt_basis = df_weighted_average(df_sapt_basis, weights)
            weighted_averages_sapt.append((sapt_label, weighted_averages_sapt_basis))
    weighted_averages_sapt = pd.DataFrame(OrderedDict(weighted_averages_sapt)).transpose()
    print(weighted_averages_sapt)

    print(weighted_averages_almo.columns)
    print(weighted_averages_sapt.columns)

    # For B3LYP-D2 and B3LYP-D3, the (empirical) dispersion energy can
    # be calculated by $E_{disp} = E_{frz}^{\text{emp. disp.}} -
    # E_{frz}^{\text{none}}$.
    disp_d2 = weighted_averages_almo.ix['B3LYP-D2/6-31G(d,p)', 'frz'] - \
              weighted_averages_almo.ix['B3LYP/6-31G(d,p)', 'frz']
    disp_d3 = weighted_averages_almo.ix['B3LYP-D3/6-31G(d,p)', 'frz'] - \
              weighted_averages_almo.ix['B3LYP/6-31G(d,p)', 'frz']
    print('disp D2: {}'.format(disp_d2))
    print('disp D3: {}'.format(disp_d3))

    print('percent contributions to the total SAPT repulsive/exchange term')
    sapt_disp_tot = weighted_averages_sapt['exchange'] + \
                    weighted_averages_sapt['exchange_induction'] + \
                    weighted_averages_sapt['exchange_dispersion']
    print('exchange repulsion')
    print(100 * weighted_averages_sapt['exchange'] / sapt_disp_tot)
    print('exchange induction')
    print(100 * weighted_averages_sapt['exchange_induction'] / sapt_disp_tot)
    print('exchange dispersion')
    print(100 * weighted_averages_sapt['exchange_dispersion'] / sapt_disp_tot)

    weighted_averages_almo.to_csv('weighted_averages_almo.csv')
    weighted_averages_sapt.to_csv('weighted_averages_sapt.csv')
