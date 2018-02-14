"""helpers.py: Functions for parsing data from SAPT0 and ALMO output
files."""

from __future__ import division
from __future__ import print_function

import os.path
import re

import numpy as np
import pandas as pd
from pandas.util.testing import assert_frame_equal
import scipy.stats as sps

from scripts.utils import make_file_iterator


KJ_TO_KCAL = 4.184
rsq_cutoff = 0.50

re_number_str = "(-?[0-9]*\.[0-9]*)"
re_number = re.compile(re_number_str)


def fit_line(x, y):
    """Return slope, intercept of best fit line."""
    # Remove entries where either x or y is NaN.
    clean_data = pd.concat([x, y], 1).dropna(0) # row-wise
    (_, x), (_, y) = clean_data.iteritems()
    slope, intercept, r, p, stderr = sps.linregress(x, y)
    return slope, intercept, r**2


def assertFrameEqual(df1, df2, **kwds):
    """Assert that two dataframes are equal, ignoring ordering of
    columns.
    """
    return assert_frame_equal(df1.sort(axis=1), df2.sort(axis=1), check_names=True, **kwds)


def linear_regression_2df(df1, df2):
    """Perform a linear regression between all pairs of columns in two
    DataFrames.
    """
    header = '               slope            intercept    rsq'
    print(header)
    template = '{:20.10f} {:20.10f} {:.4f} ({}, {})'
    columns = df1.columns

    # Are the two DataFrames identical?
    try:
        assertFrameEqual(df1, df2)
    except AssertionError:
        samedf = False
    else:
        samedf = True

    rsq_list = []

    # If the two DataFrames are the same, avoid double-counting and
    # "self-regression".
    if samedf:
        for xi in range(len(columns)):
            for yi in range(xi):
                xlabel = columns[xi]
                ylabel = columns[yi]
                if xlabel != ylabel:
                    slope, intercept, rsq = fit_line(df1.loc[:, xlabel],
                                                     df2.loc[:, ylabel])
                    rsq_list.append((rsq, (xlabel, ylabel)))
                    if rsq >= rsq_cutoff:
                        print(template.format(slope, intercept, rsq, xlabel, ylabel))

    # If the two DataFrames aren't the same, go over all columns in both.
    else:
        for xlabel in columns:
            for ylabel in columns:
                slope, intercept, rsq = fit_line(df1.loc[:, xlabel],
                                                 df2.loc[:, ylabel])
                rsq_list.append((rsq, (xlabel, ylabel)))
                if rsq >= rsq_cutoff:
                    print(template.format(slope, intercept, rsq, xlabel, ylabel))

    rsq_list = sorted(rsq_list)
    for entry in rsq_list:
        print(entry)

    return


def _add_axis_labels(pg):
    """Add labels to all possible left and bottom Axes.

    This doesn't work the way I expect it to?
    """
    # for ax, label in zip(self.axes[-1, :], self.x_vars):
    #     ax.set_xlabel(label)
    # for ax, label in zip(self.axes[:, 0], self.y_vars):
    #     ax.set_ylabel(label)
    for i, j in zip(*np.tril_indices_from(pg.axes, -1)):
        ax = pg.axes[i, j]
        # WHY ARE THINGS INVERTED
        xlabel = pg.x_vars[j]
        ylabel = pg.y_vars[i]
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    return


def read_qchem_eda(filename, is_cp=False):
    """Units of returned values: kcal/mol"""

    almo_data_snap = dict()

    with open(filename) as fh:
        for line in fh:
            if 'Frozen Density ( FRZ )' in line:
                frz = float(line.split()[-1])
                almo_data_snap['frz'] = frz
            if 'Polarization ( POL )' in line:
                pol = float(line.split()[-1])
                almo_data_snap['pol'] = pol
            if 'RS Delocalization ( P-DEL )' in line:
                rs_del = float(line.split()[-1])
                almo_data_snap['del_rs'] = rs_del
            if 'RS Basis set superposition error ( P-BSSE )' in line:
                assert is_cp
                rs_bsse = float(line.split()[-1])
                almo_data_snap['bsse_rs'] = rs_bsse
            if 'RS Charge-transfer ( P-CT = P-DEL + P-BSSE )' in line:
                assert is_cp
                rs_ct = float(line.split()[-1])
                assert abs(rs_del + rs_bsse - rs_ct) < 1.0e-4
                almo_data_snap['ct_rs'] = rs_ct
            if 'SCF Delocalization ( V-DEL )' in line:
                scf_del = float(line.split()[-1])
                almo_data_snap['del_scf'] = scf_del
            if 'SCF Basis set superposition error ( V-BSSE )' in line:
                assert is_cp
                scf_bsse = float(line.split()[-1])
                almo_data_snap['bsse_scf'] = scf_bsse
            if 'SCF Charge-transfer ( V-CT = V-DEL + V-BSSE )' in line:
                assert is_cp
                scf_ct = float(line.split()[-1])
                assert abs(scf_del + scf_bsse - scf_ct) < 1.0e-4
                almo_data_snap['ct_scf'] = scf_ct
            if 'RS Total ( P-TOT = FRZ + POL + P-DEL )' in line:
                assert not is_cp
                rs_tot = float(line.split()[-1])
                almo_data_snap['tot_rs'] = rs_tot
            if 'RS Total ( P-TOT = FRZ + POL + P-CT )' in line:
                assert is_cp
                rs_tot = float(line.split()[-1])
                almo_data_snap['tot_rs'] = rs_tot
            if 'SCF Total ( V-TOT = FRZ + POL + V-DEL )' in line:
                assert not is_cp
                scf_tot = float(line.split()[-1])
                almo_data_snap['tot_scf'] = scf_tot
            if 'SCF Total ( V-TOT = FRZ + POL + V-CT )' in line:
                assert is_cp
                scf_tot = float(line.split()[-1])
                almo_data_snap['tot_scf'] = scf_tot
            if 'Higher order relaxation ( HO = V-TOT - P-TOT )' in line:
                scf_ho = float(line.split()[-1])
                assert abs(scf_tot - rs_tot - scf_ho) < 1.0e-4
                almo_data_snap['ho_scf'] = scf_ho

    for k in almo_data_snap:
        almo_data_snap[k] /= KJ_TO_KCAL

    return almo_data_snap


def make_snapnum_to_bin_map():
    snapshots_filename = '/home/eric/Chemistry/calc.sgr/paper_02_CD_SC/inputs_freq/representative_snapshots_1qm'
    snapnum_to_bin_map = dict()
    with open(snapshots_filename) as fh:
        for line in fh:
            if line[0] == '#':
                binnum = int(line.split()[2])
            else:
                snapnum = int(line.strip())
                snapnum_to_bin_map[snapnum] = binnum
    return snapnum_to_bin_map

snapnum_to_bin_map = make_snapnum_to_bin_map()

# B3LYP/6-31G(d,p)
bin_to_weight_map = {
    1: 0.06060606,
    2: 0.24747475,
    3: 0.4010101,
    4: 0.22626263,
    5: 0.06464646,
}

sapt_headers_monomer = [
    '//         Monomer Basis SAPT        //',
]

sapt_headers_dimer = [
    '//               SAPT0               //',
    '//          Dimer Basis SAPT         //',
]

sapt_headers = sapt_headers_monomer + sapt_headers_dimer

sapt_bases = ('monomer', 'dimer')

def read_psi4_sapt_section(fi):
    sapt_single_basis_data = dict()

    line = ''
    while 'SAPT Results' not in line:
        line = next(fi)
    line = next(fi)
    assert '------' in line
    line = next(fi)
    assert 'Electrostatics' in line
    val_kcal_mol_electrostatics = float(re_number.findall(line)[1])
    sapt_single_basis_data['el'] = val_kcal_mol_electrostatics
    line = next(fi)
    assert 'Elst10,r' in line
    line = next(fi)
    assert line.strip() == ''
    line = next(fi)
    assert 'Exchange' in line
    val_kcal_mol_exchange = float(re_number.findall(line)[1])
    sapt_single_basis_data['exch'] = val_kcal_mol_exchange
    line = next(fi)
    assert 'Exch10' in line
    line = next(fi)
    assert 'Exch10(S^2)' in line
    line = next(fi)
    assert line.strip() == ''
    line = next(fi)
    assert 'Induction' in line
    line = next(fi)
    assert 'Ind20,r' in line
    val_kcal_mol_induction = float(re_number.findall(line)[1])
    sapt_single_basis_data['ind'] = val_kcal_mol_induction
    line = next(fi)
    assert 'Exch-Ind20,r' in line
    val_kcal_mol_exchange_induction = float(re_number.findall(line)[1])
    sapt_single_basis_data['exch-ind'] = val_kcal_mol_exchange_induction
    line = next(fi)
    assert 'delta HF,r (2)' in line
    val_kcal_mol_induction_delta_hf = float(re_number.findall(line)[1])
    sapt_single_basis_data['ind_HO'] = val_kcal_mol_induction_delta_hf
    line = next(fi)
    assert line.strip() == ''
    line = next(fi)
    assert 'Dispersion' in line
    line = next(fi)
    assert 'Disp20' in line
    val_kcal_mol_dispersion = float(re_number.findall(line)[1])
    sapt_single_basis_data['disp'] = val_kcal_mol_dispersion
    line = next(fi)
    assert 'Exch-Disp20' in line
    val_kcal_mol_exchange_dispersion = float(re_number.findall(line)[1])
    sapt_single_basis_data['exch-disp'] = val_kcal_mol_exchange_dispersion

    while 'Total SAPT0' not in line:
        line = next(fi)
    sapt0_total_kcal_mol_calculated = float(re_number.findall(line)[1])

    sapt0_total_kcal_mol = val_kcal_mol_electrostatics + \
                           val_kcal_mol_exchange + \
                           val_kcal_mol_induction + \
                           val_kcal_mol_exchange_induction + \
                           val_kcal_mol_induction_delta_hf + \
                           val_kcal_mol_dispersion + \
                           val_kcal_mol_exchange_dispersion

    # print(abs(sapt0_total_kcal_mol - sapt0_total_kcal_mol_calculated))
    assert abs(sapt0_total_kcal_mol - sapt0_total_kcal_mol_calculated) < 1.0e-7

    sapt_single_basis_data['total'] = sapt0_total_kcal_mol

    return sapt_single_basis_data


def read_psi4_sapt0(filename):
    sapt_data = dict()
    stub = os.path.splitext(os.path.basename(filename))[0]
    stub_tokens = stub.split('_')
    snapnum = int(stub_tokens[1])
    binnum = snapnum_to_bin_map[snapnum]
    weight = bin_to_weight_map[binnum]
    sapt_data['snapnum'] = snapnum
    sapt_data['weight'] = weight
    fi = make_file_iterator(filename)

    # Collect both dimer-centered and monomer-centered SAPT basis
    # data.
    for line in fi:
        # Dimer results always come before monomer results.
        if any(sapt_header in line for sapt_header in sapt_headers_dimer):
            sapt_data['dimer'] = read_psi4_sapt_section(fi)
        if any(sapt_header in line for sapt_header in sapt_headers_monomer):
            sapt_data['monomer'] = read_psi4_sapt_section(fi)
            break

    # Finally, check to see if a charge transfer (CT) calculation has
    # been performed.
    for line in fi:
        if 'SAPT Charge Transfer Analysis' in line:
            line = next(fi)
            assert list(set(line.strip())) == ['-']
            line = next(fi)
            # Asserts here comparing to induction values that were
            # parsed earlier?
            assert 'SAPT Induction (Dimer Basis)' in line
            line = next(fi)
            assert 'SAPT Induction (Monomer Basis)' in line
            line = next(fi)
            assert 'SAPT Charge Transfer' in line
            ct_kcal_mol = float(re_number.findall(line)[1])
            sapt_data['ct'] = ct_kcal_mol

    return sapt_data


def df_weighted_average(df, weights):
    # Have to do things manually with Pandas, it seems.
    weights = np.array(weights).reshape((len(weights), 1))
    weights_sum = np.sum(weights)
    return (weights * df.select_dtypes(include=[np.number])).sum(axis=0) / weights_sum
