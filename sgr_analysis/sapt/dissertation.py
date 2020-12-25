"""dissertation.py

This isn't how it was originally done, but this recreates the tables
and figures from jp6b09489_si_013.xlsx for use in my dissertation
(LaTeX).

Reading/"writing" code adapted from ./sapt_vs_almo.py. The tables
there were modified into an Excel file and a bunch of correlation
plots were added by hand. Reproduce them.
"""


from glob import glob
import os.path

import numpy as np
import pandas as pd

from helpers import (bin_to_weight_map, read_psi4_sapt0,
                     read_qchem_eda, snapnum_to_bin_map)

if __name__ == '__main__':

    root_dir = '/home/eric/Chemistry/calc.sgr/paper_02_CD_SC/sapt/sapt_vs_almo/'
    root_dir_sapt = '/home/eric/Chemistry/calc.sgr/paper_02_CD_SC/sapt/631gdp/ct/'

    filenames_hf_cp = sorted(glob(root_dir + '*hf_6-31gss_cp.out'))
    filenames_hf_nocp = sorted(glob(root_dir + '*hf_6-31gss_nocp.out'))
    filenames_b3lyp_cp = sorted(glob(root_dir + '*b3lyp_6-31gss_cp.out'))
    filenames_b3lyp_nocp = sorted(glob(root_dir + '*b3lyp_6-31gss_nocp.out'))

    filenames_sapt = sorted(glob(root_dir_sapt + '*.out'))

    batches = (
        filenames_hf_cp,
        filenames_hf_nocp,
        filenames_b3lyp_cp,
        filenames_b3lyp_nocp,
    )

    map_is_cp = {
        'cp': True,
        'nocp': False,
    }

    # Read in ALMO calculations.
    snapnums = set()
    almo_data = dict()

    for batch in batches:
        for filename in batch:

            stub = os.path.splitext(os.path.basename(filename))[0]
            print(stub)
            stub_tokens = stub.split('_')
            snapnum = int(stub_tokens[1])
            snapnums.add(snapnum)
            method = stub_tokens[6] + '_' + stub_tokens[7]
            cp_flag = stub_tokens[-1]
            is_cp = map_is_cp[cp_flag]

            if method not in almo_data:
                almo_data[method] = dict()
            if cp_flag not in almo_data[method]:
                almo_data[method][cp_flag] = dict()

            almo_data_snap = read_qchem_eda(filename, is_cp=is_cp)

            almo_data[method][cp_flag][snapnum] = almo_data_snap

    snapnums = sorted(snapnums)
    weights = [bin_to_weight_map[snapnum_to_bin_map[snapnum]]
               for snapnum in snapnums]

    # Read in SAPT calculations.
    sapt_data = dict()
    snapnums_sapt = set()

    for filename in filenames_sapt:
        stub = os.path.splitext(os.path.basename(filename))[0]
        print(stub)
        stub_tokens = stub.split('_')
        snapnum = int(stub_tokens[1])
        assert snapnum in snapnums
        snapnums_sapt.add(snapnum)
        sapt_data_snap = read_psi4_sapt0(filename)
        sapt_data[snapnum] = sapt_data_snap

    snapnums_sapt = sorted(snapnums_sapt)
    weights_sapt = [bin_to_weight_map[snapnum_to_bin_map[snapnum]]
                    for snapnum in snapnums_sapt]

    # header_not_data = [
    #     'snapnum',
    #     'qc_method',
    #     # 'bsse_corr',
    # ]
    # header = [
    #     'frz',
    #     'pol',
    #     # 'del_rs',
    #     # 'bsse_rs',
    #     # 'ct_rs',
    #     'del_scf',
    #     'bsse_scf',
    #     'ct_scf',
    #     'tot_scf',
    #     'ho_scf',
    # ]
    # total_header = header_not_data + header

    header = (
        'snapnum',
        'ALMO/B3LYP/6-31G(d,p) E_frz',
        'ALMO/B3LYP/6-31G(d,p) E_pol',
        'ALMO/B3LYP/6-31G(d,p) E_del (SCF)',
        'ALMO/B3LYP/6-31G(d,p) E_CT (SCF)',
        'ALMO/B3LYP/6-31G(d,p) E_tot (SCF)',
        'ALMO/HF/6-31G(d,p) E_frz',
        'ALMO/HF/6-31G(d,p) E_pol',
        'ALMO/HF/6-31G(d,p) E_del (SCF)',
        'ALMO/HF/6-31G(d,p) E_CT (SCF)',
        'ALMO/HF/6-31G(d,p) E_tot (SCF)',
        'SAPT/6-31G(d,p)/MCBS E_el',
        'SAPT/6-31G(d,p)/MCBS E_exch',
        'SAPT/6-31G(d,p)/MCBS "E_frz"',
        'SAPT/6-31G(d,p)/MCBS E_ind',
        'SAPT/6-31G(d,p)/MCBS E_ind-exch',
        'SAPT/6-31G(d,p)/MCBS E_ind_HO',
        'SAPT/6-31G(d,p)/MCBS "E_pol"',
        'SAPT/6-31G(d,p)/MCBS E_disp',
        'SAPT/6-31G(d,p)/MCBS E_disp-exch',
        'SAPT/6-31G(d,p)/MCBS E_totdisp',
        'SAPT/6-31G(d,p)/MCBS E_CT',
        'SAPT/6-31G(d,p)/MCBS E_tot',
        'SAPT/6-31G(d,p)/DCBS E_el',
        'SAPT/6-31G(d,p)/DCBS E_exch',
        'SAPT/6-31G(d,p)/DCBS "E_frz"',
        'SAPT/6-31G(d,p)/DCBS E_ind',
        'SAPT/6-31G(d,p)/DCBS E_ind-exch',
        'SAPT/6-31G(d,p)/DCBS E_ind_HO',
        'SAPT/6-31G(d,p)/DCBS "E_pol"',
        'SAPT/6-31G(d,p)/DCBS E_disp',
        'SAPT/6-31G(d,p)/DCBS E_disp-exch',
        'SAPT/6-31G(d,p)/DCBS E_totdisp',
        'SAPT/6-31G(d,p)/DCBS E_CT',
        'SAPT/6-31G(d,p)/DCBS E_tot',
    )

    rows = []

    for snapnum in snapnums_sapt:
        monomer_frz = sapt_data[snapnum]['monomer']['el'] + sapt_data[snapnum]['monomer']['exch']
        monomer_pol = sapt_data[snapnum]['monomer']['ind'] + sapt_data[snapnum]['monomer']['exch-ind'] + sapt_data[snapnum]['monomer']['ind_HO']
        monomer_tot_disp = sapt_data[snapnum]['monomer']['disp'] + sapt_data[snapnum]['monomer']['exch-disp']
        dimer_frz = sapt_data[snapnum]['dimer']['el'] + sapt_data[snapnum]['dimer']['exch']
        dimer_pol = sapt_data[snapnum]['dimer']['ind'] + sapt_data[snapnum]['dimer']['exch-ind'] + sapt_data[snapnum]['dimer']['ind_HO']
        dimer_tot_disp = sapt_data[snapnum]['dimer']['disp'] + sapt_data[snapnum]['dimer']['exch-disp']
        line = [
            snapnum,
            almo_data['b3lyp_6-31gss']['cp'][snapnum]['frz'],
            almo_data['b3lyp_6-31gss']['cp'][snapnum]['pol'],
            almo_data['b3lyp_6-31gss']['cp'][snapnum]['del_scf'],
            almo_data['b3lyp_6-31gss']['cp'][snapnum]['ct_scf'],
            almo_data['b3lyp_6-31gss']['cp'][snapnum]['tot_scf'],
            almo_data['hf_6-31gss']['cp'][snapnum]['frz'],
            almo_data['hf_6-31gss']['cp'][snapnum]['pol'],
            almo_data['hf_6-31gss']['cp'][snapnum]['del_scf'],
            almo_data['hf_6-31gss']['cp'][snapnum]['ct_scf'],
            almo_data['hf_6-31gss']['cp'][snapnum]['tot_scf'],
            sapt_data[snapnum]['monomer']['el'],
            sapt_data[snapnum]['monomer']['exch'],
            monomer_frz,
            sapt_data[snapnum]['monomer']['ind'],
            sapt_data[snapnum]['monomer']['exch-ind'],
            sapt_data[snapnum]['monomer']['ind_HO'],
            monomer_pol,
            sapt_data[snapnum]['monomer']['disp'],
            sapt_data[snapnum]['monomer']['exch-disp'],
            monomer_tot_disp,
            sapt_data[snapnum]['ct'],
            sapt_data[snapnum]['monomer']['total'],
            sapt_data[snapnum]['dimer']['el'],
            sapt_data[snapnum]['dimer']['exch'],
            dimer_frz,
            sapt_data[snapnum]['dimer']['ind'],
            sapt_data[snapnum]['dimer']['exch-ind'],
            sapt_data[snapnum]['dimer']['ind_HO'],
            dimer_pol,
            sapt_data[snapnum]['dimer']['disp'],
            sapt_data[snapnum]['dimer']['exch-disp'],
            dimer_tot_disp,
            sapt_data[snapnum]['ct'],
            sapt_data[snapnum]['dimer']['total'],
        ]
        rows.append(line)

    start = 1
    vals = np.array(rows, dtype=float)[:, start:]
    average_unweighted = np.average(vals, axis=0, weights=None)
    average_weighted = np.average(vals, axis=0, weights=weights_sapt)
    # csvwriter.writerow(['average (unweighted)'] + average_unweighted.tolist())
    # csvwriter.writerow(['average (weighted)'] + average_weighted.tolist())
