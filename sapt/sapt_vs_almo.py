"""sapt_vs_almo.py: Parse a set of ALMO output files for their
interaction energy components, writing to a CSV file."""

from __future__ import division
from __future__ import print_function

from collections import OrderedDict
from glob import glob
import csv
import os.path

import numpy as np

from helpers import (bin_to_weight_map, read_psi4_sapt0,
                     read_qchem_eda, snapnum_to_bin_map)

from summary import method_basis_outer_product


if __name__ == '__main__':

    root_dir = '/home/eric/Chemistry/calc.sgr/paper_02_CD_SC/sapt/sapt_vs_almo/'
    root_dir_sapt = '/home/eric/Chemistry/calc.sgr/paper_02_CD_SC/sapt/631gdp/ct/'
    root_dir_paper1 = '/home/eric/Chemistry/calc.sgr/paper_02_CD_SC/sapt/almo_eda_paper1/'

    filenames_hf_cp = sorted(glob(root_dir + '*hf_6-31gss_cp.out'))
    filenames_hf_nocp = sorted(glob(root_dir + '*hf_6-31gss_nocp.out'))
    filenames_b3lyp_cp = sorted(glob(root_dir + '*b3lyp_6-31gss_cp.out'))
    filenames_b3lyp_nocp = sorted(glob(root_dir + '*b3lyp_6-31gss_nocp.out'))
    filenames_b3lypd2_cp = sorted(glob(root_dir + '*b3lyp-d2_6-31gss_cp.out'))
    filenames_b3lypd2_nocp = sorted(glob(root_dir + '*b3lyp-d2_6-31gss_nocp.out'))
    filenames_b3lypd3_cp = sorted(glob(root_dir + '*b3lyp-d3_6-31gss_cp.out'))
    filenames_b3lypd3_nocp = sorted(glob(root_dir + '*b3lyp-d3_6-31gss_nocp.out'))
    filenames_wb97xd_pople_cp = sorted(glob(root_dir + '*wb97x-d_6-31gss_cp.out'))
    filenames_wb97xd_pople_nocp = sorted(glob(root_dir + '*wb97x-d_6-31gss_nocp.out'))
    filenames_wb97mv_pople_cp = sorted(glob(root_dir + '*wb97m-v_6-31gss_cp.out'))
    filenames_wb97mv_pople_nocp = sorted(glob(root_dir + '*wb97m-v_6-31gss_nocp.out'))
    filenames_wb97xd_dunning_cp = sorted(glob(root_dir + '*wb97x-d_cc-pvtz_cp.out'))
    filenames_wb97xd_dunning_nocp = sorted(glob(root_dir + '*wb97x-d_cc-pvtz_nocp.out'))
    filenames_wb97mv_dunning_cp = sorted(glob(root_dir + '*wb97m-v_cc-pvtz_cp.out'))
    filenames_wb97mv_dunning_nocp = sorted(glob(root_dir + '*wb97m-v_cc-pvtz_nocp.out'))

    filenames_sapt = sorted(glob(root_dir_sapt + '*.out'))
    filenames_paper1 = sorted(glob(root_dir_paper1 + '*.out'))

    batches = (
        filenames_hf_cp,
        filenames_hf_nocp,
        filenames_b3lyp_cp,
        filenames_b3lyp_nocp,
        filenames_b3lypd2_cp,
        filenames_b3lypd2_nocp,
        filenames_b3lypd3_cp,
        filenames_b3lypd3_nocp,
        filenames_wb97xd_pople_cp,
        filenames_wb97xd_pople_nocp,
        filenames_wb97mv_pople_cp,
        filenames_wb97mv_pople_nocp,
        filenames_wb97xd_dunning_cp,
        filenames_wb97xd_dunning_nocp,
        filenames_wb97mv_dunning_cp,
        filenames_wb97mv_dunning_nocp,
    )

    map_is_cp = {
        'cp': True,
        'nocp': False,
    }

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

    ## Start by writing all the ALMO data to a CSV file.

    csvfh = open('data_almo.csv', 'w')
    csvwriter = csv.writer(csvfh)

    header_not_data = [
        'snapnum',
        'qc_method',
        # 'bsse_corr',
    ]
    header = [
        'frz',
        'pol',
        # 'del_rs',
        # 'bsse_rs',
        # 'ct_rs',
        'del_scf',
        'bsse_scf',
        'ct_scf',
        'tot_scf',
        'ho_scf',
    ]
    total_header = header_not_data + header
    csvwriter.writerow(total_header)

    # Write only the CP-corrected values.
    for method_basis_pair in method_basis_outer_product:
        qc_method = '_'.join(method_basis_pair)
        if qc_method in almo_data:
            for snapnum in sorted(almo_data[qc_method]['cp']):
                row = [almo_data[qc_method]['cp'][snapnum][column_title]
                       for column_title in header]
                row = [snapnum, qc_method] + row
                csvwriter.writerow(row)

    csvfh.close()

    ##########

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

    csvfh = open('data_all.csv', 'w')
    csvwriter = csv.writer(csvfh)

    # 20 fields
    header = (
        'snapnum',
        'ALMO/wB97M-V/cc-pVTZ E_frz',
        'ALMO/wB97M-V/cc-pVTZ E_pol',
        'ALMO/wB97M-V/cc-pVTZ E_del (SCF)',
        'ALMO/wB97M-V/cc-pVTZ E_CT (SCF)',
        'ALMO/wB97X-D/cc-pVTZ E_frz',
        'ALMO/wB97X-D/cc-pVTZ E_pol',
        'ALMO/wB97X-D/cc-pVTZ E_del (SCF)',
        'ALMO/wB97X-D/cc-pVTZ E_CT (SCF)',
        'ALMO/wB97M-V/6-31G(d,p) E_frz',
        'ALMO/wB97M-V/6-31G(d,p) E_pol',
        'ALMO/wB97M-V/6-31G(d,p) E_del (SCF)',
        'ALMO/wB97M-V/6-31G(d,p) E_CT (SCF)',
        'ALMO/wB97X-D/6-31G(d,p) E_frz',
        'ALMO/wB97X-D/6-31G(d,p) E_pol',
        'ALMO/wB97X-D/6-31G(d,p) E_del (SCF)',
        'ALMO/wB97X-D/6-31G(d,p) E_CT (SCF)',
        'ALMO/B3LYP-D2/6-31G(d,p) E_frz',
        'ALMO/B3LYP-D2/6-31G(d,p) E_pol',
        'ALMO/B3LYP-D2/6-31G(d,p) E_del (SCF)',
        'ALMO/B3LYP-D2/6-31G(d,p) E_CT (SCF)',
        'ALMO/B3LYP-D3/6-31G(d,p) E_frz',
        'ALMO/B3LYP-D3/6-31G(d,p) E_pol',
        'ALMO/B3LYP-D3/6-31G(d,p) E_del (SCF)',
        'ALMO/B3LYP-D3/6-31G(d,p) E_CT (SCF)',
        'ALMO/B3LYP/6-31G(d,p) E_frz',
        'ALMO/B3LYP/6-31G(d,p) E_pol',
        'ALMO/B3LYP/6-31G(d,p) E_del (SCF)',
        'ALMO/B3LYP/6-31G(d,p) E_CT (SCF)',
        'ALMO/HF/6-31G(d,p) E_frz',
        'ALMO/HF/6-31G(d,p) E_pol',
        'ALMO/HF/6-31G(d,p) E_del (SCF)',
        'ALMO/HF/6-31G(d,p) E_CT (SCF)',
        'SAPT/6-31G(d,p)/DCBS E_el',
        'SAPT/6-31G(d,p)/DCBS E_exch',
        'SAPT/6-31G(d,p)/DCBS E_ind',
        'SAPT/6-31G(d,p)/DCBS E_ind-exch',
        'SAPT/6-31G(d,p)/DCBS E_ind_HO',
        'SAPT/6-31G(d,p)/DCBS E_disp',
        'SAPT/6-31G(d,p)/DCBS E_disp-exch',
        'SAPT/6-31G(d,p)/DCBS E_CT',
    )
    csvwriter.writerow(header)

    rows = []

    for snapnum in snapnums_sapt:
        line = [
            snapnum,
            almo_data['wb97m-v_cc-pvtz']['cp'][snapnum]['frz'],
            almo_data['wb97m-v_cc-pvtz']['cp'][snapnum]['pol'],
            almo_data['wb97m-v_cc-pvtz']['cp'][snapnum]['del_scf'],
            almo_data['wb97m-v_cc-pvtz']['cp'][snapnum]['ct_scf'],
            almo_data['wb97x-d_cc-pvtz']['cp'][snapnum]['frz'],
            almo_data['wb97x-d_cc-pvtz']['cp'][snapnum]['pol'],
            almo_data['wb97x-d_cc-pvtz']['cp'][snapnum]['del_scf'],
            almo_data['wb97x-d_cc-pvtz']['cp'][snapnum]['ct_scf'],
            almo_data['wb97m-v_6-31gss']['cp'][snapnum]['frz'],
            almo_data['wb97m-v_6-31gss']['cp'][snapnum]['pol'],
            almo_data['wb97m-v_6-31gss']['cp'][snapnum]['del_scf'],
            almo_data['wb97m-v_6-31gss']['cp'][snapnum]['ct_scf'],
            almo_data['wb97x-d_6-31gss']['cp'][snapnum]['frz'],
            almo_data['wb97x-d_6-31gss']['cp'][snapnum]['pol'],
            almo_data['wb97x-d_6-31gss']['cp'][snapnum]['del_scf'],
            almo_data['wb97x-d_6-31gss']['cp'][snapnum]['ct_scf'],
            almo_data['b3lyp-d2_6-31gss']['cp'][snapnum]['frz'],
            almo_data['b3lyp-d2_6-31gss']['cp'][snapnum]['pol'],
            almo_data['b3lyp-d2_6-31gss']['cp'][snapnum]['del_scf'],
            almo_data['b3lyp-d2_6-31gss']['cp'][snapnum]['ct_scf'],
            almo_data['b3lyp-d3_6-31gss']['cp'][snapnum]['frz'],
            almo_data['b3lyp-d3_6-31gss']['cp'][snapnum]['pol'],
            almo_data['b3lyp-d3_6-31gss']['cp'][snapnum]['del_scf'],
            almo_data['b3lyp-d3_6-31gss']['cp'][snapnum]['ct_scf'],
            almo_data['b3lyp_6-31gss']['cp'][snapnum]['frz'],
            almo_data['b3lyp_6-31gss']['cp'][snapnum]['pol'],
            almo_data['b3lyp_6-31gss']['cp'][snapnum]['del_scf'],
            almo_data['b3lyp_6-31gss']['cp'][snapnum]['ct_scf'],
            almo_data['hf_6-31gss']['cp'][snapnum]['frz'],
            almo_data['hf_6-31gss']['cp'][snapnum]['pol'],
            almo_data['hf_6-31gss']['cp'][snapnum]['del_scf'],
            almo_data['hf_6-31gss']['cp'][snapnum]['ct_scf'],
            sapt_data[snapnum]['dimer']['el'],
            sapt_data[snapnum]['dimer']['exch'],
            sapt_data[snapnum]['dimer']['ind'],
            sapt_data[snapnum]['dimer']['exch-ind'],
            sapt_data[snapnum]['dimer']['ind_HO'],
            sapt_data[snapnum]['dimer']['disp'],
            sapt_data[snapnum]['dimer']['exch-disp'],
            sapt_data[snapnum]['ct'],
        ]
        rows.append(line)
        csvwriter.writerow(line)

    start = 1
    vals = np.array(rows, dtype=float)[:, start:]
    average_unweighted = np.average(vals, axis=0, weights=None)
    average_weighted = np.average(vals, axis=0, weights=weights_sapt)
    csvwriter.writerow(['average (unweighted)'] + average_unweighted.tolist())
    csvwriter.writerow(['average (weighted)'] + average_weighted.tolist())

    csvfh.close()

    # ##########

    # almo_data_paper1 = dict()

    # csvfh = open('almo_eda_paper1_all.csv', 'w')
    # csvwriter = csv.writer(csvfh)

    # # 5 fields
    # header = (
    #     'system',
    #     'ALMO/B3LYP E_frz',
    #     'ALMO/B3LYP E_pol',
    #     'ALMO/B3LYP E_del (SCF)',
    #     'ALMO/B3LYP E_CT (SCF)',
    # )
    # csvwriter.writerow(header)

    # for filename in filenames_paper1:
    #     stub = os.path.splitext(os.path.basename(filename))[0]
    #     stub_tokens = stub.split('_')
    #     system = '_'.join(stub_tokens[1:])
    #     almo_data_system = read_qchem_eda(filename, is_cp=True)
    #     line = [
    #         system,
    #         almo_data_system['frz'],
    #         almo_data_system['pol'],
    #         almo_data_system['del_scf'],
    #         almo_data_system['ct_scf'],
    #     ]
    #     csvwriter.writerow(line)

    # csvfh.close()
