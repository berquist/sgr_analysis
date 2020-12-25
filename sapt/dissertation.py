"""This isn't how it was originally done, but this recreates the
tables and figures from"""

from __future__ import division
from __future__ import print_function

from glob import glob

if __name__ == '__main__':

    root_dir = '/home/eric/Chemistry/calc.sgr/paper_02_CD_SC/sapt/sapt_vs_almo/'
    root_dir_sapt = '/home/eric/Chemistry/calc.sgr/paper_02_CD_SC/sapt/631gdp/ct/'
    root_dir_paper1 = '/home/eric/Chemistry/calc.sgr/paper_02_CD_SC/sapt/almo_eda_paper1/'

    filenames_hf_cp = sorted(glob(root_dir + '*hf_6-31gss_cp.out'))
    filenames_hf_nocp = sorted(glob(root_dir + '*hf_6-31gss_nocp.out'))
    filenames_b3lyp_cp = sorted(glob(root_dir + '*b3lyp_6-31gss_cp.out'))
    filenames_b3lyp_nocp = sorted(glob(root_dir + '*b3lyp_6-31gss_nocp.out'))

    filenames_sapt = sorted(glob(root_dir_sapt + '*.out'))
    filenames_paper1 = sorted(glob(root_dir_paper1 + '*.out'))

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
