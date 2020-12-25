#!/usr/bin/env python3


import os
import pickle

from collections import OrderedDict

from analysis_utils import get_CO2_frequencies_d
from analysis_utils import get_dipoles_d
from analysis_utils import get_outputfiles_from_path
from analysis_utils import filter_n_mm_into_dict

methods = OrderedDict([
    ('blyp', 'BLYP'),
    ('tpss', 'TPSS'),
    ('b3lyp', 'B3LYP'),
    ('wb97x-d', r'$\omega$B97X-D'),
    ('hf', 'HF'),
    ('ri-mp2', 'RI-MP2'),
])

basis_sets = OrderedDict([
    ('6-31gdp', '6-31G(d,p)'),
    ('cc-pvtz', 'cc-pVTZ'),
])

def getargs():

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("file_operation",
                        choices=("none", "save", "read"),
                        help="""What operation should be done for finding output files?""")

    parser.add_argument("parse_operation",
                        choices=("none", "save", "read"),
                        help="""What operation should be done for parsing the output files?""")

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = getargs()

    # Dictionary layering:
    # 1. method
    # 2. basis set
    # 3. n_qm
    # 4. n_mm

    if args.file_operation == "save":
        print("Trying to find output files...")
        basedir = "/home/eric/Chemistry/calc.sgr/paper_02_CD_SC/snapshot_method_dependence"
        outputfiles = dict()
        for method in methods:
            outputfiles[method] = dict()
            for basis_set in basis_sets:
                outputfiles[method][basis_set] = dict()
                # hacky hack hack
                for n_qm in (0, ):
                    outputfiles[method][basis_set][n_qm] = get_outputfiles_from_path(os.path.join(basedir, "inputs_freq_0qm_{method}_{basis_set}".format(method=method, basis_set=basis_set)))
                    outputfiles[method][basis_set][n_qm] = filter_n_mm_into_dict(outputfiles[method][basis_set][n_qm])
        with open('outputfiles.pypickle', 'wb') as picklefile:
            pickle.dump(outputfiles, picklefile)

    elif args.file_operation == "read":
        print("Reading list of output files from: {}".format(os.path.abspath("outputfiles.pypickle")))
        with open("outputfiles.pypickle", "rb") as picklefile:
            outputfiles = pickle.load(picklefile)

    elif args.file_operation == "none":
        pass
    else:
        raise Exception

    if args.parse_operation == "save":
        print("Extracting valuable information from outputs...")
        frequencies = dict()
        intensities = dict()
        snapnums_f = dict()
        dipoles = dict()
        snapnums_d = dict()
        print("Parsing frequencies/intensities...")
        for method in methods:
            frequencies[method] = dict()
            intensities[method] = dict()
            snapnums_f[method] = dict()
            for basis_set in outputfiles[method]:
                frequencies[method][basis_set] = dict()
                intensities[method][basis_set] = dict()
                snapnums_f[method][basis_set] = dict()
                for n_qm in (0, ):
                    m_bs_qm_frequencies, m_bs_qm_intensities, m_bs_qm_snapnums_f = get_CO2_frequencies_d(outputfiles[method][basis_set][n_qm])
                    frequencies[method][basis_set][n_qm] = m_bs_qm_frequencies
                    intensities[method][basis_set][n_qm] = m_bs_qm_intensities
                    snapnums_f[method][basis_set][n_qm] = m_bs_qm_snapnums_f
        print("Parsing dipoles...")
        for method in methods:
            dipoles[method] = dict()
            snapnums_d[method] = dict()
            for basis_set in outputfiles[method]:
                dipoles[method][basis_set] = dict()
                snapnums_d[method][basis_set] = dict()
                for n_qm in (0, ):
                     m_bs_qm_dipoles, m_bs_qm_snapnums_d = get_dipoles_d(outputfiles[method][basis_set][n_qm])
                     dipoles[method][basis_set][n_qm] = m_bs_qm_dipoles
                     snapnums_d[method][basis_set][n_qm] = m_bs_qm_snapnums_d
        with open('frequencies.pypickle', 'wb') as picklefile:
            pickle.dump(frequencies, picklefile)
        with open('intensities.pypickle', 'wb') as picklefile:
            pickle.dump(intensities, picklefile)
        with open('dipoles.pypickle', 'wb') as picklefile:
            pickle.dump(dipoles, picklefile)
        with open('snapnums_frequencies.pypickle', 'wb') as picklefile:
            pickle.dump(snapnums_f, picklefile)
        with open('snapnums_dipoles.pypickle', 'wb') as picklefile:
            pickle.dump(snapnums_d, picklefile)

    elif args.parse_operation == "read":
        print("Reading frequency data from: {}".format(os.path.abspath('frequencies.pypickle')))
        with open('frequencies.pypickle', 'rb') as picklefile:
            frequencies = pickle.load(picklefile)
        print("Reading intensity data from: {}".format(os.path.abspath('intensities.pypickle')))
        with open('intensities.pypickle', 'rb') as picklefile:
            intensities = pickle.load(picklefile)
        print("Reading dipole data from: {}".format(os.path.abspath('dipoles.pypickle')))
        with open('dipoles.pypickle', 'rb') as picklefile:
            dipoles = pickle.load(picklefile)
        print("Reading snapshot number data from: {}".format(os.path.abspath('snapnums_frequencies.pypickle')))
        with open('snapnums_frequencies.pypickle', 'rb') as picklefile:
            snapnums_f = pickle.load(picklefile)
        print("Reading snapshot number data from: {}".format(os.path.abspath('snapnums_dipoles.pypickle')))
        with open('snapnums_dipoles.pypickle', 'rb') as picklefile:
            snapnums_d = pickle.load(picklefile)

    elif args.parse_operation == "none":
        pass
    else:
        raise Exception
