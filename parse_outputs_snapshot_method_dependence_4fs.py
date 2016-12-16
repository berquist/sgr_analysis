#!/usr/bin/env python3

from __future__ import print_function
from __future__ import division

import os
import pickle

from collections import OrderedDict

from analysis_utils import get_CO2_frequencies_d
from analysis_utils import get_dipoles_d
from analysis_utils import get_outputfiles_from_path
from analysis_utils import filter_n_mm_into_dict

methods = OrderedDict([
    # ('blyp', 'BLYP'),
    # ('tpss', 'TPSS'),
    ('b3lyp', 'B3LYP'),
    # ('wb97x-d', r'$\omega$B97X-D'),
    # ('hf', 'HF'),
    # ('ri-mp2', 'RI-MP2'),
])

basis_sets = OrderedDict([
    # ('6-31gdp', '6-31G(d,p)'),
    # ('cc-pvtz', 'cc-pVTZ'),
    ('lp', '6-311++G(d,p)'),
])

CO2_types = OrderedDict([
    ('flex', 'flex'),
    ('rigid', 'rigid'),
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
    # 1. rigid / flexible CO2
    # 2. method
    # 3. basis set
    # 4. n_qm
    # 5. n_mm

    BASEDIR = "/home/eric/Chemistry/calc.sgr/paper_02_CD_SC/snapshot_method_dependence_4fs"

    if args.file_operation == "save":
        print("Trying to find output files...")
        outputfiles = dict()
        for CO2_type in CO2_types:
            outputfiles[CO2_type] = dict()
            for method in methods:
                outputfiles[CO2_type][method] = dict()
                for basis_set in basis_sets:
                    outputfiles[CO2_type][method][basis_set] = dict()
                    # hacky hack hack
                    for n_qm in (0, ):
                        outputfiles[CO2_type][method][basis_set][n_qm] = get_outputfiles_from_path(os.path.join(BASEDIR, "inputs_freq_{CO2_type}_0qm_{method}_{basis_set}".format(**locals())))
                        outputfiles[CO2_type][method][basis_set][n_qm] = filter_n_mm_into_dict(outputfiles[CO2_type][method][basis_set][n_qm])
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
        for CO2_type in CO2_types:
            frequencies[CO2_type] = dict()
            intensities[CO2_type] = dict()
            snapnums_f[CO2_type] = dict()
            for method in methods:
                frequencies[CO2_type][method] = dict()
                intensities[CO2_type][method] = dict()
                snapnums_f[CO2_type][method] = dict()
                for basis_set in outputfiles[CO2_type][method]:
                    frequencies[CO2_type][method][basis_set] = dict()
                    intensities[CO2_type][method][basis_set] = dict()
                    snapnums_f[CO2_type][method][basis_set] = dict()
                    for n_qm in (0, ):
                        m_bs_qm_frequencies, m_bs_qm_intensities, m_bs_qm_snapnums_f = get_CO2_frequencies_d(outputfiles[CO2_type][method][basis_set][n_qm], do_4fs=True)
                        frequencies[CO2_type][method][basis_set][n_qm] = m_bs_qm_frequencies
                        intensities[CO2_type][method][basis_set][n_qm] = m_bs_qm_intensities
                        snapnums_f[CO2_type][method][basis_set][n_qm] = m_bs_qm_snapnums_f
        print("Parsing dipoles...")
        for CO2_type in CO2_types:
            dipoles[CO2_type] = dict()
            snapnums_d[CO2_type] = dict()
            for method in methods:
                dipoles[CO2_type][method] = dict()
                snapnums_d[CO2_type][method] = dict()
                for basis_set in outputfiles[CO2_type][method]:
                    dipoles[CO2_type][method][basis_set] = dict()
                    snapnums_d[CO2_type][method][basis_set] = dict()
                    for n_qm in (0, ):
                        m_bs_qm_dipoles, m_bs_qm_snapnums_d = get_dipoles_d(outputfiles[CO2_type][method][basis_set][n_qm])
                        dipoles[CO2_type][method][basis_set][n_qm] = m_bs_qm_dipoles
                        snapnums_d[CO2_type][method][basis_set][n_qm] = m_bs_qm_snapnums_d
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
