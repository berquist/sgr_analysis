#!/usr/bin/env python3

"""parse_outputs.py: Parse outputfiles for certain values, like
frequencies and dipole moments, and store them in pickled
dictionaries.
"""

from __future__ import print_function
from __future__ import division

import os
import pickle

from analysis_utils import get_outputfiles_from_path
from analysis_utils import filter_n_mm_into_dict
from analysis_utils import get_CO2_frequencies_d
from analysis_utils import get_dipoles_d


def getargs():

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("file_operation",
                        choices=("none", "save", "read"),
                        help="""What operation should be done for finding output files?""")
    parser.add_argument("--dir-to-search", default=".")
    parser.add_argument("--logfilename", default="outputfiles.log")

    parser.add_argument("parse_operation",
                        choices=("none", "save", "read"),
                        help="""What operation should be done for parsing the output files?""")
    parser.add_argument("--parse-frequencies", action="store_true")
    parser.add_argument("--parse-dipoles", action="store_true")
    parser.add_argument("--parse-eda", action="store_true")

    parser.add_argument("--debug", action="store_true")

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = getargs()

    outputfiles = []

    # We'll sort outputfiles according to these types and the number of IL
    # pairs treated quantum mechanically (as opposed to using point
    # charges).
    outputs_freq_0qm = []
    outputs_freq_1qm = []
    outputs_freq_2qm = []
    outputs_freq_3qm = []
    outputs_freq_noCT_1qm = []
    outputs_freq_noCT_2qm = []
    outputs_freq_noCT_3qm = []
    outputs_eda_covp_1qm = []
    outputs_eda_covp_2qm = []
    outputs_eda_covp_3qm = []

    if args.file_operation == "save":
        print("Trying to find output files...")
        outputfiles = get_outputfiles_from_path(args.dir_to_search)
        with open(args.logfilename, 'w') as logfile:
            for outputfilename in outputfiles:
                logfile.write("{}\n".format(outputfilename))
    elif args.file_operation == "read":
        print("Reading list of output files from: {}".format(os.path.abspath(args.logfilename)))
        with open(args.logfilename) as logfile:
            content = logfile.read()
        outputfiles = content.splitlines()
    elif args.file_operation == "none":
        pass
    else:
        raise Exception

    if args.debug:
        print("len(outputfiles)")
        print(len(outputfiles))

    print("Sorting output files into categories...")
    for outputfile in outputfiles:

        if 'freq_0qm' in outputfile:
            outputs_freq_0qm.append(outputfile)
        elif 'freq_1qm' in outputfile:
            outputs_freq_1qm.append(outputfile)
        elif 'freq_2qm' in outputfile:
            outputs_freq_2qm.append(outputfile)
        elif 'freq_3qm' in outputfile:
            outputs_freq_3qm.append(outputfile)
        elif 'freq_noCT_1qm' in outputfile:
            outputs_freq_noCT_1qm.append(outputfile)
        elif 'freq_noCT_2qm' in outputfile:
            outputs_freq_noCT_2qm.append(outputfile)
        elif 'freq_noCT_3qm' in outputfile:
            outputs_freq_noCT_3qm.append(outputfile)
        elif 'eda_covp_1qm' in outputfile:
            outputs_eda_covp_1qm.append(outputfile)
        elif 'eda_covp_2qm' in outputfile:
            outputs_eda_covp_2qm.append(outputfile)
        elif 'eda_covp_3qm' in outputfile:
            outputs_eda_covp_3qm.append(outputfile)
        else:
            continue

    if args.debug:
        print("Filtered lengths:")
        print(len(outputs_freq_0qm))
        print(len(outputs_freq_1qm))
        print(len(outputs_freq_2qm))
        print(len(outputs_freq_3qm))
        print(len(outputs_freq_noCT_1qm))
        print(len(outputs_freq_noCT_2qm))
        print(len(outputs_freq_noCT_3qm))
        print(len(outputs_eda_covp_1qm))
        print(len(outputs_eda_covp_2qm))
        print(len(outputs_eda_covp_3qm))

    print("Filtering filenames by # of MM pairs...")
    filenames_freq_0qm = filter_n_mm_into_dict(outputs_freq_0qm)
    filenames_freq_1qm = filter_n_mm_into_dict(outputs_freq_1qm)
    filenames_freq_2qm = filter_n_mm_into_dict(outputs_freq_2qm)
    filenames_freq_3qm = filter_n_mm_into_dict(outputs_freq_3qm)
    filenames_freq_noCT_1qm = filter_n_mm_into_dict(outputs_freq_noCT_1qm)
    filenames_freq_noCT_2qm = filter_n_mm_into_dict(outputs_freq_noCT_2qm)
    filenames_freq_noCT_3qm = filter_n_mm_into_dict(outputs_freq_noCT_3qm)
    filenames_eda_covp_1qm = filter_n_mm_into_dict(outputs_eda_covp_1qm)
    filenames_eda_covp_2qm = filter_n_mm_into_dict(outputs_eda_covp_2qm)
    filenames_eda_covp_3qm = filter_n_mm_into_dict(outputs_eda_covp_3qm)

    if args.parse_operation == "save":

        print("Extracting valuable information from outputs...")

        print("Parsing frequencies/intensities...")
        frequencies_CO2_0qm, intensities_CO2_0qm, snapnums_frequencies_CO2_0qm = get_CO2_frequencies_d(filenames_freq_0qm)
        frequencies_CO2_1qm, intensities_CO2_1qm, snapnums_frequencies_CO2_1qm = get_CO2_frequencies_d(filenames_freq_1qm)
        frequencies_CO2_2qm, intensities_CO2_2qm, snapnums_frequencies_CO2_2qm = get_CO2_frequencies_d(filenames_freq_2qm)
        frequencies_CO2_3qm, intensities_CO2_3qm, snapnums_frequencies_CO2_3qm = get_CO2_frequencies_d(filenames_freq_3qm)
        print("Parsing frequencies/intensities (w/o charge transfer)...")
        frequencies_noCT_CO2_1qm, intensities_noCT_CO2_1qm, snapnums_frequencies_noCT_CO2_1qm = get_CO2_frequencies_d(filenames_freq_noCT_1qm)
        frequencies_noCT_CO2_2qm, intensities_noCT_CO2_2qm, snapnums_frequencies_noCT_CO2_2qm = get_CO2_frequencies_d(filenames_freq_noCT_2qm)
        frequencies_noCT_CO2_3qm, intensities_noCT_CO2_3qm, snapnums_frequencies_noCT_CO2_3qm = get_CO2_frequencies_d(filenames_freq_noCT_3qm)
        frequencies_CO2_d = {
            0: frequencies_CO2_0qm,
            1: frequencies_CO2_1qm,
            2: frequencies_CO2_2qm,
            3: frequencies_CO2_3qm,
        }
        intensities_CO2_d = {
            0: intensities_CO2_0qm,
            1: intensities_CO2_1qm,
            2: intensities_CO2_2qm,
            3: intensities_CO2_3qm,
        }
        frequencies_noCT_CO2_d = {
            1: frequencies_noCT_CO2_1qm,
            2: frequencies_noCT_CO2_2qm,
            3: frequencies_noCT_CO2_3qm,
        }
        intensities_noCT_CO2_d = {
            1: intensities_noCT_CO2_1qm,
            2: intensities_noCT_CO2_2qm,
            3: intensities_noCT_CO2_3qm,
        }
        snapnums_frequencies_d = {
            0: snapnums_frequencies_CO2_0qm,
            1: snapnums_frequencies_CO2_1qm,
            2: snapnums_frequencies_CO2_2qm,
            3: snapnums_frequencies_CO2_3qm,
        }
        snapnums_frequencies_noCT_d = {
            1: snapnums_frequencies_noCT_CO2_1qm,
            2: snapnums_frequencies_noCT_CO2_2qm,
            3: snapnums_frequencies_noCT_CO2_3qm,
        }
        with open('frequencies.pypickle', 'wb') as picklefile:
            pickle.dump(frequencies_CO2_d, picklefile)
        with open('intensities.pypickle', 'wb') as picklefile:
            pickle.dump(intensities_CO2_d, picklefile)
        with open('frequencies_noCT.pypickle', 'wb') as picklefile:
            pickle.dump(frequencies_noCT_CO2_d, picklefile)
        with open('intensities_noCT.pypickle', 'wb') as picklefile:
            pickle.dump(intensities_noCT_CO2_d, picklefile)
        with open('snapnums_frequencies.pypickle', 'wb') as picklefile:
            pickle.dump(snapnums_frequencies_d, picklefile)
        with open('snapnums_frequencies_noCT.pypickle', 'wb') as picklefile:
            pickle.dump(snapnums_frequencies_noCT_d, picklefile)

        print("Parsing dipoles...")
        dipoles_0qm, snapnums_dipoles_0qm = get_dipoles_d(filenames_freq_0qm)
        dipoles_1qm, snapnums_dipoles_1qm = get_dipoles_d(filenames_freq_1qm)
        dipoles_2qm, snapnums_dipoles_2qm = get_dipoles_d(filenames_freq_2qm)
        dipoles_3qm, snapnums_dipoles_3qm = get_dipoles_d(filenames_freq_3qm)
        dipoles_d = {
            0: dipoles_0qm,
            1: dipoles_1qm,
            2: dipoles_2qm,
            3: dipoles_3qm,
        }
        snapnums_dipoles_d = {
            0: snapnums_dipoles_0qm,
            1: snapnums_dipoles_1qm,
            2: snapnums_dipoles_2qm,
            3: snapnums_dipoles_3qm,
        }
        with open('dipoles.pypickle', 'wb') as picklefile:
            pickle.dump(dipoles_d, picklefile)
        with open('snapnums_dipoles.pypickle', 'wb') as picklefile:
            pickle.dump(snapnums_dipoles_d, picklefile)

    elif args.parse_operation == "read":
        with open('frequencies.pypickle', 'rb') as picklefile:
            frequencies_CO2_d = pickle.load(picklefile)
        with open('intensities.pypickle', 'rb') as picklefile:
            intensities_CO2_d = pickle.load(picklefile)
        with open('frequencies_noCT.pypickle', 'rb') as picklefile:
            frequencies_noCT_CO2_d = pickle.load(picklefile)
        with open('intensities_noCT.pypickle', 'rb') as picklefile:
            intensities_noCT_CO2_d = pickle.load(picklefile)
        with open('dipoles.pypickle', 'rb') as picklefile:
            dipoles_d = pickle.load(picklefile)
        with open('snapnums_frequencies.pypickle', 'rb') as picklefile:
            snapnums_frequencies_d = pickle.load(picklefile)
        with open('snapnums_frequencies_noCT.pypickle', 'rb') as picklefile:
            snapnums_frequencies_noCT_d = pickle.load(picklefile)
        with open('snapnums_dipoles.pypickle', 'rb') as picklefile:
            snapnums_dipoles_d = pickle.load(picklefile)
    elif args.parse_operation == "none":
        pass
    else:
        raise Exception
