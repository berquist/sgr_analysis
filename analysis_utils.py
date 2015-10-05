from __future__ import print_function
from __future__ import division

import os
import re

import scipy.stats as sps

from cclib.parser import ccopen

from scripts.find_CO2_frequencies import find_CO2_mode_indices


slice_lambda = lambda x, start, end: x >= start and x < end


def make_n_mm_dict():
    """So we can avoid having to check every loop iteration later on. Make
    all possible keys even if we won't use the bigger values.
    """

    d = dict()

    keys = list(range(0, 18, 2)) + [32, 64, 128] + [256, 255, 254, 253]
    for k in keys:
        d[k] = list()

    return d


def filter_n_mm_into_dict(outputfilenames):

    d = make_n_mm_dict()

    for outputfilename in outputfilenames:
        try:
            n_mm = int(re.search(r'_(\d+)mm', outputfilename).groups()[0])
        # Due to a silly shell mishap: `rename _0 _ ./*` is a bad idea, kids.
        except AttributeError:
            n_mm = 0
        d[n_mm].append(outputfilename)

    return d


def get_CO2_frequencies(outputfilenames):

    snapnums = []
    CO2_frequencies = []
    CO2_intensities = []

    for outputfilename in outputfilenames:
        print("Parsing frequencies from {}".format(outputfilename))

        job = ccopen(outputfilename)
        try:
            data = job.parse()
        except:
            # Is this the right control flow statement?
            continue
        # geometry = data.atomcoords[-1]
        # atoms = data.atomnos
        # start_indices = find_CO2_atom_indices(atoms, geometry)
        # assert isinstance(start_indices, list)
        try:
            vibfreqs = data.vibfreqs
            vibdisps = data.vibdisps
            vibirs = data.vibirs
        except AttributeError:
            # Is this the correct control flow statement?
            continue
        # Assumption!
        # start_index = start_indices[0]
        # Assumption?
        start_index = len(data.atomnos) - 3
        mode_indices = find_CO2_mode_indices(start_index, vibdisps, thresh=0.50)
        # mode_indices = [2]
        # freqs = [vibfreqs[modeidx] for modeidx in mode_indices]
        # freqs = filter(lambda x: x > 0.0, freqs)
        # print(freqs)
        # Let's only take the last one...
        # print(outputfilename)
        CO2_frequencies.append(vibfreqs[mode_indices[-1]])
        CO2_intensities.append(vibirs[mode_indices[-1]])
        snapnum = int(re.search(r'drop_(\d+)_', outputfilename).groups()[0])
        snapnums.append(snapnum)

    return CO2_frequencies, CO2_intensities, snapnums


def get_dipoles(outputfilenames):

    dipoles = []
    snapnums = []

    for outputfilename in outputfilenames:

        with open(outputfilename) as outputfile:
            print("Parsing dipole from {}".format(outputfilename))
            for line in outputfile:
                if 'Dipole Moment (Debye)' in line:
                    line = next(outputfile)
                    dipole = list(map(float, line.split()[1::2]))
                    dipoles.append(dipole)
                    snapnum = int(re.search(r'drop_(\d+)_', outputfilename).groups()[0])
                    snapnums.append(snapnum)
                    # Only take the first one! This avoids problems
                    # when parsing numerical frequency runs.
                    break

    return dipoles, snapnums


def get_dipoles_supersystem(outputfilenames):

    pass


def mangle_dict_keys(d):
    """Not all 'maximum MM' calculations are going to have 256 pairs, if
    there are ionic liquid pairs treated explicitly; they'll take away
    from that number, because there are 256 *total* pairs in a box.

    Make it appear in a dictionary that these are all the same by
    'having' 256 MM pairs.
    """

    nd = dict()

    for k in d:
        if d[k] == []:
            d.pop(k)

    bad_keys = (253, 254, 255)
    for n_qm in d:
        nd[n_qm] = dict()
        for n_mm in d[n_qm]:
            if n_mm in bad_keys and len(d[n_qm][n_mm]) > 0:
                nd[n_qm][256] = d[n_qm][n_mm]
            else:
                nd[n_qm][n_mm] = d[n_qm][n_mm]

    return nd


def get_CO2_frequencies_d(filename_dict):

    frequencies_dict = make_n_mm_dict()
    intensities_dict = make_n_mm_dict()
    snapnums_dict = make_n_mm_dict()

    for n_mm in filename_dict:
        if len(filename_dict[n_mm]) > 0:
            frequencies_dict[n_mm], intensities_dict[n_mm], snapnums_dict[n_mm] = get_CO2_frequencies(filename_dict[n_mm])

    return frequencies_dict, intensities_dict, snapnums_dict


def get_dipoles_d(filename_dict):

    dipoles_dict = make_n_mm_dict()
    snapnums_dict = make_n_mm_dict()

    for n_mm in filename_dict:
        if len(filename_dict[n_mm]) > 0:
            dipoles_dict[n_mm], snapnums_dict[n_mm] = get_dipoles(filename_dict[n_mm])

    return dipoles_dict, snapnums_dict


def get_single_snapshot_results(snapnum, snapnums_dict, results_dict):

    snap_results_dict = dict()

    for n_qm in sorted(results_dict.keys()):
        snap_results_dict[n_qm] = dict()
        for n_mm in sorted(results_dict[n_qm].keys()):
            try:
                snapidx = snapnums_dict[n_qm][n_mm].index(snapnum)
                snap_results_dict[n_qm][n_mm] = [results_dict[n_qm][n_mm][snapidx]]
            # not in list?
            except ValueError:
                snap_results_dict[n_qm][n_mm] = []

    return snap_results_dict


def get_eda_covp_totals(outputfilepath):
    """Given a path to an output file, return the totals for each fragment from the COVP analysis.
    The first element of the tuple is the energy contribution, the second element is the
    number of millielectrons transferred.
    """

    searchstr = "#   Delta E(Alpha)    Delta E(Beta)  Delta Q(Alpha)   Delta Q(Beta)"
    remove_paren_stuff = lambda x: float(x[:-8])

    with open(outputfilepath) as outputfile:
        for line in outputfile:
            if searchstr in line:
                while line[:4] != " Tot":
                    line = next(outputfile)
                f_1_to_2 = tuple(map(remove_paren_stuff, line.split()[1::2]))
                line = next(outputfile)
                while line[:4] != " Tot":
                    line = next(outputfile)
                f_2_to_1 = tuple(map(remove_paren_stuff, line.split()[1::2]))

    return f_1_to_2, f_2_to_1


def sort(snapnums_dict, results_dict):

    assert snapnums_dict.keys() == results_dict.keys()
    for k in results_dict:
        sorting_indices = [i[0] for i in sorted(enumerate(snapnums_dict[k]),
                                                key=lambda x: x[1])]
        sorted_results = [i[1] for i in sorted(zip(sorting_indices, results_dict[k]),
                                               key=lambda x: x[0])]
        sorted_snapnums = [i[1] for i in sorted(zip(sorting_indices, snapnums_dict[k]),
                                                key=lambda x: x[0])]
        # Why is this commented out?
        # assert sorted_snapnums == list(range(min(snapnums_dict[k]), max(snapnums_dict[k]) + 1))
        snapnums_dict[k] = sorted_snapnums
        results_dict[k] = sorted_results

    return


def get_outputfiles_from_path(path, ext=".out"):
    """Walk the directory tree to find all potential output files.
    """

    outputfiles = []

    for (root, dirs, files) in os.walk(path, followlinks=True):
        for f in files:
            if f.endswith(ext):
                outputfiles.append(os.path.join(root, f))

    return sorted(outputfiles)


def pprint_lengths(d):

    for upper_key in sorted(d.keys()):
        print(upper_key,
              [len(d[upper_key][lower_key])
               for lower_key in sorted(d[upper_key].keys())])

    return


def pprint_linregress(x, y):

    slope, intercept, rval, pval, stderr = sps.linregress(x, y)
    print(" slope:     {:f}".format(slope))
    print(" intercept: {:f}".format(intercept))
    rsq = rval**2
    print(" rval2:     {:f}".format(rsq))

    return slope, intercept, rsq
