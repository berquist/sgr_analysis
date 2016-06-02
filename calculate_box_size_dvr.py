from __future__ import division
from __future__ import print_function

import subprocess as sp

import numpy as np

from scripts.utils import make_file_iterator
import scripts.periodic_table as pt

from mbe.examples.droplet import \
    (determine_fragment_grouping, make_fragments_from_grouping,
     distance_atomic_shortest, distance_atomic_longest,
     distance_twopoint, fragment_centerofmass)

def parse_molecule_from_inputfile(inputfilename):
    all_symbols = []
    all_coords = []
    fi = make_file_iterator(inputfilename)
    line = ''
    while '$molecule' not in line.strip().lower():
        line = next(fi)
    line = next(fi)
    assert len(line.split()) == 2
    line = next(fi)
    while '$end' not in line.strip().lower():
        chomp = line.split()
        assert len(chomp) == 4
        symbol = chomp[0]
        coords = [float(x) for x in chomp[1:]]
        all_symbols.append(symbol)
        all_coords.append(coords)
        line = next(fi)
    return all_symbols, all_coords


def distance_point_1_furthest_atom_2(point, fragment2):
    """Return the largest distance between a given point and all the atoms
    in a given fragment.
    """
    distances = [distance_twopoint(point, c2) for c2 in fragment2.coords]
    return max(distances)


if __name__ == '__main__':

    find_output = sp.check_output('find . -wholename "*Snap*.in"', shell=True).decode()
    inputfilenames = sorted(find_output.splitlines())

    possible_n_qm = list(range(2, 6 + 1))
    n_qm_to_distance_from_COM_COM = {n_qm: [] for n_qm in possible_n_qm}
    n_qm_to_distance_from_COM_COM_cation = {n_qm: [] for n_qm in possible_n_qm}
    n_qm_to_distance_from_COM_COM_anion = {n_qm: [] for n_qm in possible_n_qm}

    for inputfilename in inputfilenames:
        print(inputfilename)
        elements, coords = parse_molecule_from_inputfile(inputfilename)
        assert elements[:3] == ['C', 'O', 'O']
        grouping_anions, grouping_cations, grouping_CO2 = determine_fragment_grouping(elements)
        grouping = grouping_CO2 + grouping_anions + grouping_cations
        fragments = make_fragments_from_grouping(elements, coords, grouping)
        fragments_anions = make_fragments_from_grouping(elements, coords, grouping_anions)
        fragments_cations = make_fragments_from_grouping(elements, coords, grouping_cations)
        fragment_CO2 = make_fragments_from_grouping(elements, coords, grouping_CO2)
        assert len(fragment_CO2) == 1
        fragment_CO2 = fragment_CO2[0]
        # fragment_CO2 = fragments[0]
        n_qm = len(fragments[1:]) // 2
        print(n_qm)
        max_distance_from_COM_COM = -1.0e30
        max_distance_from_COM_COM_cation = -1.0e30
        max_distance_from_COM_COM_anion = -1.0e30
        COM_CO2 = fragment_centerofmass(fragment_CO2)
        # print('looping over cations + anions')
        for fragment_IL in fragments[1:]:
            COM_IL = fragment_centerofmass(fragment_IL)
            # print(COM_IL)
            max_distance_from_COM_COM = max(max_distance_from_COM_COM,
                                            distance_twopoint(COM_CO2, COM_IL))
        # print('looping over cations')
        for fragment_cation in fragments_cations:
            COM_cation = fragment_centerofmass(fragment_cation)
            # print(COM_cation)
            max_distance_from_COM_COM_cation = max(max_distance_from_COM_COM_cation,
                                                   distance_twopoint(COM_CO2, COM_cation))
        # print('looping over anions')
        for fragment_anion in fragments_anions:
            COM_anion = fragment_centerofmass(fragment_anion)
            # print(COM_anion)
            max_distance_from_COM_COM_anion = max(max_distance_from_COM_COM_anion,
                                                  distance_twopoint(COM_CO2, COM_anion))
        n_qm_to_distance_from_COM_COM[n_qm].append(max_distance_from_COM_COM)
        n_qm_to_distance_from_COM_COM_cation[n_qm].append(max_distance_from_COM_COM_cation)
        n_qm_to_distance_from_COM_COM_anion[n_qm].append(max_distance_from_COM_COM_anion)

    import csv
    csvfh = open('mean_box_sizes.csv', 'w')
    csvwriter = csv.writer(csvfh)
    row_header = [
        'n_qm',
        'mean_COM_COM',
        'mean_COM_COM_cation',
        'mean_COM_COM_anion',
    ]
    csvwriter.writerow(row_header)
    print(row_header)
    for n_qm in possible_n_qm:
        mean_COM_COM = np.mean(n_qm_to_distance_from_COM_COM[n_qm])
        mean_COM_COM_cation = np.mean(n_qm_to_distance_from_COM_COM_cation[n_qm])
        mean_COM_COM_anion = np.mean(n_qm_to_distance_from_COM_COM_anion[n_qm])
        row = [
            n_qm,
            mean_COM_COM,
            mean_COM_COM_cation,
            mean_COM_COM_anion
        ]
        csvwriter.writerow(row)
        print(row)
    csvfh.close()

    row_header = [
        'distances_COM_COM',
        'distances_COM_COM_cation',
        'distances_COM_COM_anion',
    ]
    dim = 1440
    for n_qm in possible_n_qm:
        csvfh = open('raw_distances_{}QM.csv'.format(n_qm), 'w')
        csvwriter = csv.writer(csvfh)
        csvwriter.writerow(row_header)
        assert len(n_qm_to_distance_from_COM_COM[n_qm]) == dim
        assert len(n_qm_to_distance_from_COM_COM_cation[n_qm]) == dim
        assert len(n_qm_to_distance_from_COM_COM_anion[n_qm]) == dim
        n_qm_array = np.empty(shape=(dim, 3))
        n_qm_array[:, 0] = n_qm_to_distance_from_COM_COM[n_qm]
        n_qm_array[:, 1] = n_qm_to_distance_from_COM_COM_cation[n_qm]
        n_qm_array[:, 2] = n_qm_to_distance_from_COM_COM_anion[n_qm]
        csvwriter.writerows(n_qm_array.tolist())
        csvfh.close()
