from __future__ import division
from __future__ import print_function

import csv
import os

import subprocess as sp

import numpy as np

from mbe.xyz_operations import read_xyz

from mbe.examples.droplet import \
    (determine_fragment_grouping, make_fragments_from_grouping,
     distance_twopoint, fragment_centerofmass,
     get_n_closest_fragments)


if __name__ == '__main__':

    find_output = sp.check_output('find /home/eric/Chemistry/calc.sgr/paper_02_CD_SC/NewXYZFiles/ -name "*.xyz" | sort', shell=True).decode()
    xyzfilenames = sorted(find_output.splitlines())
    dim = 1001
    assert len(xyzfilenames) == dim

    possible_n_qm = list(range(1, 64 + 1))
    dropnums = {n_qm: [] for n_qm in possible_n_qm}
    n_qm_to_distance_from_COM_COM_combined = {n_qm: [] for n_qm in possible_n_qm}
    n_qm_to_distance_from_COM_COM_cation = {n_qm: [] for n_qm in possible_n_qm}
    n_qm_to_distance_from_COM_COM_anion = {n_qm: [] for n_qm in possible_n_qm}

    for xyzfilename in xyzfilenames:
        dropnum = os.path.splitext(xyzfilename.split('_')[-1])[0]
        print(dropnum)
        _, _, elements, coords = read_xyz(xyzfilename)
        assert elements[-3:] == ['C', 'O', 'O']
        grouping_anions, grouping_cations, grouping_CO2 = determine_fragment_grouping(elements)
        fragments_anions = make_fragments_from_grouping(elements, coords, grouping_anions)
        fragments_cations = make_fragments_from_grouping(elements, coords, grouping_cations)
        fragment_CO2 = make_fragments_from_grouping(elements, coords, grouping_CO2)
        fragments = fragments_anions + fragments_cations + fragment_CO2
        assert len(fragment_CO2) == 1
        fragment_CO2 = fragment_CO2[0]
        for n_qm in possible_n_qm:
            # Wouldn't it be nice if we didn't calculate the distance twice?
            closest_cations = get_n_closest_fragments(n_qm, fragment_CO2, fragments_cations, method='centerofmass')
            closest_anions = get_n_closest_fragments(n_qm, fragment_CO2, fragments_anions, method='centerofmass')
            max_distance_from_COM_COM_combined = -1.0e30
            max_distance_from_COM_COM_cation = -1.0e30
            max_distance_from_COM_COM_anion = -1.0e30
            COM_CO2 = fragment_centerofmass(fragment_CO2)
            for fragment_cation in closest_cations:
                COM_cation = fragment_centerofmass(fragment_cation)
                max_distance_from_COM_COM_cation = max(max_distance_from_COM_COM_cation,
                                                       distance_twopoint(COM_CO2, COM_cation))
            for fragment_anion in closest_anions:
                COM_anion = fragment_centerofmass(fragment_anion)
                max_distance_from_COM_COM_anion = max(max_distance_from_COM_COM_anion,
                                                      distance_twopoint(COM_CO2, COM_anion))
            max_distance_from_COM_COM_combined = max(max_distance_from_COM_COM_cation,
                                                     max_distance_from_COM_COM_anion)
            dropnums[n_qm].append(dropnum)
            n_qm_to_distance_from_COM_COM_combined[n_qm].append(max_distance_from_COM_COM_combined)
            n_qm_to_distance_from_COM_COM_cation[n_qm].append(max_distance_from_COM_COM_cation)
            n_qm_to_distance_from_COM_COM_anion[n_qm].append(max_distance_from_COM_COM_anion)

    row_header = [
        'dropnum',
        'distances_COM_COM_combined',
        'distances_COM_COM_cation',
        'distances_COM_COM_anion',
    ]
    for n_qm in possible_n_qm:
        csvfh = open('raw_distances_{}QM.csv'.format(n_qm), 'w')
        csvwriter = csv.writer(csvfh)
        csvwriter.writerow(row_header)
        assert len(dropnums[n_qm]) == len(n_qm_to_distance_from_COM_COM_combined[n_qm]) == len(n_qm_to_distance_from_COM_COM_cation[n_qm]) == len(n_qm_to_distance_from_COM_COM_anion[n_qm])
        dim = len(dropnums[n_qm])
        n_qm_array = np.empty(shape=(dim, 4))
        n_qm_array[:, 0] = dropnums[n_qm]
        n_qm_array[:, 1] = n_qm_to_distance_from_COM_COM_combined[n_qm]
        n_qm_array[:, 2] = n_qm_to_distance_from_COM_COM_cation[n_qm]
        n_qm_array[:, 3] = n_qm_to_distance_from_COM_COM_anion[n_qm]
        csvwriter.writerows(n_qm_array.tolist())
        csvfh.close()

    csvfh = open('mean_box_sizes.csv', 'w')
    csvwriter = csv.writer(csvfh)
    row_header = [
        'n_qm',
        'mean_COM_COM_combined',
        'mean_COM_COM_cation',
        'mean_COM_COM_anion',
    ]
    csvwriter.writerow(row_header)
    print(row_header)
    for n_qm in possible_n_qm:
        mean_COM_COM_combined = np.mean(n_qm_to_distance_from_COM_COM_combined[n_qm])
        mean_COM_COM_cation = np.mean(n_qm_to_distance_from_COM_COM_cation[n_qm])
        mean_COM_COM_anion = np.mean(n_qm_to_distance_from_COM_COM_anion[n_qm])
        row = [
            n_qm,
            mean_COM_COM_combined,
            mean_COM_COM_cation,
            mean_COM_COM_anion
        ]
        csvwriter.writerow(row)
        print(row)
    csvfh.close()
