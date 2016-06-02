from __future__ import division
from __future__ import print_function

import csv
import os

import subprocess as sp

import numpy as np

if __name__ == '__main__':

    possible_n_qm = list(range(1, 64 + 1))
    dropnums = np.arange(1, 1001 + 1)
    distances = dict()

    for n_qm in possible_n_qm:
        with open('raw_distances_{}QM.csv'.format(n_qm)) as csvfile:
            csvreader = csv.reader(csvfile)
            rows = []
            header = next(csvreader)
            assert header == [
                'dropnum',
                'distances_COM_COM_combined',
                'distances_COM_COM_cation',
                'distances_COM_COM_anion',
            ]
            for row in csvreader:
                rows.append(row)
        rows = np.array(sorted(rows), dtype=float)
        distances[n_qm] = rows

    for n_qm in possible_n_qm:
        assert dropnums.all() == distances[n_qm][:, 0].all()
        assert distances[n_qm].shape == (1001, 4)
    dim = rows.shape[0]
    distances_COM_COM_combined = np.empty(shape=(dim, len(possible_n_qm) + 1))
    distances_COM_COM_cation = np.empty(shape=(dim, len(possible_n_qm) + 1))
    distances_COM_COM_anion = np.empty(shape=(dim, len(possible_n_qm) + 1))
    distances_COM_COM_combined[:, 0] = dropnums
    distances_COM_COM_cation[:, 0] = dropnums
    distances_COM_COM_anion[:, 0] = dropnums
    for n_qm in possible_n_qm:
        distances_COM_COM_combined[:, n_qm] = distances[n_qm][:, 1]
        distances_COM_COM_cation[:, n_qm] = distances[n_qm][:, 2]
        distances_COM_COM_anion[:, n_qm] = distances[n_qm][:, 3]

    means_combined = [np.mean(distances_COM_COM_combined[:, n_qm]) for n_qm in possible_n_qm]
    means_cation = [np.mean(distances_COM_COM_cation[:, n_qm]) for n_qm in possible_n_qm]
    means_anion = [np.mean(distances_COM_COM_anion[:, n_qm]) for n_qm in possible_n_qm]
    print(means_combined)
    print(means_cation)
    print(means_anion)

    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    ax.plot(means_combined, marker='o', label='combined')
    ax.plot(means_cation, marker='o', label='cation')
    ax.plot(means_anion, marker='o', label='anion')

    ax.set_xticklabels(possible_n_qm)
    ax.set_xlabel('# QM pairs')
    ax.set_ylabel(r'distance ($\AA$)')

    ax.legend(loc='best', fancybox=True, framealpha=0.50)

    fig.savefig('box_sizes.pdf', bbox_inches='tight')

    plt.close('all')

    fig, ax = plt.subplots()

    ax.plot(means_combined, possible_n_qm, marker='o', label='combined')
    ax.plot(means_cation, possible_n_qm, marker='o', label='cation')
    ax.plot(means_anion, possible_n_qm, marker='o', label='anion')

    ax.set_xlabel(r'distance ($\AA$)')
    ax.set_ylabel('# QM pairs')

    ax.legend(loc='best', fancybox=True, framealpha=0.50)

    fig.savefig('box_sizes_inverted.pdf', bbox_inches='tight')
