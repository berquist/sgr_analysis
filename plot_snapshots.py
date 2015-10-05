#!/usr/bin/env python3

from __future__ import print_function
from __future__ import division

import pickle

import numpy as np

from analysis_utils import mangle_dict_keys
# from analysis_utils import sort


if __name__ == "__main__":

    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt

    with open('frequencies.pypickle', 'rb') as picklefile:
        frequencies_CO2_d = pickle.load(picklefile)
    with open('snapnums_frequencies.pypickle', 'rb') as picklefile:
        snapnums_frequencies_d = pickle.load(picklefile)

    frequencies_CO2_d = mangle_dict_keys(frequencies_CO2_d)
    snapnums_frequencies_d = mangle_dict_keys(snapnums_frequencies_d)

    assert frequencies_CO2_d.keys() == snapnums_frequencies_d.keys()

    # sort(snapnums_frequencies_d, frequencies_CO2_d)

    # n_qm, then n_mm
    # print(snapnums_frequencies_d[2][0])
    # print(frequencies_CO2_d[2][0])

    z2 = sorted([(sn, f) for (sn, f) in zip(snapnums_frequencies_d[2][0], frequencies_CO2_d[2][0])])
    z2_sn, z2_f = [p[0] for p in z2], [p[1] for p in z2]
    z1 = sorted([(sn, f) for (sn, f) in zip(snapnums_frequencies_d[1][0], frequencies_CO2_d[1][0])
                 if sn in z2_sn])
    z1_sn, z1_f = [p[0] for p in z1], [p[1] for p in z1]

    print("np.mean(z2_f)")
    print(np.mean(z2_f))
    print("np.mean(z1_f)")
    print(np.mean(z1_f))

    fig, ax = plt.subplots()

    ax.plot(z1_sn, z1_f, marker='o', label='1 QM')
    ax.plot(z2_sn, z2_f, marker='o', label='2 QM')

    ax.set_xlabel('snapshot #')
    ax.set_ylabel(r'$\nu_3$ frequency (cm$^{-1}$)')

    ax.legend(fancybox=True, loc='best', framealpha=0.50)

    fig.savefig('plot_snapshots.pdf', bbox_inches='tight')
