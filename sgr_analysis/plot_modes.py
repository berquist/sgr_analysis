#!/usr/bin/env python3

from __future__ import print_function

import os

from glob import glob

from cclib.parser import ccopen

from model_hamiltonian_frequencies import \
    (distance, bond_angle)
from analysis_utils import pad_left_zeros



def rename_step_files_clyde(filenames):
    maxlen = 0
    for oldfilename in filenames:
        stub = os.path.basename(oldfilename)
        newlen = len(stub[4:-4])
        if newlen > maxlen:
            maxlen = newlen
    newfilenames = []
    for oldfilename in filenames:
        directory, basename = os.path.split(oldfilename)
        stub, ext = os.path.splitext(basename)
        filenum = pad_left_zeros(stub[4:], maxlen)
        newfilename = os.path.join(directory, 'step_{}{}'.format(filenum, ext))
        os.rename(oldfilename, newfilename)
        print(oldfilename + ' -> ' + newfilename)
        newfilenames.append(newfilename)
    return newfilenames


def rename_step_files_eric(filenames):
    maxlen = 0
    for oldfilename in filenames:
        stub = os.path.basename(oldfilename)
        newlen = len(stub[5:-4])
        if newlen > maxlen:
            maxlen = newlen
    newfilenames = []
    for oldfilename in filenames:
        directory, basename = os.path.split(oldfilename)
        stub, ext = os.path.splitext(basename)
        filenum = pad_left_zeros(stub[5:], maxlen)
        newfilename = os.path.join(directory, 'step_{}{}'.format(filenum, ext))
        # os.rename(oldfilename, newfilename)
        print(oldfilename + ' -> ' + newfilename)
        newfilenames.append(newfilename)
    return newfilenames


if __name__ == "__main__":

    root_directory = "/home/eric/Chemistry/calc.sgr/paper_02_CD_SC/clyde_histogram/freq_dvr_no_almo_0qm/Snap1"

    # filenames = sorted(glob(os.path.join(root_directory, "step*")))
    # rename_step_files_eric(filenames)

    filenames = sorted(glob(os.path.join(root_directory, "step_*.out")))

    frequencies = []
    intensities = []
    geometries = []
    geometry_params = dict()
    geometry_params['l1'] = []
    geometry_params['l2'] = []
    geometry_params['l12'] = []

    for filename in filenames:

        job = ccopen(filename)
        data = job.parse()

        vibfreqs = data.vibfreqs
        vibirs = data.vibirs
        geometry = data.atomcoords[0]

        assert len(vibfreqs) == len(vibirs) == 3

        frequencies.append(vibfreqs)
        intensities.append(vibirs)
        geometries.append(geometry)

        C, O1, O2 = 0, 1, 2

        d_C_O1 = distance(geometry[C], geometry[O1])
        d_C_O2 = distance(geometry[C], geometry[O2])
        d_O1_O2 = distance(geometry[O1], geometry[O2])

        bond_sum = d_C_O1 + d_C_O2
        bond_difference = abs(d_C_O1 - d_C_O2)

        angle = bond_angle(geometry[O1], geometry[C], geometry[O2])
        theta = 180.0 - angle

        l1 = d_C_O1
        l2 = d_C_O2
        o12 = d_O1_O2
        l12 = bond_sum
        dl = bond_difference

        geometry_params['l1'].append(l1)
        geometry_params['l2'].append(l2)
        geometry_params['l12'].append(l12)

    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt

    ticks = list(range(1, len(frequencies) + 1))

    fig, ax_f = plt.subplots()

    line1 = ax_f.plot([f[0] for f in frequencies], linewidth=0.5, label='mode 1')
    line2 = ax_f.plot([f[1] for f in frequencies], linewidth=0.5, label='mode 2')
    line3 = ax_f.plot([f[2] for f in frequencies], linewidth=0.5, label='mode 3')

    ax_g = ax_f.twinx()

    line4 = ax_g.plot(geometry_params['l12'], linestyle=':', color='black', label='total bond length')

    ax_f.set_xlim((min(ticks), max(ticks)))

    ax_f.set_xlabel('step #')
    ax_f.set_ylabel(r'frequency (cm$^{-1}$)')
    ax_g.set_ylabel(r'bond length ($\AA$)')

    lines = line1 + line2 + line3 + line4
    labels = [line.get_label() for line in lines]
    ax_f.legend(lines, labels,
                loc='best', fancybox=True, framealpha=0.50, fontsize='x-small')

    fig.savefig('modes.pdf', bbox_inches='tight')

    # decorate, sort, undecorate
    sorted_zip = sorted(list(zip(geometry_params['l12'],
                                 geometry_params['l1'],
                                 geometry_params['l2'],
                                 frequencies)))
    geometry_params['l12'] = [e[0] for e in sorted_zip]
    geometry_params['l1'] = [e[1] for e in sorted_zip]
    geometry_params['l2'] = [e[2] for e in sorted_zip]
    frequencies = [e[3] for e in sorted_zip]

    fig, ax_f = plt.subplots()

    line1 = ax_f.plot([f[0] for f in frequencies], linewidth=0.5, label='mode 1')
    line2 = ax_f.plot([f[1] for f in frequencies], linewidth=0.5, label='mode 2')
    line3 = ax_f.plot([f[2] for f in frequencies], linewidth=0.5, label='mode 3')

    ax_g = ax_f.twinx()

    line4 = ax_g.plot(geometry_params['l12'], linestyle=':', color='black', label='total bond length')

    ax_f.set_xlim((min(ticks), max(ticks)))

    ax_f.set_xlabel('step #')
    ax_f.set_ylabel(r'frequency (cm$^{-1}$)')
    ax_g.set_ylabel(r'bond length ($\AA$)')

    lines = line1 + line2 + line3 + line4
    labels = [line.get_label() for line in lines]
    ax_f.legend(lines, labels,
                loc='best', fancybox=True, framealpha=0.50, fontsize='x-small')

    fig.savefig('modes_sorted.pdf', bbox_inches='tight')

    # decorate, sort, undecorate
    # sorted_zip = sorted(list(zip(geometry_params['l1'],
    #                              geometry_params['l2'],
    #                              geometry_params['l12'],
    #                              frequencies)))
    # geometry_params['l1'] = [e[0] for e in sorted_zip]
    # geometry_params['l2'] = [e[1] for e in sorted_zip]
    # geometry_params['l12'] = [e[2] for e in sorted_zip]
    # frequencies = [e[3] for e in sorted_zip]

    # fig, ax_f = plt.subplots()

    # line1 = ax_f.plot([f[0] for f in frequencies], linewidth=0.5, label='mode 1')
    # line2 = ax_f.plot([f[1] for f in frequencies], linewidth=0.5, label='mode 2')
    # line3 = ax_f.plot([f[2] for f in frequencies], linewidth=0.5, label='mode 3')

    # ax_g = ax_f.twinx()

    # line4 = ax_g.plot(geometry_params['l1'], linestyle=':', color='black', label='C-O length 1')
    # line5 = ax_g.plot(geometry_params['l2'], linestyle=':', color='orange', label='C-O length 2')

    # ax_f.set_xlim((min(ticks), max(ticks)))

    # ax_f.set_xlabel('step #')
    # ax_f.set_ylabel(r'frequency (cm$^{-1}$)')
    # ax_g.set_ylabel(r'bond length ($\AA$)')

    # lines = line1 + line2 + line3 + line4 + line5
    # labels = [line.get_label() for line in lines]
    # ax_f.legend(lines, labels,
    #             loc='best', fancybox=True, framealpha=0.50, fontsize='x-small')

    # fig.savefig('modes_sorted2.pdf', bbox_inches='tight')
