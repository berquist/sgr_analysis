import os

import numpy as np

from glob import glob

from cclib.parser import ccopen

from sgr_analysis.model_hamiltonian_frequencies import (distance, bond_angle)


def pad_left_zeros(num, maxwidth):
    """Pad the given number with zeros on the left until the
    total length is maxwidth, returning it as a string.
    """

    numwidth = len(str(num))
    if numwidth < maxwidth:
        numzeros = maxwidth - numwidth
        numstr = (numzeros * '0') + str(num)
    else:
        numstr = str(num)
    return numstr


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


def extract(filenames, equilibrium_length=None):

    frequencies = []
    intensities = []
    geometry_params = dict()
    geometry_params['geometry'] = []
    geometry_params['l1'] = []
    geometry_params['l2'] = []
    geometry_params['l12'] = []
    geometry_params['theta'] = []
    if equilibrium_length:
        geometry_params['l1_diff'] = []
        geometry_params['l2_diff'] = []
        geometry_params['l12_diff'] = []
        geometry_params['o12_diff'] = []

    for filename in filenames:

        try:

            job = ccopen(filename)
            data = job.parse()

            vibfreqs = data.vibfreqs
            vibirs = data.vibirs
            geometry = data.atomcoords[0]

            assert len(vibfreqs) == len(vibirs) >= 3

            frequencies.append(vibfreqs)
            intensities.append(vibirs)

            C, O1, O2 = 0, 1, 2

            d_C_O1 = distance(geometry[C], geometry[O1])
            d_C_O2 = distance(geometry[C], geometry[O2])
            d_O1_O2 = distance(geometry[O1], geometry[O2])

            bond_sum = d_C_O1 + d_C_O2
            # bond_difference = abs(d_C_O1 - d_C_O2)

            angle = bond_angle(geometry[O1], geometry[C], geometry[O2])
            theta = 180.0 - angle

            l1 = d_C_O1
            l2 = d_C_O2
            o12 = d_O1_O2
            l12 = bond_sum
            # dl = bond_difference

            if equilibrium_length:
                l1_diff = l1 - equilibrium_length
                l2_diff = l2 - equilibrium_length
                o12_diff = o12 - (2 * equilibrium_length)
                l12_diff = l12 - (2 * equilibrium_length)
                geometry_params['l1_diff'].append(l1_diff)
                geometry_params['l2_diff'].append(l2_diff)
                geometry_params['l12_diff'].append(l12_diff)
                geometry_params['o12_diff'].append(o12_diff)

            geometry_params['geometry'].append(geometry)
            geometry_params['l1'].append(l1)
            geometry_params['l2'].append(l2)
            geometry_params['l12'].append(l12)
            geometry_params['theta'].append(theta)

        except:
            pass

    frequencies = make_numpy_array_from_ragged_list(frequencies)
    intensities = make_numpy_array_from_ragged_list(intensities)
    for k in geometry_params:
        geometry_params[k] = np.array(geometry_params[k])

    print("frequencies.shape:", frequencies.shape)
    print("intensities.shape:", intensities.shape)
    print("geometry_params['geometry'].shape:", geometry_params['geometry'].shape)

    return frequencies, intensities, geometry_params


def make_numpy_array_from_ragged_list(l, fill_value=np.nan):
    maxlen = 0
    for item in l:
        if len(item) > maxlen:
            maxlen = len(item)
    arr = np.full(shape=(len(l), maxlen), fill_value=fill_value)
    for idx, item in enumerate(l):
        arr[idx, :len(item)] = item
    return arr


if __name__ == "__main__":

    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--equilibrium-length', type=float, default=1.160791)

    args = parser.parse_args()

    directory_parent = "/home/eric/Chemistry/calc.sgr/paper_02_CD_SC/scan"
    directory_stretch_1_r = os.path.join(directory_parent, "scan_restricted/scan_stretch_1_6-311++gdp")
    directory_stretch_2_r = os.path.join(directory_parent, "scan_restricted/scan_stretch_2_6-311++gdp")
    directory_bend_r = os.path.join(directory_parent, "scan_unrestricted/scan_bend_6-311++gdp")
    directory_stretch_1_u_s = os.path.join(directory_parent, "scan_unrestricted/scan_stretch_1_6-311++gdp")
    directory_stretch_2_u_s = os.path.join(directory_parent, "scan_unrestricted/scan_stretch_2_6-311++gdp")
    directory_bend_u_s = os.path.join(directory_parent, "scan_unrestricted/scan_bend_6-311++gdp")
    directory_stretch_1_u_t = os.path.join(directory_parent, "scan_unrestricted/scan_stretch_1_6-311++gdp")
    directory_stretch_2_u_t = os.path.join(directory_parent, "scan_unrestricted/scan_stretch_2_6-311++gdp")
    directory_bend_u_t = os.path.join(directory_parent, "scan_unrestricted/scan_bend_6-311++gdp")

    filenames_stretch_1_r = glob(os.path.join(directory_stretch_1_r, "step*"))
    filenames_stretch_2_r = glob(os.path.join(directory_stretch_2_r, "step*"))
    filenames_bend_r = glob(os.path.join(directory_bend_r, "step*"))
    filenames_stretch_1_u_s = glob(os.path.join(directory_stretch_1_u_s, "step*"))
    filenames_stretch_2_u_s = glob(os.path.join(directory_stretch_2_u_s, "step*"))
    filenames_bend_u_s = glob(os.path.join(directory_bend_u_s, "step*"))
    filenames_stretch_1_u_t = glob(os.path.join(directory_stretch_1_u_t, "step*"))
    filenames_stretch_2_u_t = glob(os.path.join(directory_stretch_2_u_t, "step*"))
    filenames_bend_u_t = glob(os.path.join(directory_bend_u_t, "step*"))
    # rename_step_files_eric(filenames_stretch_1_r)
    # rename_step_files_eric(filenames_stretch_2_r)
    # rename_step_files_eric(filenames_bend_r)

    filenames_stretch_1_r = sorted(glob(os.path.join(directory_stretch_1_r, "step_*.out")))
    filenames_stretch_2_r = sorted(glob(os.path.join(directory_stretch_2_r, "step_*.out")))
    filenames_bend_r = sorted(glob(os.path.join(directory_bend_r, "step_*.out")))
    filenames_stretch_1_u_s = sorted(glob(os.path.join(directory_stretch_1_u_s, "step_*.out")))
    filenames_stretch_2_u_s = sorted(glob(os.path.join(directory_stretch_2_u_s, "step_*.out")))
    filenames_bend_u_s = sorted(glob(os.path.join(directory_bend_u_s, "step_*.out")))
    filenames_stretch_1_u_t = sorted(glob(os.path.join(directory_stretch_1_u_t, "step_*.out")))
    filenames_stretch_2_u_t = sorted(glob(os.path.join(directory_stretch_2_u_t, "step_*.out")))
    filenames_bend_u_t = sorted(glob(os.path.join(directory_bend_u_t, "step_*.out")))

    print('stretch_1_r')
    (frequencies_stretch_1_r,
     intensities_stretch_1_r,
     geometry_params_stretch_1_r) = extract(filenames_stretch_1_r,
                                          args.equilibrium_length)
    print('stretch_2_r')
    (frequencies_stretch_2_r,
     intensities_stretch_2_r,
     geometry_params_stretch_2_r) = extract(filenames_stretch_2_r,
                                          args.equilibrium_length)
    print('bend_r')
    (frequencies_bend_r,
     intensities_bend_r,
     geometry_params_bend_r) = extract(filenames_bend_r,
                                       args.equilibrium_length)
    print('stretch_1_u_s')
    (frequencies_stretch_1_u_s,
     intensities_stretch_1_u_s,
     geometry_params_stretch_1_u_s) = extract(filenames_stretch_1_u_s,
                                          args.equilibrium_length)
    print('stretch_2_u_s')
    (frequencies_stretch_2_u_s,
     intensities_stretch_2_u_s,
     geometry_params_stretch_2_u_s) = extract(filenames_stretch_2_u_s,
                                          args.equilibrium_length)
    print('bend_u_s')
    (frequencies_bend_u_s,
     intensities_bend_u_s,
     geometry_params_bend_u_s) = extract(filenames_bend_u_s,
                                       args.equilibrium_length)
    print('stretch_1_u_t')
    (frequencies_stretch_1_u_t,
     intensities_stretch_1_u_t,
     geometry_params_stretch_1_u_t) = extract(filenames_stretch_1_u_t,
                                          args.equilibrium_length)
    print('stretch_2_u_t')
    (frequencies_stretch_2_u_t,
     intensities_stretch_2_u_t,
     geometry_params_stretch_2_u_t) = extract(filenames_stretch_2_u_t,
                                          args.equilibrium_length)
    print('bend_u_t')
    (frequencies_bend_u_t,
     intensities_bend_u_t,
     geometry_params_bend_u_t) = extract(filenames_bend_u_t,
                                       args.equilibrium_length)

    l2_min = 0.75
    l2_max = 2
    l12_min = 2
    l12_max = 4

    # rename the mask variables!
    mask_1_l2 = np.logical_and(geometry_params_stretch_1_r['l2'] >= l2_min,
                               geometry_params_stretch_1_r['l2'] <= l2_max)
    mask_1_l12 = np.logical_and(geometry_params_stretch_1_r['l12'] >= l12_min,
                                geometry_params_stretch_1_r['l12'] <= l12_max)
    mask_2_l12 = np.logical_and(geometry_params_stretch_2_r['l12'] >= l12_min,
                                geometry_params_stretch_2_r['l12'] <= l12_max)

    fig, ax = plt.subplots()

    ax.scatter(geometry_params_stretch_1_r['l2_diff'][mask_1_l2],
               frequencies_stretch_1_r[:, 0][mask_1_l2],
               label='mode 1', linewidth=0.5, marker='o', facecolors='none', edgecolors='blue')
    ax.scatter(geometry_params_stretch_1_r['l2_diff'][mask_1_l2],
               frequencies_stretch_1_r[:, 1][mask_1_l2],
               label='mode 2', linewidth=0.5, marker='o', facecolors='none', edgecolors='green')
    ax.scatter(geometry_params_stretch_1_r['l2_diff'][mask_1_l2],
               frequencies_stretch_1_r[:, 2][mask_1_l2],
               label='mode 3', linewidth=0.5, marker='o', facecolors='none', edgecolors='red')
    ax.scatter(geometry_params_stretch_1_r['l2_diff'][mask_1_l2],
               frequencies_stretch_1_r[:, 3][mask_1_l2],
               label='mode 4', linewidth=0.5, marker='o', facecolors='none', edgecolors='orange')

    # ax.line([equilibrium_length_double, equilibrium_length_double],
    #         [0, ax.get_ylim()[1]],
    #         linestyle='dotted',
    #         color='black')

    ax.autoscale(axis='both', tight=True)

    ax.set_xlabel(r'$l_{2}$ difference from equilibrium ($\AA$)')
    ax.set_ylabel(r'harmonic frequency (cm$^{-1}$)')
    ax.legend(loc='best', fancybox=True, framealpha=0.50, fontsize='x-small')

    fig.savefig('scan_stretch_1_l2.pdf', bbox_inches='tight')

    fig, ax = plt.subplots()

    ax.scatter(geometry_params_stretch_1_r['l12_diff'][mask_1_l12],
               frequencies_stretch_1_r[:, 0][mask_1_l12],
               label='mode 1', linewidth=0.5, marker='o', facecolors='none', edgecolors='blue')
    ax.scatter(geometry_params_stretch_1_r['l12_diff'][mask_1_l12],
               frequencies_stretch_1_r[:, 1][mask_1_l12],
               label='mode 2', linewidth=0.5, marker='o', facecolors='none', edgecolors='green')
    ax.scatter(geometry_params_stretch_1_r['l12_diff'][mask_1_l12],
               frequencies_stretch_1_r[:, 2][mask_1_l12],
               label='mode 3', linewidth=0.5, marker='o', facecolors='none', edgecolors='red')
    ax.scatter(geometry_params_stretch_1_r['l12_diff'][mask_1_l12],
               frequencies_stretch_1_r[:, 3][mask_1_l12],
               label='mode 4', linewidth=0.5, marker='o', facecolors='none', edgecolors='orange')

    ax.autoscale(axis='both', tight=True)

    ax.set_xlabel(r'$l_{1} + l_{2}$ difference from equilibrium ($\AA$)')
    ax.set_ylabel(r'harmonic frequency (cm$^{-1}$)')
    ax.legend(loc='best', fancybox=True, framealpha=0.50, fontsize='x-small')

    fig.savefig('scan_stretch_1_l12.pdf', bbox_inches='tight')

    fig, ax = plt.subplots()

    ax.scatter(geometry_params_stretch_2['l12_diff'][mask_2_l12], frequencies_stretch_2[:, 0][mask_2_l12], label='mode 1', linewidth=0.5, marker='o', facecolors='none', edgecolors='blue')
    ax.scatter(geometry_params_stretch_2['l12_diff'][mask_2_l12], frequencies_stretch_2[:, 1][mask_2_l12], label='mode 2', linewidth=0.5, marker='o', facecolors='none', edgecolors='green')
    ax.scatter(geometry_params_stretch_2['l12_diff'][mask_2_l12], frequencies_stretch_2[:, 2][mask_2_l12], label='mode 3', linewidth=0.5, marker='o', facecolors='none', edgecolors='red')
    ax.scatter(geometry_params_stretch_2['l12_diff'][mask_2_l12], frequencies_stretch_2[:, 3][mask_2_l12], label='mode 4', linewidth=0.5, marker='o', facecolors='none', edgecolors='orange')

    ax.autoscale(axis='both', tight=True)

    ax.set_xlabel(r'$l_{1} + l_{2}$ difference from equilibrium ($\AA$)')
    ax.set_ylabel(r'harmonic frequency (cm$^{-1}$)')
    ax.legend(loc='best', fancybox=True, framealpha=0.50, fontsize='x-small')

    fig.savefig('scan_stretch_2_l12.pdf', bbox_inches='tight')

    fig, ax = plt.subplots()

    ax.scatter(geometry_params_bend['theta'], frequencies_bend[:, 0], label='mode 1', linewidth=0.5, marker='o', facecolors='none', edgecolors='blue')
    ax.scatter(geometry_params_bend['theta'], frequencies_bend[:, 1], label='mode 2', linewidth=0.5, marker='o', facecolors='none', edgecolors='green')
    ax.scatter(geometry_params_bend['theta'], frequencies_bend[:, 2], label='mode 3', linewidth=0.5, marker='o', facecolors='none', edgecolors='red')
    ax.scatter(geometry_params_bend['theta'], frequencies_bend[:, 3], label='mode 4', linewidth=0.5, marker='o', facecolors='none', edgecolors='orange')

    ax.autoscale(axis='both', tight=True)

    ax.set_xlabel(r'$\theta$ (degrees)')
    ax.set_ylabel(r'harmonic frequency (cm$^{-1}$)')
    ax.legend(loc='best', fancybox=True, framealpha=0.50, fontsize='x-small')

    fig.savefig('scan_bend_theta.pdf', bbox_inches='tight')
