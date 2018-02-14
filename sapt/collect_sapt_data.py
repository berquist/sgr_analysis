"""collect_sapt_data.py: Parse a set of SAPT0 output files for their
interaction energy components, writing to a CSV file."""

from __future__ import division
from __future__ import print_function

import os.path

from helpers import (read_psi4_sapt0, sapt_bases)


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('filename', nargs='+')
    parser.add_argument('--averages', action='store_true')
    parser.add_argument('--plot', action='store_true')

    args = parser.parse_args()

    import csv
    import xlsxwriter

    csvfh = open('data_sapt.csv', 'w')
    csvwriter = csv.writer(csvfh)

    workbook = xlsxwriter.Workbook('data_sapt.xlsx')
    worksheet = workbook.add_worksheet()
    rowcounter = 0

    # This isn't exactly the border one would get if you selected a
    # range of values and changed the stype to "Output", but it's
    # close enough.
    style_output = workbook.add_format({'bg_color': '#F2F2F2', 'font_color': '#3F3F3F', 'bold': True, 'border': 1})
    style_green = workbook.add_format({'bg_color': '#C6EFCE', 'font_color': '#006100'})
    style_yellow = workbook.add_format({'bg_color': '#FFEB9C', 'font_color': '#9C5700'})
    style_red = workbook.add_format({'bg_color': '#FFC7CE', 'font_color': '#9C0006'})

    # csvwriter.writerow(['all values in kcal/mol'])
    header = [
        'filename',
        'snapnum',
        'weight',
        'sapt_basis',
        'electrostatics',
        'exchange',
        '"frz"',
        'induction',
        'exchange_induction',
        'induction_delta_hf',
        '"pol"',
        '"frz" + "pol"',
        'dispersion',
        'exchange_dispersion',
        '"disp"',
        'total',
        # '"tot"',
        'ct',
    ]
    csvwriter.writerow(header)
    worksheet.write_row(rowcounter, 0, header)
    rowcounter += 1

    idx_total = header.index('total')

    rows_monomer = []
    rows_dimer = []

    for filename in sorted(args.filename):

        filename = os.path.abspath(filename)

        print(filename)

        sapt_data = read_psi4_sapt0(filename)

        for sapt_basis in sapt_bases:

            if sapt_basis in sapt_data:

                # Approximate the ALMO frozen density term
                # (electrostatics + Pauli repulsion).
                frz = sapt_data[sapt_basis]['el'] + \
                      sapt_data[sapt_basis]['exch']

                # Approximate the ALMO polarization term as SAPT
                # induction "to all orders".
                pol = sapt_data[sapt_basis]['ind'] + \
                      sapt_data[sapt_basis]['exch-ind'] + \
                      sapt_data[sapt_basis]['ind_HO']

                # Combine the pure dispersion term with
                # exchange-dispersion.
                disp = sapt_data[sapt_basis]['disp'] + \
                       sapt_data[sapt_basis]['exch-disp']

                # Define the total SAPT interaction energy as
                # everything (frz + pol + disp) *except* the CT term.
                # Unnecessary since we already have the total column.
                # tot = frz + pol + disp

                row = [
                    filename,
                    sapt_data['snapnum'],
                    sapt_data['weight'],
                    sapt_basis,
                    sapt_data[sapt_basis]['el'],
                    sapt_data[sapt_basis]['exch'],
                    frz,
                    sapt_data[sapt_basis]['ind'],
                    sapt_data[sapt_basis]['exch-ind'],
                    sapt_data[sapt_basis]['ind_HO'],
                    pol,
                    frz + pol,
                    sapt_data[sapt_basis]['disp'],
                    sapt_data[sapt_basis]['exch-disp'],
                    disp,
                    sapt_data[sapt_basis]['total'],
                    # tot,
                ]

                if 'ct' in sapt_data:
                    row.append(sapt_data['ct'])
                else:
                    row.append('')
                if sapt_basis == 'monomer':
                    rows_monomer.append(row)
                if sapt_basis == 'dimer':
                    rows_dimer.append(row)

    rows = rows_monomer + rows_dimer
    csvwriter.writerows(rows)
    for row in rows:
        worksheet.write_row(rowcounter, 0, row)
        rowcounter += 1

    if args.averages:
        start = 4
        import numpy as np
        for sapt_basis_rows in (rows_monomer, rows_dimer):
            if sapt_basis_rows:
                sapt_basis = sapt_basis_rows[0][3]
                snapnums = [row[0] for row in sapt_basis_rows]
                weights = [row[2] for row in sapt_basis_rows]
                vals = [row[start:] for row in sapt_basis_rows]
                average_unweighted = np.average(vals, axis=0, weights=None)
                average_weighted = np.average(vals, axis=0, weights=weights)
                # Assumes we have CT data?
                tot_average_unweighted = average_unweighted[-2]
                tot_average_weighted = average_weighted[-2]
                frac_unweighted = 100 * average_unweighted / tot_average_unweighted
                frac_weighted = 100 * average_weighted / tot_average_weighted
                row_average_unweighted = ['average (unweighted) [{}]'.format(sapt_basis), '', '', ''] + average_unweighted.tolist()
                row_average_weighted = ['average (weighted) [{}]'.format(sapt_basis), '', '', ''] + average_weighted.tolist()
                row_frac_unweighted = ['100 * fraction of total (unweighted) [{}]'.format(sapt_basis), '', '', ''] + frac_unweighted.tolist()
                row_frac_weighted = ['100 * fraction of total (weighted) [{}]'.format(sapt_basis), '', '', ''] + frac_weighted.tolist()
                csvwriter.writerow(row_average_unweighted)
                csvwriter.writerow(row_average_weighted)
                csvwriter.writerow(row_frac_unweighted)
                csvwriter.writerow(row_frac_weighted)
                worksheet.write_row(rowcounter, 0, row_average_unweighted)
                rowcounter += 1
                worksheet.write_row(rowcounter, 0, row_average_weighted, style_output)
                rowcounter += 1
                worksheet.write_row(rowcounter, 0, row_frac_unweighted)
                rowcounter += 1
                worksheet.write_row(rowcounter, 0, row_frac_weighted)
                rowcounter += 1

        # if args.plot:

        #     from collections import OrderedDict
        #     labels = [
        #         ('electrostatics', 'electrostatics ($E_{elst}^{(10)}$)'),
        #         ('exchange', 'exchange ($E_{exch}^{(10)}$)'),
        #         ('induction', 'induction ($E_{ind,resp}^{(20)}$)'),
        #         ('exchange_induction', 'exchange-induction ($E_{exch-ind,resp}^{(20)}$)'),
        #         ('induction_delta_hf', 'delta HF ($\delta_{HF}^{(2)}$)'),
        #         ('dispersion', 'dispersion ($E_{disp}^{(20)}$)'),
        #         ('exchange_dispersion', 'exchange-dispersion ($E_{exch-disp}^{(20)}$)'),
        #         ('total', 'total ($E_{SAPT0}$)'),
        #     ]
        #     if has_ct:
        #         labels.append(('ct', 'CT'))
        #     labels = OrderedDict(labels)

        #     import matplotlib as mpl
        #     mpl.use('Agg')
        #     import matplotlib.pyplot as plt

        #     nxticks = vals.shape[0]
        #     xticks = list(range(nxticks))
        #     xticklabels = snapnums

        #     # cmap = cm.get_cmap('viridis')

        #     fig, ax = plt.subplots()

        #     # ax.plot(vals[:, 0], label='electrostatics', marker='o')
        #     # ax.plot(vals[:, 1], label='exchange', marker='o')
        #     # ax.plot(vals[:, 2], label='induction', marker='o')
        #     # ax.plot(vals[:, 3], label='exchange-induction', marker='o')
        #     # ax.plot(vals[:, 4], label='delta HF', marker='o')
        #     # ax.plot(vals[:, 5], label='dispersion', marker='o')
        #     # ax.plot(vals[:, 6], label='exchange-dispersion', marker='o')
        #     # ax.plot(vals[:, 7], label='total', marker='s', linestyle='--')
        #     # if has_ct:
        #     #     ax.plot(vals[:, 8], label='CT', marker='s', linestyle='--')
        #     for idx, label in enumerate(labels):
        #         if label == 'total' or label == 'ct':
        #             ax.plot(vals[:, idx],
        #                     label=labels[label],
        #                     linestyle='--',
        #                     marker='s')
        #         else:
        #             ax.plot(vals[:, idx],
        #                     label=labels[label],
        #                     marker='o')

        #     ax.legend(loc='best', fancybox=True, framealpha=0.50, fontsize='x-small')

        #     ax.tick_params(direction='out')

        #     ax.set_xticks(xticks)
        #     ax.set_xticklabels(xticklabels, fontsize='x-small')

        #     ax.set_xlabel('snapshot #')
        #     ax.set_ylabel('energy (kcal/mol)')

        #     # fig.savefig('sapt_data.pdf', bbox_inches='tight')

    csvfh.close()
    workbook.close()
