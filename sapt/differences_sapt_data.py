import csv

import numpy as np


HEADER = [
    'filename',
    'snapnum',
    'weight',
    'sapt_basis',
    'electrostatics',
    'exchange',
    'induction',
    'exchange_induction',
    'induction_delta_hf',
    'dispersion',
    'exchange_dispersion',
    'total',
    'ct',
]
HEADER2 = [
    'elstat',
    'exch',
    'ind',
    'exch-ind',
    'ind-del-HF',
    'disp',
    'exch-disp',
    'tot',
    'ct',
]

START = 4

def read_sapt_csv(csvfilename):
    with open(csvfilename) as csvfile:
        reader = csv.reader(csvfile)
        # row = next(reader)
        # assert row == ['all values in kcal/mol']
        row = next(reader)
        assert row == HEADER
        # for now, just look at the unweighted and weighted averages
        data = []
        for row in reader:
            if 'average (unweighted) [dimer]' in row[0]:
                assert len(row) == len(HEADER)
                data.append([float(x) for x in row[START:]])
            elif 'average (weighted) [dimer]' in row[0]:
                assert len(row) == len(HEADER)
                data.append([float(x) for x in row[START:]])
            else:
                pass
        return data


if __name__ == '__main__':

    filename_1 = './631gdp/ct/sapt_data_all.csv'
    filename_2 = './jun-cc-pvtz/ct/sapt_data_all.csv'

    data_1 = np.array(read_sapt_csv(filename_1))
    data_2 = np.array(read_sapt_csv(filename_2))

    th = ('{:>13} ' * 9).format
    t = ('{:>13f} ' * 9).format
    print(th(*HEADER2))
    print('jun-cc-pVTZ [weighted]')
    print(t(*data_2[1]))
    print('jun-cc-pVTZ [unweighted]')
    print(t(*data_2[0]))
    print('6-31G(d,p) [weighted]')
    print(t(*data_1[1]))
    print('6-31G(d,p) [unweighted]')
    print(t(*data_1[0]))
    print('[jun-cc-pVTZ] - [6-31G(d,p)] [weighted]')
    print(t(*(data_2[1] - data_1[1])))
    print('[jun-cc-pVTZ] - [6-31G(d,p)] [unweighted]')
    print(t(*(data_2[0] - data_1[0])))
