import pickle
import numpy as np
import scipy.stats as sps
from collections import OrderedDict

# from parse_outputs_snapshot_method_dependence import \
#     (methods, basis_sets)
methods = OrderedDict([
    ('blyp', 'BLYP'),
    ('tpss', 'TPSS'),
    ('b3lyp', 'B3LYP'),
    ('wb97x-d', r'$\omega$B97X-D'),
    ('hf', 'HF'),
    ('ri-mp2', 'RI-MP2'),
])

basis_sets = OrderedDict([
    ('6-31gdp', '6-31G(d,p)'),
    ('cc-pvtz', 'cc-pVTZ'),
])

def fit_line(x, y):
    """Return slope, intercept, R-squared of best fit line."""
    slope, intercept, r, p, stderr = sps.linregress(x, y)
    return slope, intercept, r**2


if __name__ == '__main__':

    # Read in the pickle files that contain all the raw data.
    with open('frequencies.pypickle', 'rb') as picklefile:
        frequencies_d = pickle.load(picklefile)
    with open('snapnums_frequencies.pypickle', 'rb') as picklefile:
        snapnums_d = pickle.load(picklefile)

    n_qm = 0
    n_mm = 256

    # Sort by increasing RI-MP2/cc-pVTZ frequency.
    frequencies_benchmark = frequencies_d['ri-mp2']['cc-pvtz'][n_qm][n_mm]
    idxsort = np.argsort(frequencies_benchmark)
    frequencies_benchmark = np.array(frequencies_benchmark)[idxsort]

    for basis_set in basis_sets:
        for idx, method in enumerate(methods):

            frequencies = np.array(frequencies_d[method][basis_set][n_qm][n_mm])[idxsort]
            print(method, basis_set)
            print(fit_line(frequencies_benchmark, frequencies))
