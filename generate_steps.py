#!/usr/bin/env python3

from __future__ import print_function

import numpy as np


def template(coord_O1, coord_O2, coord_C=0.0000000000):

    return """$molecule
0 1
C         0.0000000000   {coord_C:>13.10f}    0.0000000000
O         0.0000000000    0.0000000000   {coord_O1:>13.10f}
O         0.0000000000    0.0000000000   {coord_O2:>13.10f}
$end

$rem
SYMMETRY = FALSE
SYM_IGNORE = TRUE
BASIS  =  6-311++G(d,p)
EXCHANGE  =  B3LYP
JOB_TYPE  =  FREQ
SCF_CONVERGENCE  =  8
SCF_ALGORITHM = RCA_DIIS
THRESH  =  12
XC_SMART_GRID = True
XC_GRID = 000099000590
MAX_SCF_CYCLES = 10000
MAX_RCA_CYCLES = 1000
THRESH_RCA_SWITCH = 6
FAST_XC = TRUE
$end
""".format(**locals())


def make_steps(stepsize, stepmin=-0.5, stepmax=10.0):
    steps = np.arange(stepmin, stepmax + stepsize, stepsize)
    # I hate floating point
    for i, step in enumerate(steps):
        if abs(step) < 1.0e-12:
            steps[i] = 0.0
    return steps


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('kind', choices=('stretch_1', 'stretch_2', 'bend', 'none'))
    parser.add_argument('--stepmin', type=float, default=-0.5)
    parser.add_argument('--stepmax', type=float, default=10.0)
    parser.add_argument('--stepsize', type=float, default=0.025)

    args = parser.parse_args()

    steps = make_steps(args.stepsize, stepmin=args.stepmin, stepmax=args.stepmax)
    print(steps)

    start = -1.169160

    if args.kind == 'stretch_1':
        for count, step in enumerate(steps, start=1):
            input_file_contents = template(start, -start + step)
            with open('step_{}.in'.format(count), 'w') as input_file:
                input_file.write(input_file_contents)
    elif args.kind == 'stretch_2':
        for count, step in enumerate(steps, start=1):
            input_file_contents = template(start - step, -start + step)
            with open('step_{}.in'.format(count), 'w') as input_file:
                input_file.write(input_file_contents)
    elif args.kind == 'bend':
        for count, step in enumerate(steps, start=1):
            input_file_contents = template(coord_O1=start,
                                           coord_O2=-start,
                                           coord_C=0.0 + step)
            with open('step_{}.in'.format(count), 'w') as input_file:
                input_file.write(input_file_contents)
    else:
        pass
