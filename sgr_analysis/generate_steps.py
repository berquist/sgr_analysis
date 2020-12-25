#!/usr/bin/env python3


import math

import numpy as np


ZERO_VEC = np.array([0.0, 0.0, 0.0])


def qchem_template_restricted(coord_O1=ZERO_VEC,
                              coord_O2=ZERO_VEC,
                              coord_C=ZERO_VEC,
                              comment=None):
    # each field should be >13.10f
    assert coord_C.shape == coord_O1.shape == coord_O2.shape == (3,)
    coord_C_X, coord_C_Y, coord_C_Z = coord_C
    coord_O1_X, coord_O1_Y, coord_O1_Z = coord_O1
    coord_O2_X, coord_O2_Y, coord_O2_Z = coord_O2
    return """$comment
{comment}
$end

$molecule
0 1
C         {coord_C_X:>13.10f} {coord_C_Y:>13.10f} {coord_C_Z:>13.10f}
O         {coord_O1_X:>13.10f} {coord_O1_Y:>13.10f} {coord_O1_Z:>13.10f}
O         {coord_O2_X:>13.10f} {coord_O2_Y:>13.10f} {coord_O2_Z:>13.10f}
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


def qchem_template_unrestricted(coord_O1=ZERO_VEC,
                                coord_O2=ZERO_VEC,
                                coord_C=ZERO_VEC,
                                comment=None,
                                mult=1):

    # each field should be >13.10f
    assert coord_C.shape == coord_O1.shape == coord_O2.shape == (3,)
    coord_C_X, coord_C_Y, coord_C_Z = coord_C
    coord_O1_X, coord_O1_Y, coord_O1_Z = coord_O1
    coord_O2_X, coord_O2_Y, coord_O2_Z = coord_O2
    return """$comment
{comment}
$end

$molecule
0 {mult}
C         {coord_C_X:>13.10f} {coord_C_Y:>13.10f} {coord_C_Z:>13.10f}
O         {coord_O1_X:>13.10f} {coord_O1_Y:>13.10f} {coord_O1_Z:>13.10f}
O         {coord_O2_X:>13.10f} {coord_O2_Y:>13.10f} {coord_O2_Z:>13.10f}
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
UNRESTRICTED = TRUE
SCF_GUESS_MIX = TRUE
SCF_GUESS = GWH
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
    parser.add_argument('--r',
                        help="""The bond length to be used. Default value is from B3LYP/6-311++G(d,p).""",
                        type=float,
                        default=-1.160791)
    parser.add_argument('--stepmin', type=float, default=-0.5)
    parser.add_argument('--stepmax', type=float, default=10.0)
    parser.add_argument('--stepsize', type=float, default=0.025)
    parser.add_argument('--mult', type=int, default=1)
    parser.add_argument('--unrestricted', action='store_true')

    args = parser.parse_args()

    if args.unrestricted:
        qchem_template = qchem_template_unrestricted
    else:
        qchem_template = qchem_template_restricted

    steps = make_steps(args.stepsize, stepmin=args.stepmin, stepmax=args.stepmax)
    print(steps)

    r = args.r

    if args.kind == 'stretch_1':
        coord_C = np.array([0.0, 0.0, 0.0])
        coord_O1 = np.array([0.0, 0.0, -r])
        for count, step in enumerate(steps, start=1):
            coord_O2 = np.array([0.0, 0.0, r + step])
            comment = "r: {} step: {} final: {}".format(r, step, r + step)
            input_file_contents = qchem_template(coord_C=coord_C,
                                                 coord_O1=coord_O1,
                                                 coord_O2=coord_O2,
                                                 comment=comment,
                                                 mult=args.mult)
            with open('step_{}.in'.format(count), 'w') as input_file:
                input_file.write(input_file_contents)
    elif args.kind == 'stretch_2':
        coord_C = np.array([0.0, 0.0, 0.0])
        for count, step in enumerate(steps, start=1):
            coord_O1 = np.array([0.0, 0.0, -r - step])
            coord_O2 = np.array([0.0, 0.0, r + step])
            comment = "r: {} step: {} final: {}".format(r, step, r + step)
            input_file_contents = qchem_template(coord_C=coord_C,
                                                 coord_O1=coord_O1,
                                                 coord_O2=coord_O2,
                                                 comment=comment,
                                                 mult=args.mult)
            with open('step_{}.in'.format(count), 'w') as input_file:
                input_file.write(input_file_contents)
    elif args.kind == 'bend':
        coord_C = np.array([0.0, 0.0, 0.0])
        coord_O2 = np.array([r, 0.0, 0.0])
        for count, step in enumerate(steps, start=1):
            theta_deg = 180 - step
            coord_O1 = np.array([r * math.cos(math.radians(theta_deg)),
                                 r * math.sin(math.radians(theta_deg)),
                                 0.0])
            # abs(decimal.Decimal(str(theta_deg)).as_tuple().exponent)
            comment = "r: {} angle: {:.4f} step: {:.4f}".format(r, theta_deg, step)
            input_file_contents = qchem_template(coord_C=coord_C,
                                                 coord_O1=coord_O1,
                                                 coord_O2=coord_O2,
                                                 comment=comment,
                                                 mult=args.mult)
            with open('step_{}.in'.format(count), 'w') as input_file:
                input_file.write(input_file_contents)
    else:
        pass
