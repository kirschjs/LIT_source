import subprocess
import multiprocessing
import os, fnmatch, copy, struct, time
import numpy as np
import sympy as sy
# CG(j1, m1, j2, m2, j3, m3)
from sympy.physics.quantum.cg import CG
from datetime import datetime
import time

from parameters_and_constants import *
from rrgm_functions import *
from three_particle_functions import *

from smart_diag import *


def cartesian_coord(*arrays):
    grid = np.meshgrid(*arrays)
    coord_list = [entry.ravel() for entry in grid]
    points = np.vstack(coord_list).T
    return points
    #a = np.arange(2)  # fake data
    #print(cartesian_coord(*3 * [a]))


# float16 : Half   precision float: sign bit,  5 bits exponent, 10 bits mantissa
# float32 : Single precision float: sign bit,  8 bits exponent, 23 bits mantissa
# float64 : Double precision float: sign bit, 11 bits exponent, 52 bits mantissa
#
NEWLINE_SIZE_IN_BYTES = -1
dt = 'float64'

#1 run <PSI_parallel.py> for boundsatekanal und streukas
#2 1st <A3_lit_par.py>   run
#3 2nd <A3_lit_par.py>   run
cal = [
    'dbg',
    'einzel',
    'rhs_lu-ob-qua',
    'rhs-qual',
    'rhs-end',
    'rhs',
    'rhs-couple',
]

DC = True if time.tzname[0] == 'EST' else False
MaxProc = int(len(os.sched_getaffinity(0)) / 2)

orig_dir = os.getcwd()

MPIRUN = subprocess.run(["which", "mpirun"],
                        capture_output=True).stdout.strip().decode()

bkpdir = os.getenv("HOME") + '/compton_tmp' if os.path.isdir(
    '/scratch') == False else '/scratch/compton_tmp'
if os.path.isdir(bkpdir) == False:
    os.mkdir(bkpdir)

pathbase = '/tmp'
suffi = '/mul_helion/'
litpath3He = pathbase + suffi
respath = litpath3He + 'results/'
if os.path.isdir(litpath3He) == False:
    os.mkdir(litpath3He)
    os.mkdir(respath)
    with open(respath + 'dtype.dat', 'w') as outf:
        outf.write(dt)

helionpath = litpath3He + 'he3/'

BINBDGpath = orig_dir + '/../../src_nucl/'
BINLITpath = orig_dir + '/../../src_elma_pol/'

# NN: tnni=10   NN+NNN: tnni=11
tnni = 10
parall = -1

nnStr = 'AV18'  #'BONN'  #AV4.14'
potnn = orig_dir + '/../../data/%s' % nnStr
nnnStr = 'urbana9_AK_neu'
potnnn = orig_dir + '/../../data/%s' % nnnStr

# convention: bound-state-expanding BVs: (1-8), i.e., 8 states per rw set => nzf0*8
channels = {
    # helion
    'npp0.5^+': [
        ['000', ['he_no1', 'he_no6']],  # 1,2
        ['022', ['he_no2']],  # 3
        ['202', ['he_no2']],  # 4
        ['111', ['he_no3', 'he_no5']],  # 5,6
        ['112', ['he_no5']],  # 7
        ['220', ['he_no1', 'he_no6']],  # 8,9
        ['221', ['he_no1', 'he_no2', 'he_no6']],  # 10,11,12
        ['222', ['he_no2']],  # 13
    ],
    #          [l1l2L,[compatible (iso)spin configurations]]
    '0.5^-': [
        ['011', ['he_no1', 'he_no6']],
        ['101', ['he_no3']],
        ['211', ['he_no2', 'he_no1', 'he_no6']],
        ['212', ['he_no2']],
        ['121', ['he_no3', 'he_no5']],
        ['122', ['he_no5']],
    ],
    '1.5^-': [
        ['011', ['he_no1', 'he_no2', 'he_no6']],
        ['101', ['he_no3']],
        ['211', ['he_no1', 'he_no2', 'he_no6']],
        ['212', ['he_no2']],
        ['121', ['he_no3', 'he_no5']],
        ['122', ['he_no3', 'he_no5']],
        ['213', ['he_no2']],
    ]
}

streukas = ['0.5^-', '1.5^-']

#                  realistic    L>0 (only)         deuteron
boundstatekanal = 'npp0.5^+'

J0 = float(boundstatekanal.split('^')[0][-3:])

multipolarity = 1

anz_phot_e = 1
phot_e_0 = 0.1  #  enems_e converts to fm^-1, but HERE the value is in MeV
phot_e_d = 1.0  #  delta E

# basis ------------------------------------------------------------------

# maximal number of basis vectors per calculation block (numerical parameter)
bvma = 12

# maximal number of radial Gauss widths for the expansion of the coordinate
# between the two fragments which are associated with a basis vector (numerical
# paramter)
rwma = 45