import subprocess
import os, fnmatch, copy, struct
import numpy as np
import sympy as sy
# CG(j1, m1, j2, m2, j3, m3)
from sympy.physics.quantum.cg import CG

from parameters_and_constants import *
from rrgm_functions import *
from three_particle_functions import *

import multiprocessing
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
dt = 'float32'

# penta+ -------------------------------------------
#1 run <PSI_parallel.py> for boundsatekanal und streukas
#2 1st <A3_lit_par.py>   run
#3 2nd <A3_lit_par.py>   run
cal = [
    #'construe_fresh_helion',
    #'reduce',
    #'coeff',
    #'diag',
    #'reset',
    'dbg',
    'einzel',
    'rhs_lu-ob-qua',
    'rhs-qual',
    'rhs-end',
    'rhs',
    'rhs-couple',
    'lhs'
]

suffix = 'miwchan-v22'
anzproc = 6  #int(len(os.sched_getaffinity(0)) / 1)

nbrOfBases = 4

home = os.getenv("HOME")

pathbase = home + '/kette_repo/ComptonLIT'

litpath3He = pathbase + '/systems/mul_helion_' + suffix + '/'
basisPath = litpath3He + 'basis_struct/'
v18uixpath = litpath3He + 'LITstate/'
#
#if os.path.isdir(litpath3He) != False:
#    if 'reset' in cal:
#        os.system('rm -rf ' + litpath3He)
#        os.mkdir(litpath3He)
#    else:
#        pass
#else:
#    os.mkdir(litpath3He)
#
#helionpath = litpath3He + 'he3/'
#if os.path.isdir(helionpath) == False:
#    os.mkdir(helionpath)
#v18uixpath = litpath3He + 'LITstate/'
#respath = litpath3He + 'results/'
#if os.path.isdir(respath) == False:
#    os.mkdir(respath)

BINBDGpath = pathbase + '/source/src_nucl/'
BINLITpath = pathbase + '/source/src_elma_new/'
BINLITpathPOL = pathbase + '/source/src_elma_pol/'

# NN: tnni=10   NN+NNN: tnni=11
tnni = 10
parall = -1

mpii = '137'
potnn = pathbase + '/data/AV4.14'  #'/data/AV18'  #'/data/BONN'  #
potnnn = pathbase + '/data/urbana9_AK_neu'

new_helion = True
# convention: bound-state-expanding BVs: (1-8), i.e., 8 states per rw set => nzf0*8
channels = {
    # helion
    'npp0.5^+': [
        #['000', ['he_no1', 'he_no1y', 'he_no6', 'he_no6y']],
        ['000', ['he_no1', 'he_no6']],
        #['220', ['he_no1', 'he_no6']],
        #['221', ['he_no1', 'he_no2', 'he_no6']],
        #['111', ['he_no3', 'he_no3i', 'he_no3ii', 'he_no5', 'he_no5i']],
        #['022', ['he_no2', 'he_no2i']],
        #['202', ['he_no2', 'he_no2i']],
        #['222', ['he_no2']],
        #['112', ['he_no5', 'he_no5i']],
    ],
    #          [l1l2L,[compatible (iso)spin configurations]]
    '0.5^-': [
        #['011', ['he_no1', 'he_no1y', 'he_no6', 'he_no6y']],
        #['101', ['he_no3', 'he_no3y']],
        ['011', ['he_no1', 'he_no6']],
        ['101', ['he_no3']],
        #[
        #    '211',
        #    ['he_no1', 'he_no1i', 'he_no2', 'he_no2i', 'he_no6', 'he_no6i']
        #],
        #['121', ['he_no3', 'he_no3i', 'he_no3ii', 'he_no5', 'he_no5i']],
        #['212', ['he_no2', 'he_no2i']],
        #['122', ['he_no5', 'he_no5i']],
    ],
    '1.5^-': [
        ['011', ['he_no1', 'he_no1y', 'he_no2', 'he_no6']],
        ['011', ['he_no6', 'he_no6', 'he_no6', 'he_no6']],
        ['211', ['he_no1', 'he_no1y', 'he_no2', 'he_no6']],
        ['212', ['he_no2']],  #,'he_no1',  'he_no6']],
        ['213', ['he_no2']],  #,'he_no6',  'he_no1']],
        #['101', ['he_no3']],  #, 'he_no4i', 'he_no5'
        #['121', ['he_no3']],  #, 'he_no4i', 'he_no5'
        #['122', ['he_no5']],  #, 'he_no4i'
    ]
}

streukas = ['0.5^-']  #,'1.5^-']

#                  realistic    L>0 (only)         deuteron
boundstatekanal = 'npp0.5^+'

J0 = float(boundstatekanal.split('^')[0][-3:])

multipolarity = 1

anz_phot_e = 1
phot_e_0 = 0.1  #  enems_e converts to fm^-1, but HERE the value is in MeV
phot_e_d = 1.0  #  delta E

opME_th_low = 10**(-24)
opME_th_up = 10**24

# deuteron/initial-state basis -------------------------------------------
cluster_centers_per_zerl = 3
min_eucl_pair_dist = 0.0001
eps_up = [10.2, 10.01]
eps_low = [0.2, 0.1]

wini0 = w120  #wLAPLACE[::-1]

addw = 1
addwt = 'middle'
scale = 1.
min_spacing = 0.02
min_spacing_to_LITWs = 0.001

rw0 = wid_gen(add=addw,
              addtype=addwt,
              w0=wini0,
              ths=[1e-5, 2e2, 0.2],
              sca=scale)
rw0 = sparsify(rw0, min_spacing)

nzf0 = int(np.ceil(len(rw0) / 20.0))

# basis ------------------------------------------------------------------

bvma = 12
rwma = 45

# -- here, I allowed for an enhancement of certain operators, to bind an S-wave triton with v18/uix
costr = ''

zop = 31 if tnni == 11 else 14

for nn in range(1, zop):
    cf = 1.0 if (nn < 28) else 0.0
    cf = 0.0 if (nn == 1) else cf
    costr += '%12.7f' % cf if (nn % 7 != 0) else '%12.7f\n' % cf

#print('costr = ', costr)
#print(BINBDGpath)
# for the cleaner --------------------------------------------------------

maxCoef = 10000
minCoef = 1200
ncycl = 30
maxDiff = 0.001
delPcyc = 1