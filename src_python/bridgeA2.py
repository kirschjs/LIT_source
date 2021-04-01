import subprocess
import os, fnmatch
import numpy as np
import sympy as sy
# CG(j1, m1, j2, m2, j3, m3)
from sympy.physics.quantum.cg import CG

from parameters_and_constants import *
from rrgm_functions import *
from two_particle_functions import *

import multiprocessing
from smart_diag import *


def cartesian_coord(*arrays):
    grid = np.meshgrid(*arrays)
    coord_list = [entry.ravel() for entry in grid]
    points = np.vstack(coord_list).T
    return points
    #a = np.arange(2)  # fake data
    #print(cartesian_coord(*3 * [a]))


suffix = 'miwchan'
anzproc = 6  #int(len(os.sched_getaffinity(0)) / 1)

home = os.getenv("HOME")

pathbase = home + '/kette_repo/ComptonLIT'

litpathD = pathbase + '/systems/mul_deuteron_' + suffix + '/'

if os.path.isdir(litpathD) != False:
    os.system('rm -rf ' + litpathD)
os.mkdir(litpathD)

deuteronpath = litpathD + 'D/'
if os.path.isdir(deuteronpath) == False:
    os.mkdir(deuteronpath)
v18uixpath = litpathD + 'LITstate/'
respath = litpathD + 'results/'
if os.path.isdir(respath) == False:
    os.mkdir(respath)

BINBDGpath = pathbase + '/source/src_nucl/'
BINLITpath = pathbase + '/source/src_elma_new/'

mpii = '137'
potnn = pathbase + '/data/AV18'
potnnn = pathbase + '/data/urbana9_AK_neu'

# penta+ -------------------------------------------
#1 run <PSI_parallel.py> for boundsatekanal und streukas
#2 1st <A3_lit_par.py>   run
#3 2nd <A3_lit_par.py>   run
cal = [
    'dbg',
    'einzel',
    #'construe_fresh_deuteron',
    #'reduce',
    #'coeff',
    'rhs',
    'rhs_lu-ob-qua',
    'rhs-qual',
    'rhs-end',
    'lhs',
    'lhs_lu-ob-qua',
    #'couple',
]
new_deuteron = True
# convention: bound-state-expanding BVs: (1-8), i.e., 8 states per rw set => nzf0*8
channels = {
    # deuteron
    'np1^+': [
        ['0', ['np_S1']],
        ['2', ['np_S1']],
    ],
    #          [l1l2L,[compatible (iso)spin configurations]]
    '0^-': [
        ['1', ['np_S1']],
    ],
    '1^-': [
        ['1', ['np_S1']],
    ],
    '2^-': [
        ['1', ['np_S1']],
        #['3', ['np_S1']],
    ],
}

streukas = ['0^-', '1^-', '2^-']

#                  realistic    L>0 (only)         deuteron
boundstatekanal = 'np1^+'

J0 = float(boundstatekanal.split('^')[0][-1])

multipolarity = 1

anz_phot_e = 1
phot_e_0 = 0.01  #  enems_e converts to fm^-1, but HERE the value is in MeV
phot_e_d = 1.0  #  delta E

opME_th_low = 10**(-24)
opME_th_up = 10**24

# deuteron/initial-state basis -------------------------------------------
cluster_centers_per_zerl = 3
min_eucl_pair_dist = 0.0001
eps_up = [10.2, 10.01]
eps_low = [0.2, 0.1]

# -- here, I allowed for an enhancement of certain operators, to bind an S-wave triton with v18/uix
costr = ''
for nn in range(1, 14):
    cf = 1.0 if (nn < 28) else 0.0
    cf = 1.0 if (nn == 2) else cf
    costr += '%12.7f' % cf if (nn % 7 != 0) else '%12.7f\n' % cf

#print('costr = ', costr)

# for the cleaner --------------------------------------------------------

maxCoef = 10000
minCoef = 1200
ncycl = 30
maxDiff = 0.001
delPcyc = 1