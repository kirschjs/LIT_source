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

litpath3He = pathbase + '/systems/mul_helion_' + suffix + '/'
if os.path.isdir(litpath3He) == False:
    os.mkdir(litpath3He)

helionpath = litpath3He + 'he3/'
if os.path.isdir(helionpath) == False:
    os.mkdir(helionpath)
v18uixpath = litpath3He + 'LITstate/'
respath = litpath3He + 'results/'
if os.path.isdir(respath) == False:
    os.mkdir(respath)

BINBDGpath = pathbase + '/source/src_nucl/'
BINLITpath = pathbase + '/source/src_elma/'

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
    'construe_fresh_helion',
    'rhs',
    'rhs_lu-ob-qua',
    'rhs-end',
    'lhs',
    'lhs_lu-ob-qua',
    'couple',
    'rhs-qual',
]
new_helion = True
# convention: bound-state-expanding BVs: (1-8), i.e., 8 states per rw set => nzf0*8
channels = {
    # helion
    'npp0.5^+': [
        ['000', ['he_no1', 'he_no6']],
        ['202', ['he_no2']],
        ['022', ['he_no2']],
        ['222', ['he_no2']],
        ['221', ['he_no1']],  #, 'he_no6', 'he_no2']],
        ['220', ['he_no1']],  #, 'he_no6']],
        #['111', ['he_no3']],  #, 'he_no5'
        #['112', ['he_no5']],
    ],
    #          [l1l2L,[compatible (iso)spin configurations]]
    '0.5^-': [
        ['011', ['he_no1', 'he_no2', 'he_no6']],
        ['212', ['he_no2']],
        ['211', ['he_no2']],
        ['101', ['he_no2']],  #, 'he_no4i', 'he_no5']],
        ['121', ['he_no2']],  #, 'he_no4i', 'he_no5']],
        ['122', ['he_no2']],
    ],
    '1.5^-': [
        ['011', ['he_no1', 'he_no2', 'he_no6']],
        ['211', ['he_no1', 'he_no2', 'he_no6']],
        ['212', ['he_no1', 'he_no2', 'he_no6']],
        ['213', ['he_no2']],  #, 'he_no6''he_no1',
        #['101', ['he_no3']],  #, 'he_no4i', 'he_no5'
        #['121', ['he_no3']],  #, 'he_no4i', 'he_no5'
        #['122', ['he_no5']],  #, 'he_no4i'
    ]
}

streukas = ['0.5^-']  #, '1.5^-']

#                  realistic    L>0 (only)         deuteron
boundstatekanal = 'npp0.5^+'

J0 = float(boundstatekanal.split('^')[0][-3:])

multipolarity = 1

anz_phot_e = 4
phot_e_0 = 0.001  #  enems_e converts to fm^-1, but HERE the value is in MeV
phot_e_d = 100.0  #  delta E

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

rw0 = wid_gen(
    add=addw, addtype=addwt, w0=wini0, ths=[1e-5, 2e2, 0.2], sca=scale)
rw0 = sparsify(rw0, min_spacing)

nzf0 = int(np.ceil(len(rw0) / 20.0))

#LIT basis ---------------------------------------------------------------

basisdimLITint = 10
basisdimLITrel = 10
LD = 20

basisdimLIT = 40
wli = 'wd'

if wli == 'wd':
    #scale deuteron
    winiLIT = [ww for ww in 1.1 * rw0 if ((ww < 10) & (ww > 0.01))] + [10.01]

if wli == 'lin':
    #linspace
    w0l, dw = 0.035, 1.4
    winiLIT = np.linspace(
        start=w0l,
        stop=w0l + basisdimLIT * dw,
        num=basisdimLIT,
        endpoint=True,
        dtype=None)
if wli == 'log':
    #logspace
    exp0log, expmaxlog = -1, 1
    winiLIT = np.logspace(
        start=exp0log,
        stop=expmaxlog,
        num=basisdimLIT,
        endpoint=True,
        dtype=None)
if wli == 'geom':
    #geomspace
    wminl, wmaxl = 0.01, 20
    winiLIT = np.geomspace(
        start=wminl, stop=wmaxl, num=basisdimLIT, endpoint=True, dtype=None)
if wli == 'lap':
    #laplace space
    laplace_loc, laplace_scale = .9, .4
    winiLITlaplace = np.sort(
        np.abs(np.random.laplace(laplace_loc, laplace_scale, basisdimLIT)))
    winiLITlaplace = wid_gen(
        add=addw,
        addtype=addwt,
        w0=winiLITlaplace[::-1],
        ths=[1e-5, 2e2, 0.2],
        sca=scale)
    winiLIT = sparsify(winiLITlaplace, min_spacing)

# -- here, I allowed for an enhancement of certain operators, to bind an S-wave triton with v18/uix
costr = ''
for nn in range(1, 30):
    cf = 1.0 if (nn < 28) else 0.0
    cf = 1.0 if (nn == 2) else cf
    costr += '%12.7f' % cf if (nn % 7 != 0) else '%12.7f\n' % cf

#print('costr = ', costr)

# for the cleaner --------------------------------------------------------
streukanalweiten = range(1, len(winiLIT) + 1)

maxCoef = 10000
minCoef = 1200
ncycl = 30
maxDiff = 0.001
delPcyc = 1