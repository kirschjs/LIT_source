import os, sys
import multiprocessing
import subprocess
import shlex
from multiprocessing.pool import ThreadPool

import numpy as np
from bridge import *

anzproc = 2  #int(len(os.sched_getaffinity(0)) / 1)


def cal_rhs_compo(para, procnbr):

    slave_pit = litpath3He + 'tmp_%d' % para[0]
    cmdlu = BINLITpath + 'luise.exe > dump'
    cmdob = BINLITpath + 'obem.exe > dump'
    cmdqu = BINLITpath + 'qual.exe'
    print('%s in %s' % (cmdlu, slave_pit))
    plu = subprocess.Popen(
        shlex.split(cmdlu),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=slave_pit)
    out, err = plu.communicate()
    print('process = %d-1 : luise exits.' % para[0])

    print('%s in %s' % (cmdob, slave_pit))
    pob = subprocess.Popen(
        shlex.split(cmdob),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=slave_pit)
    out, err = pob.communicate()
    print('process = %d-1 : ober exits.' % para[0])

    print('%s in %s' % (cmdqu, slave_pit))
    pqu = subprocess.Popen(
        shlex.split(cmdqu),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=slave_pit)
    out, err = pqu.communicate()
    print('process = %d-1 : qual exits.' % para[0])

    return (out, err)


parameter_set = []

for mM in range(2):
    parameter_set.append([mM])

pool = ThreadPool(anzproc)
results = []
for procnbr in range(len(parameter_set)):
    pars = parameter_set[procnbr]
    results.append(pool.apply_async(cal_rhs_compo, (
        pars,
        procnbr,
    )))
pool.close()
pool.join()