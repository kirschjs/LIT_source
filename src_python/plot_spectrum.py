import os
import math
import numpy as np
import matplotlib.pyplot as plt
from bridge import *


def plotHspec(Jj='tmp', inf='', outf=''):

    home = litpath3He + 'lit_bas_lhs/'
    if inf == '':
        inf = home + 'OUTPUT'

    out = [line for line in open(inf)]

    ham = []
    for nj in range(1, len(out)):
        if (out[nj].strip() == 'EIGENWERTE DES HAMILTONOPERATORS'):
            nl = 3
            while out[nj + nl] != '\n':
                ham.append(out[nj + nl].strip().split())
                nl += 1

    ham = np.array(np.concatenate(ham)).astype(float)
    #print(ham, int(np.max(ham)))
    dimi = np.size(ham)

    clis = np.arange(0, int(np.max(ham) + 1), 0.25)
    #print(clis)
    hist, bin_edges = np.histogram(ham, bins=clis)

    fig = plt.figure(figsize=[12, 6])
    #fig.subplots_adjust(hspace=0.4, wspace=0.4)
    ax1 = plt.subplot(121)
    ax1.set_title(
        r'Eigenvalue density in the LIT basis of the V18+UIX Hamilationian')
    ax1.set_xlabel(r'$E_n$ [MeV]')
    ax1.set_ylabel(r'')
    ax1.plot(ham, np.ones(dimi), '.', color='black', label=r'$|d|=%d$' % dimi)
    leg = ax1.legend(loc='best')

    ax1 = plt.subplot(122)
    ax1.set_xlabel(r'$E_n$ [MeV]')
    ax1.set_ylabel(r'')
    ax1.bar(bin_edges[:-1], hist, width=0.25, color='black', alpha=0.9)

    if outf == '':
        outf = home + 'h_spect_%s_%s.pdf' % (suffix, Jj)

    fig.savefig(outf)

    print('Eigenvalue spectrum plotted in %s' % outf)