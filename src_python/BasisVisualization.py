from bridgeA3 import *
from three_particle_functions import *
from triton_width_gen import *
import operator

import matplotlib.pyplot as plt


def visbas(basispath, widthpath, distances=[], Jstrstr='0.5', btype='0.5^-'):

    curd = os.getcwd()
    os.chdir(basispath)

    litbas_full = np.loadtxt(basispath + 'LITbas_full_J%s_%s.dat' %
                             (Jstrstr, btype)).astype(int)
    try:
        litbas_red = np.loadtxt(basispath + 'LITbas_red_J%s_%s.dat' %
                                (Jstrstr, btype)).astype(int)
    except:
        litbas_red = litbas_full

    basvs = {}
    for bv in litbas_red:
        try:
            basvs[bv[0]].append(bv[1])
        except:
            basvs[bv[0]] = [bv[1]]
    sbas = sorted(basvs.items(), key=operator.itemgetter(0))

    n3_inen_bdg(
        sbas, float(Jstrstr), costr, fn='INEN', pari=0, nzop=31, tni=11)

    intwLIT = [
        np.array(ln.split(';')).astype(float).tolist()
        for ln in open(widthpath + 'intw3heLIT_J%s_%s.dat' % (Jstrstr, btype))
    ]

    relwLIT = [
        np.array(ln.split(';')).astype(float).tolist()
        for ln in open(widthpath + 'relw3heLIT_J%s_%s.dat' % (Jstrstr, btype))
    ]

    iws_full = []
    rws_full = []

    iws_red = []
    rws_red = []

    numered_widths = []
    #print(' BV REL          wi          wr')

    for bv in litbas_full:
        for fr in range(len(intwLIT)):
            if bv[0] <= sum([len(ws) for ws in intwLIT[:fr + 1]]):
                iw = bv[0] - sum([len(ws) for ws in intwLIT[:fr]])
                iws_full.append(float(intwLIT[fr][iw - 1]))
                rws_full.append(float(relwLIT[fr][bv[1] - 1]))
                break

    for bv in litbas_red:
        for fr in range(len(intwLIT)):
            if bv[0] <= sum([len(ws) for ws in intwLIT[:fr + 1]]):
                iw = bv[0] - sum([len(ws) for ws in intwLIT[:fr]])
                iws_red.append(float(intwLIT[fr][iw - 1]))
                rws_red.append(float(relwLIT[fr][bv[1] - 1]))

                numered_widths.append([bv[0], bv[1], iws_red[-1], rws_red[-1]])
                break

    combos = [[numered_widths[bv1], numered_widths[bv2]]
              for bv1 in range(len(numered_widths))
              for bv2 in range(1 + bv1, len(numered_widths))]

    distances = []
    for pair in combos:
        normdist = np.linalg.norm(
            np.array(pair[0][2:]) - np.array(pair[1][2:]))
        distances.append(normdist)

    numered_widths = np.array(numered_widths)

    numered_widths = numered_widths[np.lexsort(([
        numered_widths[:, i] for i in range(numered_widths.shape[1] - 2,
                                            numered_widths.shape[1] - 1, +1)
    ]))]

    with open(v18uixpath + 'LITbas_red_J%s_doc.dat' % Jstrstr, 'w') as f:
        np.savetxt(f, numered_widths, fmt='%5d  %5d  %12.8f  %12.8f')
    f.close()

    fig = plt.figure(figsize=[12, 6])
    #fig.subplots_adjust(hspace=0.4, wspace=0.4)
    ax = plt.subplot(121)

    ax.plot(
        iws_full,
        rws_full,
        'o',
        alpha=0.3,
        color='blue',
        label=r'full basis set',
        marker='.',
        markersize=2)
    ax.plot(
        iws_red,
        rws_red,
        'o',
        alpha=0.8,
        marker='.',
        markeredgecolor='red',
        markersize=4,
        color='k',
        label=r'reduced basis set')

    plt.xlabel(r'$\gamma(\rho_1)\;\;\; [fm^{-2}]$')
    plt.ylabel(r'$\gamma(\rho_2)\;\;\; [fm^{-2}]$')

    plt.legend(loc='best')

    plt.title(
        r'$\gamma_1<\gamma_2\Rightarrow$ 1 broader than 2   (J=%s)' % Jstrstr)

    clis = np.linspace(np.min(distances), np.max(distances), 10)

    hist, bin_edges = np.histogram(distances, bins=clis)
    ax = plt.subplot(122)
    ax.set_xlabel(r'$|.|$')
    ax.set_ylabel(r'')
    ax.bar(bin_edges[:-1], hist, color='black', alpha=0.9, width=0.001)

    fig.savefig(litpath3He + 'results/WidthXY_%s_J%s.pdf' % (suffix, Jstrstr))

    print('Width grid plotted and saved to %s' %
          (litpath3He + 'results/WidthXY_J%s.pdf' % Jstrstr))

    os.chdir(curd)


#Jstreu = float(streukas[-1].split('^')[0])
#Jstreustring = '%s' % str(Jstreu)[:3]
##
#os.chdir(litpath3He + 'lit_bas_rhs/')
##
#os.system('cp QUAOUT_J%s QUAOUT' % Jstreustring)
#os.system('cp DRQUAOUT_J%s DRQUAOUT' % Jstreustring)
##
#visbas(
#    basispath=litpath3He + 'lit_bas_rhs/',
#    widthpath=litpath3He,
#    exepath=BINBDGpath,
#    Jstrstr=Jstreustring)
#
#os.system(BINBDGpath + 'DR2END_AK.exe && grep -A 3 \'EIGENWER\' OUTPUT')