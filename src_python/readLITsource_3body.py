import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import shutil
import re

from bridge import *
# CG(j1, m1, j2, m2, j3, m3)
from sympy.physics.quantum.cg import CG

MeVfm = 197.3161329


def read_norm(streukanal, basisSET=''):
    # collects end OUTPUT files and stores them
    # *without* modifications

    # find mL, mJl s.t. (L,mL;Jr,mJl-mL|Jl,mJl) != 0 -------------------------
    # ecce: Jr = Jdeuteron = 1
    # mM[0] = m(L) ; mM[1] = m(Jlit)
    # collects end OUTPUT files and stores them
    # *without* modifications

    # find mL, mJl s.t. (L,mL;Jr,mJl-mL|Jl,mJl) != 0 -------------------------
    #
    # mM[0] = m(L) ; mM[1] = m(Jlit)
    Jstreu = float(streukanal.split('^')[0])

    Jstrstr = '%s' % str(Jstreu)[:3]

    mLmJl, mLrange, mJlrange = non_zero_couplings(multipolarity, J0, Jstreu)

    sourceRHS = {}
    basdim = len(basisSET)

    # as norm ME's must be photon-energy independent, only one line is read
    photen2write = 0

    for mM in mLmJl:
        outp = litpath3He + 'LIT_SOURCE_%s_mLmJ_%1.1f_%1.1f_NORM' % (
            streukanal, float(mM[0]), float(mM[1]))
        # output string: ME bv1\n ME bv2\n...
        outs = ''

        for bv in basisSET:
            # mM[0] = m(L) ; mM[1] = m(Jlit)

            instream = [
                line
                for line in open('inhomo%d-%d_J%f_mJ%f-mL%d.log' % (bv[0],
                                                                    bv[1],
                                                                    Jstreu,
                                                                    mM[1],
                                                                    mM[0]))
            ][photen2write + 1]

            JDEUT2 = int(instream.split()[3])
            JLIT2 = int(instream.split()[1])
            mJLIT2 = int(instream.split()[4])
            MUL2 = int(instream.split()[2])
            mMUL2 = int(instream.split()[5])
            if ((MUL2 != 0) | (mMUL2 != 0)):
                print('ECCE: MUL of Norm != 0.')
                exit()

            opME = np.array([float(instream.split()[6])])

            sourceRHS[('%d-%d' % (bv[0], bv[1]), '%d' % JLIT2, '%d' % mJLIT2,
                       '%d' % MUL2, '%d' % mMUL2)] = opME

            outs += '%12.12E  ' % opME

            outs += '\n'

        if 'purge' in cal:
            if os.path.isfile(outp):
                print('removing previous <LIT_SOURCE>')
                os.system('rm ' + outp)
        print('writing source column to <%s>' % outp)
        with open(outp, 'w') as outfile:
            #outfile.seek(0)
            outfile.write(outs)
        #print('I read r,2*(Jlit,mJlit,L,mL):', streukanalweite,
        #      JLIT2, mJLIT2, MUL2, mMUL2)

    return sourceRHS


def read_uncoupled_source(RHSofBVsk,
                          photEn,
                          streukanal,
                          basisSET='',
                          firstbv=0):
    # collects end OUTPUT files and stores them
    # *without* modifications

    # find mL, mJl s.t. (L,mL;Jr,mJl-mL|Jl,mJl) != 0 -------------------------
    #
    # mM[0] = m(L) ; mM[1] = m(Jlit)
    Jstreu = float(streukanal.split('^')[0])

    Jstrstr = '%s' % str(Jstreu)[:3]

    mLmJl, mLrange, mJlrange = non_zero_couplings(multipolarity, J0, Jstreu)

    signf_overlap_bv = []
    #sourceRHS = {}
    basdim = len(basisSET)

    bv_number = firstbv

    print(basisSET)

    for bv in basisSET:
        for mM in mLmJl:
            # mM[0] = m(L) ; mM[1] = m(Jlit)
            instream = [
                line for line in open('endlit%d-%d_J%3.1f_mJ%3.1f-mL%d.log' % (
                    bv[0], bv[1], Jstreu, mM[1], mM[0]))
            ]
            #print('endlit%d-%d_J%3.1f_mJ%3.1f-mL%d.log' %
            #      (bv[0], bv[1], Jstreu, mM[1], mM[0]))
            for ln in range(len(instream)):
                if re.search('EK2', instream[ln]):
                    JDEUT2 = int(instream[ln + 2].split()[4])
                    JLIT2 = int(instream[ln + 2].split()[2])
                    mJLIT2 = int(instream[ln + 2].split()[5])
                    MUL2 = int(instream[ln + 2].split()[3])
                    mMUL2 = int(instream[ln + 2].split()[6])
                    photon_energy = MeVfm * np.array([
                        float(instream[ln + 2 + 3 * en].split()[1])
                        for en in range(anz_phot_e)
                    ])
                    opME = np.array([
                        float(instream[ln + 2 + 3 * en].split()[7])
                        for en in range(anz_phot_e)
                    ])
                    if (((abs(opME) > opME_th_low).all()) &
                        ((abs(opME) < opME_th_up).all())):
                        signf_overlap_bv.append(bv)

                    RHSofBVsk[('%d' % (bv_number), '%d' % JLIT2, '%d' % mJLIT2,
                               '%d' % MUL2, '%d' % mMUL2)] = opME

                    #print('I read 2*(Jlit,mJlit,L,mL):', JLIT2, mJLIT2, MUL2,
                    #      mMUL2)
        bv_number += 1

    if signf_overlap_bv != []:
        with open('LITbas_red_J%s_sig.dat' % Jstrstr, 'w') as f:
            np.savetxt(f, np.unique(signf_overlap_bv, axis=0), fmt='%5d  %5d')
        f.close()
        photEn = photon_energy
    else:
        print(
            '------ ECCE: no basis vector with significant overlap with the perturbed state! EXIT '
        )
        #exit()

    return  #sourceRHS, photon_energy


# return S[ r,Jlit,m(Jlit),L, [k1,k2,...,kn] ]
def couple_source(RHSofmJsk, streukanal, sourceRHS, basisSET='', firstbv=0):

    #coupledSOURCE = {}
    basdim = len(basisSET)

    # mM[0] = m(L) ; mM[1] = m(Jlit)
    Jstreu = float(streukanal.split('^')[0])
    mLmJl, mLrange, mJlrange = non_zero_couplings(multipolarity, J0, Jstreu)

    bv_number = firstbv
    for bv in basisSET:
        for mJ in mJlrange:
            RHSofmJsk[('%d' % (bv_number), '%d' % (2 * Jstreu),
                       '%d' % (2 * mJ), '%d' % (2 * multipolarity))] = 0.0
        bv_number += 1

    bv_number = firstbv
    for bv in basisSET:
        for mM in mLmJl:
            tmp = sourceRHS[('%d' % (bv_number), '%d' % (2 * Jstreu),
                             '%d' % (2 * mM[1]), '%d' % int(2 * multipolarity),
                             '%d' % (2 * mM[0]))]
            cgtmp = CG(multipolarity, mM[0], J0, mM[1] - mM[0], Jstreu,
                       mM[1]).doit()
            RHSofmJsk[('%d' % (bv_number), '%d' % (2 * Jstreu), '%d' %
                       (2 * mM[1]), '%d' % (2 * multipolarity))] += tmp * cgtmp
        bv_number += 1

    for mJ in mJlrange:

        outs = ''
        for nMom in range(anz_phot_e):

            for bv in range(bv_number):
                outs += '%12.12E  ' % float(RHSofmJsk[(
                    '%d' % (bv), '%d' % (2 * Jstreu), '%d' % (2 * mJ),
                    '%d' % (2 * multipolarity))][nMom])

            outs += '\n'

        outp = respath + 'LIT_SOURCE_%s_mJ' % streukanal + '{:,g}'.format(
            mJ) + '_MUL%d_bv_1-%d' % (multipolarity, bv_number)

        print('writing source column to <%s>' % outp)
        with open(outp, 'w') as outfile:

            #outfile.seek(0)
            outfile.write(outs)

    return  #coupledSOURCE