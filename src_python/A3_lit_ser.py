import os, sys
import multiprocessing
import subprocess
import shlex
from multiprocessing.pool import ThreadPool

import numpy as np
from scipy.optimize import fmin
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from bridge import *
from rrgm_functions import *
from three_particle_functions import *
from triton_width_gen import *
from readLITsource_3body import *
from BasisVisualization import visbas
from plot_spectrum import plotHspec

RHSofBV = {}
RHSofmJ = {}
siffux = '_ref'

if os.path.isdir(litpath3He) == False:
    os.mkdir(litpath3He)
os.chdir(litpath3He)

with open(litpath3He + 'kRange.dat', 'wb') as f:
    np.savetxt(f, [anz_phot_e, phot_e_0, phot_e_d], fmt='%f')
f.close()

if 'construe_fresh_helion' in cal:

    os.chdir(v18uixpath)
    print('(working dir) %s' % v18uixpath)

    os.system('cp INLUCN%s INLUCN' % siffux)
    os.system(BINBDGpath + 'LUDW_CN.exe')
    os.system('cp INLU%s INLU' % siffux)
    os.system(BINBDGpath + 'DRLUD.exe')
    os.system('cp INOB%s INOB' % siffux)
    os.system(BINBDGpath + 'KOBER.exe')
    os.system(BINBDGpath + 'DROBER.exe')
    os.system('cp INQUA_N%s INQUA_N' % siffux)
    repl_line('INQUA_N', 1, potnn + '\n')
    os.system(BINBDGpath + 'QUAFL_N.exe')
    repl_line('INQUA_N', 1, potnnn + '\n')
    os.system(BINBDGpath + 'DRQUA_AK_N.exe')
    os.system('cp INEN%s INEN' % siffux)
    os.system(BINBDGpath + 'DR2END_AK.exe')

    EBDG = get_h_ev()[0]
    np.savetxt('E0.dat', np.array([EBDG]), fmt='%12.4f')

    os.system('cp OUTPUT end_out_b && cp INEN inen_b')
    os.system('cat E0.dat')

    rrgm_functions.parse_ev_coeffs(infil='end_out_b')

    print('helion ground state calculated!')

if 'rhs' in cal:

    os.chdir(litpath3He)

    for streukanal in streukas:

        Jstreu = float(streukanal.split('^')[0])
        Jstreustring = '%s' % str(Jstreu)[:3]

        mLmJl, mLrange, mJlrange = non_zero_couplings(multipolarity, J0,
                                                      Jstreu)

        # ECCE
        BUECO = [cof.strip() for cof in open(v18uixpath + 'COEFF')]
        EBDG = get_h_ev(ifi=v18uixpath + 'end_out_b')[0]

        print('(iv)    LS-scheme: B(2,%s) = %4.4f MeV [' % (boundstatekanal,
                                                            EBDG),
              get_h_ev(n=4, ifi=v18uixpath + 'end_out_b'), ']')
        print('        dim(B_0)   = %d' % len(BUECO))

        if 'rhs_lu-ob-qua' in cal:

            #            he_iw = [
            #                np.array(ln.split(';')).astype(float).tolist() for ln in open(
            #                    litpath3He + 'basis_struct/intw3heLIT_J%s_%s.dat' %
            #                    (Jstreustring, boundstatekanal))
            #            ]
            #
            #            he_rw = [
            #                np.array(ln.split(';')).astype(float).tolist() for ln in open(
            #                    litpath3He + 'basis_struct/relw3heLIT_J%s_%s.dat' %
            #                    (Jstreustring, boundstatekanal))
            #            ]
            he_iw, he_rw, frgs = retrieve_he3_widths(
                v18uixpath + 'INQUA_N_ref')
            lfrags = []
            sfrags = []

            for lcfg in range(len(channels[boundstatekanal])):
                sfrags = sfrags + channels[boundstatekanal][lcfg][1]
                for scfg in channels[boundstatekanal][lcfg][1]:
                    lfrags = lfrags + [channels[boundstatekanal][lcfg][0]]
#            fragfile = [
#                ln for ln in
#                open(litpath3He + 'basis_struct/frags_LIT_J%s_%s.dat' %
#                     (Jstreustring, boundstatekanal))
#            ]
#            lfrags = [fr.split(' ')[1].strip() for fr in fragfile]
#            sfrags = [fr.split(' ')[0] for fr in fragfile]

# read widths and frags of the LIT basis as determined via
# v18uix_LITbasis.py
            fragfile = [
                ln for ln in
                open(litpath3He + 'basis_struct/frags_LIT_J%s_%s.dat' %
                     (Jstreustring, streukanal))
            ]
            lfrags2 = [fr.split(' ')[1].strip() for fr in fragfile]
            sfrags2 = [fr.split(' ')[0] for fr in fragfile]

            intwLIT = [
                np.array(ln.split(';')).astype(float).tolist() for ln in open(
                    litpath3He + 'basis_struct/intw3heLIT_J%s_%s.dat' %
                    (Jstreustring, streukanal))
            ]

            relwLIT = [
                np.array(ln.split(';')).astype(float).tolist() for ln in open(
                    litpath3He + 'basis_struct/relw3heLIT_J%s_%s.dat' %
                    (Jstreustring, streukanal))
            ]

            if 'dbg' in cal:
                print(
                    '\n3He components (full) + LIT-basis components (bare):\n',
                    len(lfrags))
                print(sfrags)
                print('\nLIT-basis components (full):\n', len(lfrags2))
                print(sfrags2)

            for file in os.listdir(os.getcwd()):
                if fnmatch.fnmatch(file, 'endlit*'):
                    if 'dbg' in cal:
                        print('removing old <*en*lit*> files.')
                    os.system('rm *en*lit*')
                    break
            for file in os.listdir(os.getcwd()):
                if fnmatch.fnmatch(file, 'inhomo*'):
                    if 'dbg' in cal:
                        print('removing old <*inhomo*> files.')
                    os.system('rm inhomo*')
                    break

            lit_3inqua_seq(
                intwi=he_iw + intwLIT,
                relwi=he_rw + relwLIT,
                anzo=11,
                LREG='  1  0  0  0  0  0  0  0  0  1  1',
                outfile=litpath3He + 'INQUA')
            lit_3inlu(
                mul=multipolarity,
                frag=lfrags + lfrags2,
                fn=litpath3He + 'INLU')
            lit_3inob(fr=sfrags + sfrags2, fn=litpath3He + 'INOB')

            os.system(BINLITpath + 'luise.exe > dump')
            os.system(BINLITpath + 'obem.exe > dump')
            os.system(BINLITpath + 'qual.exe')

            exit()

        def cal_rhs_end(para, procnbr):

            inenf = 'inenlit%d-%d_J%3.1f_mJ%3.1f-mL%d.log' % (para[1], para[2],
                                                              Jstreu,
                                                              para[0][1],
                                                              para[0][0])

            outfseli = 'endlit%d-%d_J%3.1f_mJ%3.1f-mL%d.log' % (
                para[1],
                para[2],
                Jstreu,
                para[0][1],
                para[0][0],
            )
            outfsbare = 'inhomo%d-%d_J%3.1f_mJ%3.1f-mL%d.log' % (
                para[1],
                para[2],
                Jstreu,
                para[0][1],
                para[0][0],
            )

            leftpar = int(1 + 0.5 * (
                1 + (-1)**(int(channels[streukanal][0][0][0]
                               ) + int(channels[streukanal][0][0][1]))))
            #print('I, %d of %d, will do something.\n' % (procnbr, anzproc),
            #      leftpar)

            lit_3inen(
                MREG='  1  0  0  0  0  0  0  0  0  1  1',
                #                   (shifted) QBV                     nr.rw
                KSTREU=[para[1], para[2]],
                JWSL=Jstreu,
                JWSLM=para[0][1],
                MULM2=para[0][0],
                NPARL=leftpar,
                JWSR=J0,
                NPARR=2,
                EB=EBDG,
                BUECO=BUECO,
                NZE=anz_phot_e,
                EK0=phot_e_0,
                EKDIFF=phot_e_d,
                bnd=v18uixpath + 'INEN%s' % siffux,
                outfile=inenf)

            cmdend = BINLITpath + 'enemb.exe %s %s %s' % (inenf, outfseli,
                                                          outfsbare)

            pend = subprocess.Popen(
                shlex.split(cmdend),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=os.getcwd())
            out, err = pend.communicate()

            return (out, err)

        parameter_set_end = []

        wfn = litpath3He + 'basis_struct/LITbas_full_J%s_%s.dat' % (
            Jstreustring, streukanal)

        print('[...] reading BV-rw tupel from %s' % wfn)
        litbas = [bv for bv in np.loadtxt(wfn).astype(int) if bv[1] != 0]

        parameter_set = []

        he_iw, he_rw, frgs = retrieve_he3_widths(v18uixpath + 'INQUA_N_ref')

        anzhebv = sum([len(z) for z in he_iw])

        for mM in mLmJl:
            for bv in litbas:
                parameter_set.append([mM, bv[0] + anzhebv, bv[1]])
        parameter_set_end = parameter_set

        if 'rhs-end' in cal:

            results = []
            pool = ThreadPool(anzproc)
            parameter_set = parameter_set_end
            for procnbr in range(len(parameter_set)):
                pars = parameter_set[procnbr]
                results.append(
                    pool.apply_async(cal_rhs_end, (
                        pars,
                        procnbr,
                    )))

            pool.close()
            pool.join()

            exit()

if 'lhs' in cal:

    print('(ii)    calculating norm/ham in scattering-channel basis')

    if os.path.isdir(litpath3He + 'lit_bas_lhs/') == False:
        os.mkdir(litpath3He + 'lit_bas_lhs/')
    os.chdir(litpath3He + 'lit_bas_lhs/')

    for streukanal in streukas:

        Jstreu = float(streukanal.split('^')[0])
        Jstreustring = '%s' % str(Jstreu)[:3]

        fragfile = [
            ln
            for ln in open(litpath3He + 'basis_struct/frags_LIT_J%s_%s.dat' %
                           (Jstreustring, streukanal))
        ]

        lfrags = [fr.split(' ')[1].strip() for fr in fragfile]
        sfrags = [fr.split(' ')[0] for fr in fragfile]

        intwLIT = [
            np.array(ln.split(';')).astype(float)
            for ln in open(litpath3He + 'basis_struct/intw3heLIT_J%s_%s.dat' %
                           (Jstreustring, streukanal))
        ]

        anzLITbv = sum([len(frgm) for frgm in intwLIT])

        if (len([len(frgm) for frgm in intwLIT]) != (len(fragfile))):
            print('LIT-basis fragments inconcistent!',
                  len([len(frgm) for frgm in intwLIT]), (len(fragfile)))
            exit()

        relwLIT = [
            np.array(ln.split(';')).astype(float)
            for ln in open(litpath3He + 'basis_struct/relw3heLIT_J%s_%s.dat' %
                           (Jstreustring, streukanal))
        ]

        if 'dbg' in cal:
            print(lfrags, sfrags)
            #print(intwLIT, relwLIT)

        Jstreu = float(streukanal.split('^')[0])

        if 'lhs_lu-ob-qua' in cal:

            os.chdir(litpath3He + 'lit_bas_lhs/')

            n3_inlu(8, fn='INLU', fr=lfrags, indep=-0)
            os.system(BINBDGpath + 'DRLUD.exe')
            n3_inlu(8, fn='INLUCN', fr=lfrags, indep=-0)
            os.system(BINBDGpath + 'LUDW_CN.exe')

            n3_inob(sfrags, 8, fn='INOB', indep=-0)
            os.system(BINBDGpath + 'KOBER.exe')
            n3_inob(sfrags, 15, fn='INOB', indep=-0)
            os.system(BINBDGpath + 'DROBER.exe')

            he3inqua(intwi=intwLIT, relwi=relwLIT, potf=potnn)
            os.system(BINBDGpath + 'QUAFL_N.exe')

            he3inqua(intwi=intwLIT, relwi=relwLIT, potf=potnnn)
            os.system(BINBDGpath + 'DRQUA_AK_N.exe')

        litbas = np.loadtxt(
            litpath3He + 'basis_struct/LITbas_full_J%s_%s.dat' %
            (Jstreustring, streukanal)).astype(int)

        litbas = [bv for bv in np.unique(litbas, axis=0) if bv[1] != 0]

        n3_inen_rhs(
            litbas,
            Jstreu,
            costr,
            np.ones(len(relwLIT[0])),
            fn='INEN',
            pari=0,
            nzop=31,
            tni=11)

        os.system(
            BINBDGpath + 'DR2END_AK.exe && grep -A 3 \'EIGENWER\' OUTPUT')

        os.system('cp INEN inen-lit-%s' % streukanal)

        os.system('cp %s/MATOUT %s/norm-ham-litME-%s' %
                  (litpath3He + 'lit_bas_lhs/', litpath3He + 'lit_bas_lhs/',
                   streukanal))
        plotHspec(Jstreustring)

if 'couple' in cal:

    os.system('cp ' + v18uixpath + 'E0.dat ' + litpath3He)

    for streukanal in streukas:

        Jstreu = float(streukanal.split('^')[0])
        Jstreustring = '%s' % str(Jstreu)[:3]

        fragfile = [
            ln
            for ln in open(litpath3He + 'basis_struct/frags_LIT_J%s_%s.dat' %
                           (Jstreustring, streukanal))
        ]

        lfrags2 = [fr.split(' ')[1].strip() for fr in fragfile]

        RHSofBV[streukanal] = {}
        RHSofmJ[streukanal] = {}
        photEn = []

        # read uncoupled source ME's
        os.chdir(litpath3He)

        litbas = np.loadtxt(
            litpath3He + 'basis_struct/LITbas_full_J%s_%s.dat' %
            (Jstreustring, streukanal)).astype(int)

        litbas = [bv for bv in np.unique(litbas, axis=0) if bv[1] != 0]

        read_uncoupled_source(
            RHSofBV[streukanal],
            photEn,
            streukanal,
            basisSET=litbas,
            firstbv=0)
        print(RHSofBV[streukanal])
        exit()
        #RHSofBV[streukanal] = read_norm(streukanal, basisSET=litbas)
        #photEn = MeVfm * np.array(
        #    [phot_e_0 + en * phot_e_d for en in range(anz_phot_e)])
        # couple incoming state with photon multipole to Jlit
        couple_source(
            RHSofmJ[streukanal],
            streukanal,
            RHSofBV[streukanal],
            basisSET=litbas,
            firstbv=bv_offset)
        bv_offset += len(litbas)
        if 'plt' in cal:
            fig = plt.figure(figsize=(12, 6))
            #fig.subplots_adjust(hspace=1.4, wspace=0.4)
            #for i in range(len(streukas)):
            i = 0
            ax1 = fig.add_subplot(1, 2, 1)
            ax1.set_title(r'$J^\pi=%s^%s$' % (Jstreustring, streukanal[-1]))
            ax1.set_xlabel('photon momentum [MeV]')
            #ax1.set_title(r'$J^\pi=%d^%s$' % (Jstreu, streukas[i][-1]))
            mLmJl, mLrange, mJlrange = non_zero_couplings(
                multipolarity, J0, Jstreu)
            mM = mLmJl[0]
            #    for bv in litbas:
            #        print(photEn)
            #        print(RHSofmJ[streukas[i]][('%d-%d' % (bv[0], bv[1]),
            #                                    '%d' % (2 * Jstreu), '%d' % (2 * mM[1]),
            #                                    '%d' % (2 * multipolarity),
            #                                    '%d' % (2 * mM[0]))].astype(float), bv[0],
            #              bv[1], Jstreu, mM)
            #        exit()
            [
                ax1.plot(
                    photEn,
                    RHSofBV[streukanal][('%d-%d' % (bv[0], bv[1]), '%d' %
                                         (2 * Jstreu), '%d' % (2 * mM[1]),
                                         '%d' % (2 * multipolarity),
                                         '%d' % (2 * mM[0]))].astype(float))
                for bv in litbas
            ]
            ax2 = fig.add_subplot(1, 2, 2)
            ax2.set_xlabel('photon momentum [MeV]')
            ax2.set_ylabel(r'$\left\langle\,Jm\,\vert\,Jm\,\right\rangle$ [-]')
            ax2.set_title(r'$J$-coupled RHS')
            [
                ax2.plot(
                    photEn,
                    np.array(RHSofmJ[streukanal][(
                        '%d-%d' % (bv[0], bv[1]), '%d' % (2 * Jstreu),
                        '%d' % (2 * mM[1]),
                        '%d' % (2 * multipolarity))].astype(float)),
                    label=r'$BV_{%d}W^{%d}$' % (int(bv[0]), int(bv[1])))
                for bv in litbas
            ]
            box1 = ax1.get_position()
            ax1.set_position([box1.x0, box1.y0, box1.width, box1.height])
            box2 = ax2.get_position()
            ax2.set_position([box2.x0, box2.y0, box2.width, box2.height])
            axins = inset_axes(
                ax2,
                width="60%",
                height="75%",
                bbox_to_anchor=(.3, .1, .95, .35),
                bbox_transform=ax2.transAxes,
                loc=3)
            [
                axins.plot(
                    photEn,
                    np.array(RHSofmJ[streukanal][(
                        '%d-%d' % (bv[0], bv[1]), '%d' % (2 * Jstreu),
                        '%d' % (2 * mM[1]),
                        '%d' % (2 * multipolarity))].astype(float)),
                    label=r'$BV_{%d}W^{%d}$' % (int(bv[0]), int(bv[1])))
                for bv in litbas
            ]
            # sub region of the original image
            x1, x2, y1, y2 = photEn[0], photEn[-1], -0.1, 0.1
            axins.set_xlim(x1, x2)
            axins.set_ylim(y1, y2)
            #axins.set_xticklabels('')
            #axins.set_yticklabels('')
            #ax2.indicate_inset_zoom(axins)
            fig.savefig('LITrhs_J%s.pdf' % (Jstreustring))
            fig_leg = plt.figure(figsize=(10, 8), dpi=95)
            ax_leg = fig_leg.add_subplot(111)
            # add the legend from the previous axes
            ax_leg.legend(
                *ax2.get_legend_handles_labels(),
                fontsize=10,
                ncol=8,
                loc='center')
            # hide the axes frame and the x/y labels
            ax_leg.axis('off')
            ax_leg.set_title(r'coupled-source legend')
            fig_leg.savefig('LITrhs_legend_J%s.pdf' % (Jstreustring))
            print('RHS vector visualized in <LITrhs_Jxx.pdf>')