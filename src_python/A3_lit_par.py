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

os.chdir(litpath3He)

if os.path.isfile(litpath3He + 'kRange.dat') == True:
    os.system('rm ' + litpath3He + 'kRange.dat')

with open(litpath3He + 'kRange.dat', 'wb') as f:
    np.savetxt(f, [anz_phot_e, phot_e_0, phot_e_d], fmt='%f')
f.close()

siffux = '_ref'
he_iw, he_rw, he_frgs = retrieve_he3_widths(
    helionpath + 'INQUA_V18%s' % siffux)

if 'construe_fresh_helion' in cal:

    os.chdir(helionpath)
    print('(working dir) %s' % helionpath)

    os.system('cp INLUCN%s INLUCN' % siffux)
    os.system(BINBDGpath + 'LUDW_CN.exe')
    os.system('cp INLU%s INLU' % siffux)
    os.system(BINBDGpath + 'DRLUD.exe')
    os.system('cp INOB%s INOB' % siffux)
    os.system(BINBDGpath + 'KOBER.exe')
    os.system('cp DRINOB%s INOB' % siffux)
    os.system(BINBDGpath + 'DROBER.exe')
    os.system('cp INQUA_V18%s INQUA_N' % siffux)
    repl_line('INQUA_N', 1, potnn + '\n')

    subprocess.run([
        'mpirun', '-np',
        '%d' % anzproc,
        home + '/kette_repo/ComptonLIT/src_nucl/V18_PAR/mpi_quaf_v6'
    ])
    subprocess.run([home + '/kette_repo/ComptonLIT/src_nucl/V18_PAR/sammel'])
    #os.system(BINBDGpath + 'QUAFL_N.exe')

    os.system('cp INQUA_UIX%s INQUA_N' % siffux)

    subprocess.run([
        'mpirun', '-np',
        '%d' % anzproc,
        home + '/kette_repo/ComptonLIT/src_nucl/UIX_PAR/mpi_drqua_uix'
    ])
    subprocess.run(
        [home + '/kette_repo/ComptonLIT/src_nucl/UIX_PAR/SAMMEL-uix'])
    #os.system(BINBDGpath + 'DRQUA_AK_N.exe')

    os.system('cp INEN%s INEN' % siffux)
    subprocess.run([home + '/kette_repo/ComptonLIT/src_nucl/TDR2END_AK.exe'])
    #os.system(BINBDGpath + 'DR2END_AK.exe')

    EBDG = get_h_ev()[0]

    np.savetxt('E0.dat', np.array([EBDG]), fmt='%12.4f')

    os.system('cp OUTPUT end_out_b && cp INEN inen_b')
    os.system('cat E0.dat')

    #rrgm_functions.parse_ev_coeffs(infil='end_out_b')
    rrgm_functions.parse_ev_coeffs_normiert(infil='end_out_b')

    subprocess.call('rm *QUAOUT*', shell=True)

    print('helion ground state calculated with B = %4.4f MeV' % EBDG)

    os.chdir(litpath3He)
    os.system('find . -name \"T*OUT.*\" -print0 | xargs -0 rm')
    os.system('find . -name \"D*OUT.*\" -print0 | xargs -0 rm')
    os.system('find . -name \"T*OUT.*\" -print0 | xargs -0 rm')
    os.system('find . -name \"OUTPUT\" -print0 | xargs -0 rm')
    os.system('find . -name \"MATOUT\" -print0 | xargs -0 rm')
    os.system('find . -name \"*OUT*\" -print0 | xargs -0 rm')

if 'rhs' in cal:

    os.chdir(litpath3He)

    for streukanal in streukas:

        Jstreu = float(streukanal.split('^')[0])
        Jstreustring = '%s' % str(Jstreu)[:3]

        mLmJl, mLrange, mJlrange = non_zero_couplings(multipolarity, J0,
                                                      Jstreu)

        # ECCE
        #parse_ev_coeffs(
        #    mult=0, infil=helionpath + 'OUTPUT', outf=helionpath + 'COEFF')
        #parse_ev_coeffs_normiert(infil=helionpath + 'OUTPUT', )

        #BUECO = np.array(
        #    [float(cof.strip()) for cof in open(helionpath + 'COEFF')])

        BUECO = (10**-4) * np.array(
            [float(cof.strip()) for cof in open(helionpath + 'COEFF_NORMAL')])

        EBDG = get_h_ev(ifi=helionpath + 'end_out_b')[0]

        print('(iv)    LS-scheme: B(2,%s) = %4.4f MeV [' % (boundstatekanal,
                                                            EBDG),
              get_h_ev(n=4, ifi=helionpath + 'end_out_b'), ']')
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

            lfrags = []
            sfrags = []

            #for lcfg in range(len(channels[boundstatekanal])):
            #    sfrags = sfrags + channels[boundstatekanal][lcfg][1]
            #    for scfg in channels[boundstatekanal][lcfg][1]:
            #        lfrags = lfrags + [channels[boundstatekanal][lcfg][0]]
            fragfile = [
                ln for ln in
                open(litpath3He + 'basis_struct/frags_LIT_J%s_%s.dat' %
                     (J0, boundstatekanal))
            ]
            lfrags = [fr.split(' ')[1].strip() for fr in fragfile]
            sfrags = [fr.split(' ')[0] for fr in fragfile]

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

            for lit_zerl in range(len(lfrags2)):

                if os.path.isdir(litpath3He + 'tmp_%d' % lit_zerl) == False:
                    os.mkdir(litpath3He + 'tmp_%d' % lit_zerl)
                os.chdir(litpath3He + 'tmp_%d' % lit_zerl)

                for file in os.listdir(os.getcwd()):
                    if fnmatch.fnmatch(file, '*J%s*.log' % Jstreu):
                        print('I would have deleted .log files...')
                        #exit()
                        if 'dbg' in cal:
                            print('removing old <*.log> files.')
                        os.system('rm *.log')
                        break

                lit_3inqua_seq(
                    intwi=he_iw + [intwLIT[lit_zerl]],
                    relwi=he_rw + [relwLIT[lit_zerl]],
                    anzo=11,
                    LREG='  1  0  0  0  0  0  0  0  0  1  1',
                    outfile=litpath3He + 'tmp_%d/INQUA' % (lit_zerl))
                lit_3inlu(
                    mul=multipolarity,
                    frag=lfrags + [lfrags2[lit_zerl]],
                    fn=litpath3He + 'tmp_%d/INLU' % (lit_zerl))
                lit_3inob(
                    fr=sfrags + [sfrags2[lit_zerl]],
                    fn=litpath3He + 'tmp_%d/INOB' % (lit_zerl))

        leftpar = int(1 + 0.5 * (1 + (-1)**
                                 (int(channels[streukanal][0][0][0]
                                      ) + int(channels[streukanal][0][0][1]))))

        def cal_rhs_lu_ob_qua(para, procnbr):

            slave_pit = litpath3He + 'tmp_%d' % para
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
            print('process = %d-1 : luise exits.' % para)

            print('%s in %s' % (cmdob, slave_pit))
            pob = subprocess.Popen(
                shlex.split(cmdob),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=slave_pit)
            out, err = pob.communicate()
            print('process = %d-1 : ober exits.' % para)

            print('%s in %s' % (cmdqu, slave_pit))
            pqu = subprocess.Popen(
                shlex.split(cmdqu),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=slave_pit)
            out, err = pqu.communicate()
            print('process = %d-1 : qual exits.' % para)

        def cal_rhs_end(para, procnbr):

            slave_pit = litpath3He + 'tmp_%d/' % para[3]

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
                bnd=helionpath + 'INEN',
                outfile=slave_pit + inenf)

            cmdend = BINLITpath + 'enemb.exe %s %s %s' % (inenf, outfseli,
                                                          outfsbare)

            pend = subprocess.Popen(
                shlex.split(cmdend),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=slave_pit)
            out, err = pend.communicate()

            return (out, err)

        parameter_set_lu_ob_qua = range(len(lfrags2))
        parameter_set_end = []

        wfn = litpath3He + 'basis_struct/LITbas_full_J%s_%s.dat' % (
            Jstreustring, streukanal)

        print('[...] reading BV-rw tupel from %s' % wfn)
        litbas = np.loadtxt(wfn).astype(int)

        for lit_zerl in range(len(lfrags2)):
            bsbv = sum([len(b) for b in he_iw])
            parameter_set = []
            bvrange = range(
                sum([len(z) for z in intwLIT[:lit_zerl]]) + 1,
                sum([len(z) for z in intwLIT[:lit_zerl]]) + 1 +
                len(intwLIT[lit_zerl]))

            litbas3 = []
            for bv in filter(lambda x: (x[0] in bvrange), litbas):
                litbas3.append([int(bv[0] - (bvrange[0] - 1) + bsbv), bv[1]])

            with open(litpath3He + 'tmp_%d/LITbas_full_J%s.dat' %
                      (lit_zerl, Jstreustring), 'wb') as f:
                np.savetxt(f, [[jj[0], jj[1]] for jj in litbas3], fmt='%d')
            f.close()

            for mM in mLmJl:
                for bv in filter(lambda x: (x[0] in bvrange), litbas):
                    parameter_set.append([
                        mM,
                        int(bv[0] - (bvrange[0] - 1) + bsbv), bv[1], lit_zerl
                    ])

            parameter_set_end.append(parameter_set)

        if 'rhs-qual' in cal:
            parameter_set = parameter_set_lu_ob_qua
            results = []
            pool = ThreadPool(anzproc)
            for procnbr in range(len(parameter_set)):
                pars = parameter_set[procnbr]
                results.append(
                    pool.apply_async(cal_rhs_lu_ob_qua, (
                        pars,
                        procnbr,
                    )))
            pool.close()
            pool.join()

            for lit_zerl in range(len(lfrags2)):
                os.system('mv ' + litpath3He + 'tmp_%d/QUAOUT ' % lit_zerl +
                          litpath3He + 'tmp_%d/QUAOUT_J%3.1f' % (lit_zerl,
                                                                 Jstreu))

        if 'rhs-end' in cal:
            for lit_zerl in range(len(lfrags2)):
                print('(J=%s)  werkle in %d' % (Jstreu, lit_zerl))

                try:
                    os.system('cp ' + litpath3He + 'tmp_%d/QUAOUT_J%3.1f ' % (
                        lit_zerl,
                        Jstreu) + litpath3He + 'tmp_%d/QUAOUT' % (lit_zerl))
                except:
                    print('<QUAOUT> na for this channel.')
                    exit()

                results = []
                pool = ThreadPool(anzproc)
                parameter_set = parameter_set_end[lit_zerl]
                for procnbr in range(len(parameter_set)):

                    pars = parameter_set[procnbr]

                    results.append(
                        pool.apply_async(cal_rhs_end, (
                            pars,
                            procnbr,
                        )))

                pool.close()
                pool.join()

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

            os.chdir(pathbase + '/data/eob/')
            n3_inob(
                [
                    'he_no1', 'he_no1i', 'he_no2', 'he_no2i', 'he_no6',
                    'he_no6i', 'he_no3i', 'he_no5i'
                ],
                8,
                fn='INOB',
                indep=+1)
            os.system(BINBDGpath + 'KOBER.exe')
            os.chdir(pathbase + '/data/eob-tni/')
            n3_inob(
                [
                    'he_no1', 'he_no1i', 'he_no2', 'he_no2i', 'he_no6',
                    'he_no6i', 'he_no3i', 'he_no5i'
                ],
                15,
                fn='INOB',
                indep=+1)
            os.system(BINBDGpath + 'DROBER.exe')

            os.chdir(pathbase + '/data/elu/')
            n3_inlu(
                8,
                fn='INLUCN',
                fr=[
                    '000', '202', '022', '110', '101', '011', '111', '112',
                    '211', '121', '122', '212', '222', '221', '220'
                ],
                indep=+1)
            os.system(BINBDGpath + 'LUDW_CN.exe')
            os.chdir(pathbase + '/data/elu-tni/')
            n3_inlu(
                8,
                fn='INLU',
                fr=[
                    '000', '202', '022', '110', '101', '011', '111', '112',
                    '211', '121', '122', '212', '222', '221', '220'
                ],
                indep=+1)
            os.system(BINBDGpath + 'DRLUD.exe')

            os.chdir(litpath3He + 'lit_bas_lhs/')

            n3_inlu(8, fn='INLU', fr=lfrags, indep=-0)
            os.system(BINBDGpath + 'DRLUD.exe')
            n3_inlu(8, fn='INLUCN', fr=lfrags, indep=-0)
            os.system(BINBDGpath + 'LUDW_CN.exe')

            n3_inob(sfrags, 8, fn='INOB', indep=-0)
            os.system(BINBDGpath + 'KOBER.exe')
            n3_inob(sfrags, 15, fn='INOB', indep=-0)
            os.system(BINBDGpath + 'DROBER.exe')

            insam(len(lfrags))

            he3inqua(intwi=intwLIT, relwi=relwLIT, potf=potnn)

            os.system('cp INQUA_N INQUA_N%s' % boundstatekanal)

            parallel_mod_of_3inqua(
                lfrags, sfrags, infile='INQUA_N', outfile='INQUA_N')

            subprocess.run([
                'mpirun', '-np',
                '%d' % anzproc, pathbase + '/src_nucl/V18_PAR/mpi_quaf_v6'
            ])

            subprocess.run([pathbase + '/src_nucl/V18_PAR/sammel'])

            he3inqua(intwi=intwLIT, relwi=relwLIT, potf=potnnn)

            parallel_mod_of_3inqua(
                lfrags, sfrags, infile='INQUA_N', outfile='INQUA_N', tni=1)

            subprocess.run([
                'mpirun', '-np',
                '%d' % anzproc, pathbase + '/src_nucl/UIX_PAR/mpi_drqua_uix'
            ])

            subprocess.run([pathbase + '/src_nucl/UIX_PAR/SAMMEL-uix'])

        litbas = np.loadtxt(
            litpath3He + 'basis_struct/LITbas_full_J%s_%s.dat' %
            (Jstreustring, streukanal)).astype(int)
        litbas = [bv for bv in np.unique(litbas, axis=0) if bv[1] != 0]

        anzbs = []

        for lit_zerl in range(len(lfrags)):
            bvmax = 0
            ntmp = 0
            for mm in range(lit_zerl + 1):
                bvmax += len(intwLIT[mm])

            for bv in litbas:
                if bv[0] <= bvmax:
                    ntmp += 1
                else:
                    continue
            anzbs.append(ntmp)

        if len(anzbs) != len(lfrags):
            print('Basis blocks inconsistent. Exiting...')
            exit()

        lit_zerl = 0

        for anzbtmp in anzbs:

            mypath = litpath3He + 'tmp_%d/' % lit_zerl
            lit_zerl += 1
            n3_inen_rhs(
                litbas,
                Jstreu,
                costr,
                np.ones(len(relwLIT[0])),
                fn='INEN',
                pari=0,
                nzop=31,
                tni=11,
                anzb=anzbtmp)
            subprocess.run([pathbase + '/src_nucl/TDR2END_AK.exe'])
            os.system('cp INEN ' + mypath + 'inen-lit-%s_1-%d' % (streukanal,
                                                                  anzbtmp))
            os.system('cp %s/MATOUT ' % (litpath3He + 'lit_bas_lhs/') + respath
                      + 'norm-ham-litME-%s_1-%d' % (streukanal, anzbtmp))

        #plotHspec(Jstreustring)

if 'couple' in cal:

    os.system('cp ' + helionpath + 'E0.dat ' + litpath3He)

    for streukanal in streukas:

        Jstreu = float(streukanal.split('^')[0])
        Jstreustring = '%s' % str(Jstreu)[:3]

        intwLIT = [
            np.array(ln.split(';')).astype(float).tolist()
            for ln in open(litpath3He + 'basis_struct/intw3heLIT_J%s_%s.dat' %
                           (Jstreustring, streukanal))
        ]

        bv_pro_zerl = [len(zset) for zset in intwLIT]
        first_bv_pro_zerl = [
            sum(bv_pro_zerl[:n]) + 1 for n in range(len(bv_pro_zerl))
        ]

        relwLIT = [
            np.array(ln.split(';')).astype(float).tolist()
            for ln in open(litpath3He + 'basis_struct/relw3heLIT_J%s_%s.dat' %
                           (Jstreustring, streukanal))
        ]

        fragfile = [
            ln
            for ln in open(litpath3He + 'basis_struct/frags_LIT_J%s_%s.dat' %
                           (Jstreustring, streukanal))
        ]

        lfrags = [fr.split(' ')[1].strip() for fr in fragfile]

        litbas = np.loadtxt(
            litpath3He + 'basis_struct/LITbas_full_J%s_%s.dat' %
            (Jstreustring, streukanal)).astype(int)
        litbas = [bv for bv in np.unique(litbas, axis=0) if bv[1] != 0]

        print('He3 structure: anz[int,rel]weiten:', np.sum(he_frgs, axis=0))

        bv_offset = 0

        RHSofBV[streukanal] = {}
        RHSofmJ[streukanal] = {}

        for lit_zerl in range(len(lfrags)):

            # read uncoupled source ME's
            os.chdir(litpath3He + 'tmp_%d' % lit_zerl)

            litbas = np.loadtxt(
                litpath3He + 'basis_struct/LITbas_full_J%s_%s.dat' %
                (Jstreustring, streukanal)).astype(int)
            litbas = [bv for bv in np.unique(litbas, axis=0) if bv[1] != 0]

            print('werkle in %s' % os.getcwd())
            tmp_iw, tmp_rw, tmp_frgs = retrieve_he3_widths('INQUA')

            he_bv = np.sum(he_frgs, axis=0)[0]
            print('bv in He3 : %d' % he_bv)
            bv_in_zerl = tmp_frgs[-1][0]

            print('bv in zerl : %d' % bv_in_zerl)

            zerl_bas = [
                bv for bv in litbas if first_bv_pro_zerl[lit_zerl] <= bv[0] <
                first_bv_pro_zerl[lit_zerl + 1]
            ] if (lit_zerl < len(lfrags) - 1) else [
                bv for bv in litbas if first_bv_pro_zerl[lit_zerl] <= bv[0]
            ]

            zerl_bas = [[
                bv[0] + 1 + he_bv - first_bv_pro_zerl[lit_zerl], bv[1]
            ] for bv in zerl_bas]

            photEn = []

            read_uncoupled_source(
                RHSofBV[streukanal],
                photEn,
                streukanal,
                basisSET=zerl_bas,
                firstbv=bv_offset)
            #RHSofBV[streukanal] = read_norm(streukanal, basisSET=litbas)
            #photEn = MeVfm * np.array(
            #    [phot_e_0 + en * phot_e_d for en in range(anz_phot_e)])
            # couple incoming state with photon multipole to Jlit
            couple_source(
                RHSofmJ[streukanal],
                streukanal,
                RHSofBV[streukanal],
                basisSET=zerl_bas,
                firstbv=bv_offset)

            bv_offset += len(zerl_bas)

            if 'plt' in cal:
                fig = plt.figure(figsize=(12, 6))

                #fig.subplots_adjust(hspace=1.4, wspace=0.4)
                #for i in range(len(streukas)):
                i = 0

                ax1 = fig.add_subplot(1, 2, 1)
                ax1.set_title(r'$J^\pi=%s^%s$' % (Jstreustring,
                                                  streukanal[-1]))
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
                                             '%d' % (2 * multipolarity), '%d' %
                                             (2 * mM[0]))].astype(float))
                    for bv in litbas
                ]

                ax2 = fig.add_subplot(1, 2, 2)
                ax2.set_xlabel('photon momentum [MeV]')
                ax2.set_ylabel(
                    r'$\left\langle\,Jm\,\vert\,Jm\,\right\rangle$ [-]')
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

os.chdir(litpath3He)

os.system('find . -name \"T*OUT.*\" -print0 | xargs -0 rm')