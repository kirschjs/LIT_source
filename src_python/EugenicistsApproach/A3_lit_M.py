import os, sys
import multiprocessing
import subprocess
import shlex
import glob
from multiprocessing.pool import ThreadPool

import numpy as np

from bridgeA3 import *
from rrgm_functions import *
from three_particle_functions import *

# CG(j1, m1, j2, m2, j3, m3)
from sympy.physics.quantum.cg import CG
from scipy.io import FortranFile

#import PSI_parallel_M

RHSofBV = {}
RHSofmJ = {}

os.chdir(litpath3He)

if os.path.isfile(respath + 'kRange.dat') == True:
    os.system('rm ' + respath + 'kRange.dat')

with open(respath + 'kRange.dat', 'wb') as f:
    np.savetxt(f, [anz_phot_e, phot_e_0, phot_e_d], fmt='%f')
    f.seek(NEWLINE_SIZE_IN_BYTES, 2)
    f.truncate()
f.close()

siffux = '_ref'
he_iw, he_rw, he_frgs = retrieve_he3_M(helionpath + 'basis_struct/' +
                                       'INQUA_V18%s' % siffux)
HelBasDimRef = len(sum(sum(he_rw, []), []))

anzStreuBases = len([f for f in glob.glob(litpath3He + 'results/*_BasNR-*')])
finalStatePaths = [
    litpath3He[:-1] + '-%d/' % nB for nB in range(anzStreuBases)
]

for nB in range(anzStreuBases):

    wrkDir = finalStatePaths[nB]
    basisPath = wrkDir + 'basis_struct/'

    siffux = '_fin'
    final_iw, final_rw, final_frgs = retrieve_he3_M(basisPath +
                                                    'INQUA_V18%s' % siffux)
    FinBasDimRef = len(sum(sum(final_rw, []), []))

    with open(respath + 'BareBasDims_%d.dat' % nb, 'wb') as f:
        np.savetxt(f, [HelBasDimRef, FinBasDimRef], fmt='%d')
        f.seek(NEWLINE_SIZE_IN_BYTES, 2)
        f.truncate()
    f.close()

    if 'rhs' in cal:

        os.chdir(wrkDir)

        for streukanal in streukas:

            Jstreu = float(streukanal.split('^')[0])
            Jstreustring = '%s' % str(Jstreu)[:3]

            mLmJl, mLrange, mJlrange = non_zero_couplings(
                multipolarity, J0, Jstreu)

            HelBasDim = sum([
                len(ln.split()[1:])
                for ln in open(helionpath +
                               'basis_struct/SLITbas_full_J%s_%s.dat' %
                               (J0, boundstatekanal))
            ])

            lfrags = []
            sfrags = []
            #for lcfg in range(len(channels[boundstatekanal])):
            #    sfrags = sfrags + channels[boundstatekanal][lcfg][1]
            #    for scfg in channels[boundstatekanal][lcfg][1]:
            #        lfrags = lfrags + [channels[boundstatekanal][lcfg][0]]
            fragfile = [
                ln for ln in open(helionpath +
                                  'basis_struct/Sfrags_LIT_J%s_%s.dat' %
                                  (J0, boundstatekanal))
            ]

            lfrags = [fr.split(' ')[1].strip() for fr in fragfile]
            sfrags = [fr.split(' ')[0] for fr in fragfile]
            # read widths and frags of the LIT basis as determined via
            # v18uix_LITbasis.py
            fragfile = [
                ln
                for ln in open(wrkDir + 'basis_struct/Sfrags_LIT_J%s_%s.dat' %
                               (Jstreustring, streukanal))
            ]

            lfrags2 = [fr.split(' ')[1].strip() for fr in fragfile]
            sfrags2 = [fr.split(' ')[0] for fr in fragfile]
            intwLIT = [
                np.array(ln.split()).astype(float).tolist()
                for ln in open(wrkDir + 'basis_struct/Sintw3heLIT_J%s_%s.dat' %
                               (Jstreustring, streukanal))
            ]
            relwLIT = [
                np.array(ln.split()).astype(float).tolist()
                for ln in open(wrkDir + 'basis_struct/Srelw3heLIT_J%s_%s.dat' %
                               (Jstreustring, streukanal))
            ]

            if 'dbg' in cal:
                print(
                    '\n3He components (full) + LIT-basis components (bare):\n',
                    len(lfrags))
                print(sfrags)
                print('\nLIT-basis components (full):\n', len(lfrags2))
                print(sfrags2)

            if 'rhs_lu-ob-qua' in cal:

                he_iw_2 = [
                    np.array(ln.split()).astype(float).tolist()
                    for ln in open(helionpath +
                                   'basis_struct/Sintw3heLIT_J%s_%s.dat' %
                                   (Jstreustring, boundstatekanal))
                ]

                he_rw_2 = [
                    np.array(ln.split()).astype(float).tolist()
                    for ln in open(helionpath +
                                   'basis_struct/Srelw3heLIT_J%s_%s.dat' %
                                   (Jstreustring, boundstatekanal))
                ]

                for lit_zerl in range(len(lfrags2)):

                    if os.path.isdir(wrkDir + 'tmp_%d' % lit_zerl) == False:
                        os.mkdir(wrkDir + 'tmp_%d' % lit_zerl)
                    os.chdir(wrkDir + 'tmp_%d' % lit_zerl)

                    for file in os.listdir(os.getcwd()):
                        if fnmatch.fnmatch(file, '*J%s*.log' % Jstreu):
                            if 'dbg' in cal:
                                print('removing old <*.log> files.')
                            os.system('rm *.log')
                            break

                    rwtttmp = he_rw + [
                        relwLIT[sum([len(fgg) for fgg in intwLIT[:lit_zerl]]):
                                sum([len(fgg) for fgg in intwLIT[:lit_zerl]]) +
                                len(intwLIT[lit_zerl])]
                    ]
                    lit_3inqua_M(intwi=he_iw + [intwLIT[lit_zerl]],
                                 relwi=rwtttmp,
                                 anzo=11,
                                 LREG='  1  0  0  0  0  0  0  0  0  1  1',
                                 outfile=wrkDir + 'tmp_%d/INQUA' % (lit_zerl))
                    lit_3inlu(mul=multipolarity,
                              frag=lfrags + [lfrags2[lit_zerl]],
                              fn=wrkDir + 'tmp_%d/INLU' % (lit_zerl))
                    lit_3inob(fr=sfrags + [sfrags2[lit_zerl]],
                              fn=wrkDir + 'tmp_%d/INOB' % (lit_zerl))

            leftpar = int(1 + 0.5 *
                          (1 + (-1)**(int(channels[streukanal][0][0][0]) +
                                      int(channels[streukanal][0][0][1]))))

            def cal_rhs_lu_ob_qua(para, procnbr):

                slave_pit = wrkDir + 'tmp_%d' % para
                cmdlu = BINLITpath + 'juelma.exe'
                cmdob = BINLITpath + 'jobelma.exe'
                cmdqu = BINLITpath + 'jquelma.exe'
                print('%s in %s' % (cmdlu, slave_pit))
                plu = subprocess.Popen(shlex.split(cmdlu),
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       cwd=slave_pit)
                out, err = plu.communicate()
                print('process = %d-1 : luise exits.' % para)

                print('%s in %s' % (cmdob, slave_pit))
                pob = subprocess.Popen(shlex.split(cmdob),
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       cwd=slave_pit)
                out, err = pob.communicate()
                print('process = %d-1 : ober exits.' % para)

                print('%s in %s' % (cmdqu, slave_pit))
                pqu = subprocess.Popen(shlex.split(cmdqu),
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       cwd=slave_pit)
                out, err = pqu.communicate()
                print('process = %d-1 : qual exits.' % para)

            def cal_rhs_end(para, procnbr):

                slave_pit = wrkDir + 'tmp_%d/' % para[3]

                inenf = 'inenlit%d-%d_J%3.1f_mJ%3.1f-mL%d.log' % (
                    para[1], para[2], Jstreu, para[0][1], para[0][0])
                outfseli = 'endlit%d-%d_J%3.1f_mJ%3.1f-mL%d.log' % (
                    para[1],
                    para[2],
                    Jstreu,
                    para[0][1],
                    para[0][0],
                )

                # rhs matrix (LMJ0m-M|Jm)*<J_lit m|LM|J0 m-M>
                #  <component>_S_<J>_<mJ>_<M>
                outfsbare = '%s_S_%s_%s_%s.lit' % (
                    str(para[3]).replace('.', ''),
                    #para[1],
                    #para[2],
                    str(Jstreu).replace('.', ''),
                    str(para[0][1]).replace('.', ''),
                    str(para[0][0]).replace('.', ''),
                )

                lit_3inen_bare(MREG='  1  0  0  0  0  0  0  0  0  1  1',
                               JWSL=Jstreu,
                               JWSLM=para[0][1],
                               MULM2=para[0][0],
                               JWSR=J0,
                               outfile=slave_pit + inenf)

                cmdend = BINLITpath + 'jenelmas.exe %s %s %s' % (
                    inenf, outfseli, outfsbare)

                pend = subprocess.Popen(shlex.split(cmdend),
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        cwd=slave_pit)
                out, err = pend.communicate()

                return (out, err)

            parameter_set_lu_ob_qua = range(len(lfrags2))
            parameter_set_end = []

            #        wfn = wrkDir + 'basis_struct/LITbas_full_J%s_%s.dat' % (
            #            Jstreustring, streukanal)
            #
            #        print('[...] reading BV-rw tupel from %s' % wfn)
            #        litbas = [np.loadtxt(wfn).astype(int)[0]]
            #        print(litbas)
            for lit_zerl in range(len(lfrags2)):
                bsbv = sum([len(b) for b in he_iw])
                parameter_set = []
                bvrange = range(
                    sum([len(z) for z in intwLIT[:lit_zerl]]) + 1,
                    sum([len(z) for z in intwLIT[:lit_zerl]]) + 1 +
                    len(intwLIT[lit_zerl]))

                #            litbas3 = []
                #            for bv in filter(lambda x: (x[0] in bvrange), litbas):
                #                litbas3.append([int(bv[0] - (bvrange[0] - 1) + bsbv), bv[1]])
                #
                #            with open(
                #                    wrkDir + 'tmp_%d/LITbas_full_J%s.dat' %
                #                (lit_zerl, Jstreustring), 'wb') as f:
                #                np.savetxt(f, [[jj[0], jj[1]] for jj in litbas3], fmt='%d')
                #                #f.seek(NEWLINE_SIZE_IN_BYTES, 2)
                #                #f.truncate()
                #            f.close()

                for mM in mLmJl:
                    parameter_set.append([mM, 1, 1, lit_zerl])
                    #for bv in filter(lambda x: (x[0] in bvrange), litbas):
                    #    parameter_set.append([
                    #        mM,
                    #        int(bv[0] - (bvrange[0] - 1) + bsbv), bv[1], lit_zerl
                    #    ])

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
                    os.system('cp ' + wrkDir + 'tmp_%d/QUAOUT ' % lit_zerl +
                              wrkDir + 'tmp_%d/QUAOUT_J%3.1f' %
                              (lit_zerl, Jstreu))

            if 'rhs-end' in cal:

                for lit_zerl in range(len(lfrags2)):
                    print('(J=%s)  werkle in %d' % (Jstreu, lit_zerl))
                    try:
                        os.system('cp ' + wrkDir + 'tmp_%d/QUAOUT_J%3.1f ' %
                                  (lit_zerl, Jstreu) + wrkDir +
                                  'tmp_%d/QUAOUT' % (lit_zerl))
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

                os.system('mv ' + wrkDir + 'tmp_*/*_S_* ' + respath)

        if 'rhs-couple' in cal:
            os.chdir(respath)
            rhs = []
            #print('commencing coupling...',HelBasDim)
            for nch in range(len(streukas)):
                streukanal = streukas[nch]
                In = float(streukanal.split('^')[0])
                Jstreustring = '%s' % str(In)[:3]
                fragfile = [
                    ln for ln in open(wrkDir +
                                      'basis_struct/Sfrags_LIT_J%s_%s.dat' %
                                      (Jstreustring, streukanal))
                ]

                lfrags2 = [fr.split(' ')[1].strip() for fr in fragfile]

                mLmJl, mLrange, mJlrange = non_zero_couplings(
                    multipolarity, J0, In)
                for mJ in mJlrange:
                    firstmJ = True
                    for mL in mLrange:

                        Ins = str(In).replace('.', '').ljust(2, '0')
                        mLs = str(mL).replace('.', '').ljust(
                            3, '0') if mL < 0 else str(mL).replace(
                                '.', '').ljust(2, '0')
                        mJs = str(mJ).replace('.', '').ljust(
                            3, '0') if mJ < 0 else str(mJ).replace(
                                '.', '').ljust(2, '0')
                        clebsch = float(
                            CG(multipolarity, mL, J0, mJ - mL, In, mJ).doit())

                        if np.abs(clebsch) != 0:
                            print('(%d,%d;%s,%s|%s,%s) = %f' %
                                  (multipolarity, mL, str(J0), str(mJ - mL),
                                   str(In), str(mJ), clebsch))

                            rhstmp = []
                            for lit_zerl in range(len(lfrags2)):

                                fna = "%d_S_%s_%s_%s.lit" % (lit_zerl, Ins,
                                                             mJs, mLs)

                                kompo_vects_bare = [f for f in glob.glob(fna)]

                                if ((kompo_vects_bare == []) &
                                    (np.abs(clebsch) > 0)):
                                    print(
                                        'RHS component missing: Z,In,MIn,ML:%d,%d,%d,%d'
                                        % (lit_zerl, In, mJ, mL))
                                    print('Clebsch = ', clebsch)
                                    print('file <%s> not found.' % fna)

                                fortranIn = FortranFile(
                                    kompo_vects_bare[0], 'r').read_reals(float)

                                tDim = int(np.sqrt(np.shape(fortranIn)[0]))
                                OutBasDimFr = int(tDim - HelBasDimRef)
                                #print(
                                #    'processing final fragment: %s\ndim(he_bare) = %d ; dim(fin) = %d ; dim(total) = %d'
                                #    % (fna, HelBasDimRef, OutBasDimFr, tDim))

                                subIndices = [
                                    range((HelBasDimRef + ni) * tDim,
                                          (HelBasDimRef + ni) * tDim +
                                          HelBasDimRef)
                                    for ni in range(OutBasDimFr)
                                ]

                                test = np.take(fortranIn, subIndices)
                                test = np.reshape(test, (1, -1))

                                rhstmp = np.concatenate((rhstmp, test[0]))

                            if firstmJ == True:
                                rhsInMIn = clebsch * rhstmp
                                firstmJ = False

                            else:
                                temp = clebsch * rhstmp
                                rhsInMIn = rhsInMIn + temp

                    print('%s -- %s' % (str(In), str(mJ)))
                    outstr = "InMIn_%s_%s_BasNR-%d.%s" % (str(In), str(mJ), nB,
                                                          dt)
                    fortranOut = open(outstr, 'wb+')
                    rhsInMInF = np.asfortranarray(rhsInMIn, dt)
                    rhsInMInF.tofile(fortranOut)
                    fortranOut.close()

    if 'lhs' in cal:

        print('(ii)    calculating norm/ham in scattering-channel basis')

        if os.path.isdir(wrkDir + 'lit_bas_lhs/') == False:
            os.mkdir(wrkDir + 'lit_bas_lhs/')
        os.chdir(wrkDir + 'lit_bas_lhs/')

        for streukanal in streukas:

            Jstreu = float(streukanal.split('^')[0])
            Jstreustring = '%s' % str(Jstreu)[:3]

            fragfile = [
                ln
                for ln in open(wrkDir + 'basis_struct/Sfrags_LIT_J%s_%s.dat' %
                               (Jstreustring, streukanal))
            ]

            lfrags = [fr.split(' ')[1].strip() for fr in fragfile]
            sfrags = [fr.split(' ')[0] for fr in fragfile]

            intwLIT = [
                np.array(ln.split()).astype(float)
                for ln in open(wrkDir + 'basis_struct/Sintw3heLIT_J%s_%s.dat' %
                               (Jstreustring, streukanal))
            ]

            anzLITbv = sum([len(frgm) for frgm in intwLIT])

            if (len([len(frgm) for frgm in intwLIT]) != (len(fragfile))):
                print('LIT-basis fragments inconcistent!',
                      len([len(frgm) for frgm in intwLIT]), (len(fragfile)))
                exit()

            relwLIT = [
                np.array(ln.split()).astype(float)
                for ln in open(wrkDir + 'basis_struct/Srelw3heLIT_J%s_%s.dat' %
                               (Jstreustring, streukanal))
            ]

            if 'dbg' in cal:
                print(lfrags, sfrags)
                #print(intwLIT, relwLIT)

            Jstreu = float(streukanal.split('^')[0])

            if 'lhs_lu-ob-qua' in cal:

                os.chdir(pathbase + '/data/eob/')
                n3_inob([
                    'he_no1',
                    'he_no1y',
                    'he_no2',
                    'he_no2y',
                    'he_no3',
                    'he_no3y',
                    'he_no5',
                    'he_no5y',
                    'he_no6',
                    'he_no6y',
                ],
                        8,
                        fn='INOB',
                        indep=+1)
                os.system(BINBDGpath + 'KOBER.exe')

                os.chdir(pathbase + '/data/elu/')
                n3_inlu(8,
                        fn='INLUCN',
                        fr=[
                            '000',
                            '202',
                            '022',
                            '110',
                            '101',
                            '011',
                            '111',
                            '112',
                            '211',
                            '121',
                            '122',
                            '212',
                            '222',
                            '221',
                            '220',
                        ],
                        indep=+1)
                os.system(BINBDGpath + 'LUDW_CN.exe')

                if tnni == 11:
                    os.chdir(pathbase + '/data/eob-tni/')
                    n3_inob([
                        'he_no1',
                        'he_no1y',
                        'he_no2',
                        'he_no2y',
                        'he_no3',
                        'he_no3y',
                        'he_no5',
                        'he_no5y',
                        'he_no6',
                        'he_no6y',
                    ],
                            15,
                            fn='INOB',
                            indep=+1)
                    os.system(BINBDGpath + 'DROBER.exe')

                    os.chdir(pathbase + '/data/elu-tni/')
                    n3_inlu(8,
                            fn='INLU',
                            fr=[
                                '000',
                                '202',
                                '022',
                                '110',
                                '101',
                                '011',
                                '111',
                                '112',
                                '211',
                                '121',
                                '122',
                                '212',
                                '222',
                                '221',
                                '220',
                            ],
                            indep=+1)
                    os.system(BINBDGpath + 'DRLUD.exe')

                os.chdir(wrkDir + 'lit_bas_lhs/')

                n3_inlu(8, fn='INLUCN', fr=lfrags, indep=-0)
                os.system(BINBDGpath + 'LUDW_CN.exe')
                n3_inob(sfrags, 8, fn='INOB', indep=-0)
                os.system(BINBDGpath + 'KOBER.exe')

                if tnni == 11:
                    n3_inlu(8, fn='INLU', fr=lfrags, indep=-0)
                    os.system(BINBDGpath + 'DRLUD.exe')
                    n3_inob(sfrags, 15, fn='INOB', indep=-0)
                    os.system(BINBDGpath + 'DROBER.exe')

                insam(len(lfrags))

                rwtttmp = []
                for zerle in range(len(lfrags)):
                    rwtttmp.append(
                        relwLIT[sum([len(fgg) for fgg in intwLIT[:zerle]]):
                                sum([len(fgg) for fgg in intwLIT[:zerle]]) +
                                len(intwLIT[zerle])])
                relwLIT = rwtttmp

                he3inquaBS(intwi=intwLIT, relwi=relwLIT, potf=potnn)

                os.system('cp INQUA_M INQUA_M%s' % boundstatekanal)

                parallel_mod_of_3inqua(lfrags,
                                       sfrags,
                                       infile='INQUA_M',
                                       outfile='INQUA_M')

                subprocess.run([
                    'mpirun', '-np',
                    '%d' % anzproc, BINBDGpath + 'V18_PAR/mpi_quaf_v7'
                ])

                subprocess.run([BINBDGpath + 'V18_PAR/sammel'])

                if tnni == 11:
                    he3inquaBS(intwi=intwLIT, relwi=relwLIT, potf=potnnn)

                    parallel_mod_of_3inqua(lfrags,
                                           sfrags,
                                           infile='INQUA_M',
                                           outfile='INQUA_M',
                                           tni=1)

                    subprocess.run([
                        'mpirun', '-np',
                        '%d' % anzproc, BINBDGpath + 'UIX_PAR/mpi_drqua_v7'
                    ])

                    subprocess.run([BINBDGpath + 'UIX_PAR/SAMMEL-uix'])

            #litbas = np.loadtxt(wrkDir +
            #                    'basis_struct/LITbas_full_J%s_%s.dat' %
            #                    (Jstreustring, streukanal)).astype(int)
            litbas2 = [[
                int(line.split()[0]), [int(rwn) for rwn in line.split()[1:]]
            ] for line in open(wrkDir +
                               'basis_struct/SLITbas_full_J%s_%s.dat' %
                               (Jstreustring, streukanal))]
            lb = []
            for bv in litbas2:
                for rw in range(len(bv[1])):
                    lb.append([bv[0], int(bv[1][rw])])

            #        litbas = [bv for bv in np.unique(litbas, axis=0) if bv[1] != 0]
            litbas = np.array(lb).astype(int)
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

                mypath = wrkDir + 'tmp_%d/' % lit_zerl
                lit_zerl += 1
                n3_inen_rhs(litbas,
                            Jstreu,
                            costr,
                            np.ones(len(relwLIT[0])),
                            fn='INEN',
                            pari=0,
                            nzop=zop,
                            tni=tnni,
                            anzb=anzbtmp)
                os.system('cp INEN ' + mypath + 'inen-lit-%s_1-%d' %
                          (streukanal, anzbtmp))
                #if anzbtmp==anzbs[-1]:
                #    subprocess.run([BINBDGpath + 'TDR2END_NORMAL.exe'])
                #    os.system('cp %s/MATOUT ' % (wrkDir + 'lit_bas_lhs/') + respath
                #          + 'norm-ham-litME-%s_1-%d' % (streukanal, anzbtmp))
                #os.system('cp ' + v18uixpath + 'mat_* ' + respath)

            #plotHspec(Jstreustring)

    if 'couple' in cal:

        os.system('cp ' + helionpath + 'E0.dat ' + wrkDir)

        for streukanal in streukas:

            Jstreu = float(streukanal.split('^')[0])
            Jstreustring = '%s' % str(Jstreu)[:3]

            intwLIT = [
                np.array(ln.split()).astype(float).tolist()
                for ln in open(wrkDir + 'basis_struct/Sintw3heLIT_J%s_%s.dat' %
                               (Jstreustring, streukanal))
            ]

            bv_pro_zerl = [len(zset) for zset in intwLIT]
            first_bv_pro_zerl = [
                sum(bv_pro_zerl[:n]) + 1 for n in range(len(bv_pro_zerl))
            ]

            relwLIT = [
                np.array(ln.split()).astype(float).tolist()
                for ln in open(wrkDir + 'basis_struct/Srelw3heLIT_J%s_%s.dat' %
                               (Jstreustring, streukanal))
            ]

            fragfile = [
                ln
                for ln in open(wrkDir + 'basis_struct/Sfrags_LIT_J%s_%s.dat' %
                               (Jstreustring, streukanal))
            ]

            lfrags = [fr.split(' ')[1].strip() for fr in fragfile]

            litbas = np.loadtxt(wrkDir +
                                'basis_struct/SLITbas_full_J%s_%s.dat' %
                                (Jstreustring, streukanal)).astype(int)
            litbas = [bv for bv in np.unique(litbas, axis=0) if bv[1] != 0]

            print('He3 structure: anz[int,rel]weiten:', np.sum(he_frgs,
                                                               axis=0))

            bv_offset = 0

            RHSofBV[streukanal] = {}
            RHSofmJ[streukanal] = {}

            for lit_zerl in range(len(lfrags)):

                # read uncoupled source ME's
                os.chdir(wrkDir + 'tmp_%d' % lit_zerl)

                litbas = np.loadtxt(wrkDir +
                                    'basis_struct/SLITbas_full_J%s_%s.dat' %
                                    (Jstreustring, streukanal)).astype(int)
                litbas = [bv for bv in np.unique(litbas, axis=0) if bv[1] != 0]

                print('werkle in fragment %s' % os.getcwd())
                tmp_iw, tmp_rw, tmp_frgs = retrieve_he3_widths('INQUA')

                he_bv = np.sum(he_frgs, axis=0)[0]
                print('bv in He3 : %d' % he_bv)
                bv_in_zerl = tmp_frgs[-1][0]

                print('bv in zerl : %d' % bv_in_zerl)

                zerl_bas = [
                    bv for bv in litbas if first_bv_pro_zerl[lit_zerl] <= bv[0]
                    < first_bv_pro_zerl[lit_zerl + 1]
                ] if (lit_zerl < len(lfrags) - 1) else [
                    bv for bv in litbas if first_bv_pro_zerl[lit_zerl] <= bv[0]
                ]

                zerl_bas = [[
                    bv[0] + 1 + he_bv - first_bv_pro_zerl[lit_zerl], bv[1]
                ] for bv in zerl_bas]

                photEn = []

                read_uncoupled_source(RHSofBV[streukanal],
                                      photEn,
                                      streukanal,
                                      basisSET=zerl_bas,
                                      firstbv=bv_offset)
                #RHSofBV[streukanal] = read_norm(streukanal, basisSET=litbas)
                #photEn = MeVfm * np.array(
                #    [phot_e_0 + en * phot_e_d for en in range(anz_phot_e)])
                # couple incoming state with photon multipole to Jlit
                couple_source(RHSofmJ[streukanal],
                              streukanal,
                              RHSofBV[streukanal],
                              basisSET=zerl_bas,
                              firstbv=bv_offset)

                bv_offset += len(zerl_bas)

    os.chdir(wrkDir)

    #os.system('rm -rf ./tmp*')
    os.system('find . -name \"T*OUT.*\" -print0 | xargs -0 rm')