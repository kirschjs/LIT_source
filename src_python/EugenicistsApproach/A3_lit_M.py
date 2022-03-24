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

os.chdir(bkpdir)

respath = bkpdir + suffi + 'results/'
helionpath = bkpdir + suffi + 'he3/'

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

anzStreuBases = len(
    [f for f in glob.glob(respath + 'mat_%s_BasNR-*' % streukas[0])])

finalStatePaths = [
    bkpdir + suffi[:-1] + '-%d/' % nB for nB in range(anzStreuBases)
]

for nB in range(anzStreuBases):

    wrkDir = finalStatePaths[nB]
    basisPath = wrkDir + 'basis_struct/'

    siffux = '_fin'
    final_iw, final_rw, final_frgs = retrieve_he3_M(basisPath +
                                                    'INQUA_V18%s' % siffux)
    FinBasDimRef = len(sum(sum(final_rw, []), []))

    with open(respath + 'BareBasDims_%d.dat' % nB, 'wb') as f:
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
                               'basis_struct/SLITbas_full_%s.dat' %
                               boundstatekanal)
            ])

            lfrags = []
            sfrags = []
            #for lcfg in range(len(channels[boundstatekanal])):
            #    sfrags = sfrags + channels[boundstatekanal][lcfg][1]
            #    for scfg in channels[boundstatekanal][lcfg][1]:
            #        lfrags = lfrags + [channels[boundstatekanal][lcfg][0]]
            fragfile = [
                ln
                for ln in open(helionpath + 'basis_struct/Sfrags_LIT_%s.dat' %
                               boundstatekanal)
            ]

            lfrags = [fr.split(' ')[1].strip() for fr in fragfile]
            sfrags = [fr.split(' ')[0] for fr in fragfile]
            # read widths and frags of the LIT basis as determined via
            # v18uix_LITbasis.py
            fragfile = [
                ln
                for ln in open(wrkDir +
                               'basis_struct/Sfrags_LIT_%s.dat' % streukanal)
            ]

            lfrags2 = [fr.split(' ')[1].strip() for fr in fragfile]
            sfrags2 = [fr.split(' ')[0] for fr in fragfile]
            intwLIT = [
                np.array(ln.split()).astype(float).tolist()
                for ln in open(wrkDir +
                               'basis_struct/Sintw3heLIT_%s.dat' % streukanal)
            ]
            relwLIT = [
                np.array(ln.split()).astype(float).tolist()
                for ln in open(wrkDir +
                               'basis_struct/Srelw3heLIT_%s.dat' % streukanal)
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
                                   'basis_struct/Sintw3heLIT_%s.dat' %
                                   boundstatekanal)
                ]

                he_rw_2 = [
                    np.array(ln.split()).astype(float).tolist()
                    for ln in open(helionpath +
                                   'basis_struct/Srelw3heLIT_%s.dat' %
                                   boundstatekanal)
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
                cmdlu = BINLITpath + 'juelmanoo.exe'
                cmdob = BINLITpath + 'jobelmanoo.exe'
                cmdqu = BINLITpath + 'jquelmanoo.exe'
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

                cmdend = BINLITpath + 'jenelmasnoo.exe %s %s %s' % (
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
                pool = ThreadPool(min(MaxProc, len(parameter_set)))
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
                    pool = ThreadPool(min(MaxProc, len(parameter_set)))
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

                    #tmps = wrkDir + 'tmp_%d' % lit_zerl
                    #subprocess.call('mv %s %s' % (tmps, respath), shell=True)
                    #os.chdir(wrkDir)
                    #subprocess.call('rm  -rf %s' % tmps, shell=True)

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
                    ln
                    for ln in open(wrkDir + 'basis_struct/Sfrags_LIT_%s.dat' %
                                   streukanal)
                ]

                lfrags2 = [fr.split(' ')[1].strip() for fr in fragfile]

                mLmJl, mLrange, mJlrange = non_zero_couplings(
                    multipolarity, J0, In)
                print(multipolarity, J0, In, ':', mJlrange)
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

    os.chdir(bkpdir)
    subprocess.call('rm  -rf %s' % wrkDir, shell=True)
    #os.system('find . -name \"T*OUT.*\" -print0 | xargs -0 rm')

resdest = bkpdir + '/latestresults'
resdestbkp = resdest + '%s' % (datetime.now().strftime("%d-%b-%Y--%H-%M-%S"))
if os.path.isdir(resdest) == True:
    os.system('mv %s %s' % (resdest, resdestbkp))

os.system('cp -r %s %s' % (respath[:-1], resdest))
print('\n\nECCE! Results copied from\n%s\nto\n%s' % (respath[:-1], resdest))