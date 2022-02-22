from bridgeA3 import *
from genetic_width_growth import *

import glob, copy
import shlex
import multiprocessing
from multiprocessing.pool import ThreadPool
from itertools import permutations, product

from PSI_parallel_M import span_initial_basis

os.chdir(litpath3He)

dbg = False

bastypes = [boundstatekanal]  #streukas  #+

# > optimize the various basis types, e.g., in case of the npp system:
# > helion ground state, final J=1/2- and J=3/2- states
for bastype in bastypes:

    Jay = float(bastype.split('^')[0][-3:])
    Jaystring = '%s' % str(Jay)[:3]

    # number of final-state bases which are grown with the above-set criteria
    anzStreuBases = 1

    costr = ''
    zop = 31 if tnni == 11 else 14
    for nn in range(1, zop):
        if bastype == boundstatekanal:
            cf = 1.0 if (1 <= nn <= 28) else 0.0
        else:
            cf = 1.0 if (1 <= nn <= 28) else 0.0
        costr += '%12.7f' % cf if (nn % 7 != 0) else '%12.7f\n' % cf

    if bastype == boundstatekanal:
        if os.path.isdir(helionpath) != False:
            print('an optimized initial state is already present')
            #continue
        else:
            os.mkdir(helionpath)
            os.mkdir(helionpath + 'basis_struct/')

    else:
        finalStatePaths = [
            litpath3He[:-1] + '-%d/' % nB for nB in range(anzStreuBases)
        ]
        for finalStatePath in finalStatePaths:
            if os.path.isdir(finalStatePath) == True:
                os.system('rm -rf ' + finalStatePath)

            os.mkdir(finalStatePath)
            os.mkdir(finalStatePath + 'basis_struct/')

    # numerical stability
    minDiffwidthsINT = 10**-2
    minDiffwidthsREL = 10**-3
    maxParLen = 18
    pwpurge = 0
    purgeStr = ['Condition number', 'pulchritude',
                'Ground-state energy'][pwpurge]

    # evolution criteria
    minCond = 10**-6
    denseEVinterval = [-10., 150.0]
    removalGainFactor = 1.5
    maxOnPurge = 43
    maxOnTrail = 10**2
    muta_initial = 0.15
    pwopt = 1

    CgfCycles = 3
    # nRaces := |i|
    nRaces = 5 if bastype == boundstatekanal else 5

    cradleCapacity = 60

    # > nState > produce/optimize/grow multiple bases with pseudo-random initial seeds
    for nB in range(anzStreuBases):

        wrkDir = helionpath if bastype == boundstatekanal else finalStatePaths[
            nB]
        basisPath = wrkDir + 'basis_struct/'

        os.chdir(wrkDir)

        seedMat = span_initial_basis(
            basisType=bastype,
            ini_grid_bounds=[0.06, 6.25, 0.04, 7.5, 0.005, 8.25, 0.001, 7.5],
            ini_dims=[5, 7, 12, 4],
            coefstr=costr,
            anzOp=zop)

        dim = int(np.sqrt(len(seedMat) * 0.5))

        # read Norm and Hamilton matrices
        normat = np.reshape(
            np.array(seedMat[:dim**2]).astype(float), (dim, dim))
        hammat = np.reshape(
            np.array(seedMat[dim**2:]).astype(float), (dim, dim))
        # diagonalize normalized norm (using "eigh(ermitian)" to speed-up the computation)
        ewN, evN = eigh(normat)
        idx = ewN.argsort()[::-1]
        ewN = [eww for eww in ewN[idx]]
        evN = evN[:, idx]

        try:
            ewH, evH = eigh(hammat, normat)
            idx = ewH.argsort()[::-1]
            ewH = [eww for eww in ewH[idx]]
            evH = evH[:, idx]

        except:
            print(
                'failed to solve generalized eigenvalue problem (norm ev\'s < 0 ?) in chan %s for basis set (%d)'
                % (bastype, nB))
            exit()

        anzSigEV = len([
            bvv for bvv in ewH if denseEVinterval[0] < bvv < denseEVinterval[1]
        ])
        gsEnergy = ewH[-1]
        basCond = np.min(np.abs(ewN)) / np.max(np.abs(ewN))
        attractiveness = loveliness(gsEnergy, basCond, anzSigEV, minCond)

        print(
            '\n> basType %s > basSet %d/%d: seed basis: E0 = %f   cond=|Emin|/|Emax| = %e'
            % (bastype, nB + 1, anzStreuBases, gsEnergy, basCond))

        cfgs = [
            con.split()
            for con in open(basisPath + 'frags_LIT_%s.dat' % bastype)
        ]
        origCFGs = copy.deepcopy(cfgs)

        intwLIT = [
            np.array(ln.split()).astype(float).tolist()
            for ln in open(basisPath + 'intw3heLIT_%s.dat' % bastype)
        ]
        relwLIT = [
            np.array(ln.split()).astype(float).tolist()
            for ln in open(basisPath + 'relw3heLIT_%s.dat' % bastype)
        ]
        rws = []
        rw0 = 0
        for cfg in range(len(intwLIT)):
            rws.append([])
            for bv in range(len(intwLIT[cfg])):
                rws[-1].append(relwLIT[bv + rw0])
            rw0 += len(intwLIT[cfg])

        initialCiv = [cfgs, intwLIT, rws, []]
        nbv = 0
        for cfg in range(len(initialCiv[0])):
            nbvc = 0
            for bv in initialCiv[1][cfg]:
                nbv += 1
                nbvc += 1
                initialCiv[3] += [[
                    nbv,
                    np.array(range(
                        1, 1 + len(initialCiv[2][cfg][nbvc - 1]))).tolist()
                ]]

        # > nState > nBasis > stabilize the seed basis
        goPurge = True

        D0 = initialCiv[3]

        # purge just entire bv sets with identical internal width
        print('Stratifying the initial seed -- criterion: %s' % purgeStr)
        while goPurge:
            newpopList = []
            goPurge = False
            ParaSets = []

            ParaSets.append([
                D0, Jay, costr, zop, 10, [0, 0], BINBDGpath, minCond,
                denseEVinterval
            ])

            for bvTrail in D0:

                bvID = [int(bvTrail[0]), int(''.join(map(str, bvTrail[1])))]
                cpy = copy.deepcopy(D0)
                cpy.remove(bvTrail)

                ParaSets.append([
                    cpy, Jay, costr, zop, 10, bvID, BINBDGpath, minCond,
                    denseEVinterval
                ])

            tst = np.random.choice(np.arange(len(ParaSets)),
                                   size=min(maxOnPurge, len(ParaSets)),
                                   replace=False)
            if not 0 in tst:
                tst = tst.tolist() + [0]
            if maxOnPurge < len(ParaSets):
                tkkg = [ParaSets[t] for t in tst]
                ParaSets = tkkg
            # x) the parallel environment is set up in sets(chunks) of bases
            #    in order to limit the number of files open simultaneously
            split_points = [
                n * maxParLen
                for n in range(1 + int(len(ParaSets) / maxParLen))
            ] + [len(ParaSets) + 1024]

            Parchunks = [
                ParaSets[split_points[i]:split_points[i + 1]]
                for i in range(len(split_points) - 1)
            ]

            cand_list = []

            print(
                'checking each vector in %d-dim basis on its effect on the stability,'
                % len(ParaSets))
            for chunk in Parchunks:
                pool = ThreadPool(max(min(MaxProc, len(ParaSets)), 2))
                jobs = []
                for procnbr in range(len(chunk)):
                    recv_end, send_end = multiprocessing.Pipe(False)
                    pars = chunk[procnbr]
                    p = multiprocessing.Process(target=endmat,
                                                args=(pars, send_end))
                    jobs.append(p)
                    cand_list.append(recv_end)
                    p.start()
                for proc in jobs:
                    proc.join()
            cand_ladder = [x.recv() for x in cand_list]

            # ranking following condition-number (0) or quality (1)
            cand_ladder.sort(key=lambda tup: np.abs(tup[pwpurge]))
            reff = [[bvm[0], bvm[1], bvm[2]] for bvm in cand_ladder
                    if bvm[3] == [0, 0]][0]

            for cand in cand_ladder[-3:]:
                print(cand[:4])

            print(
                '\n> basType %s > basSet %d/%d: purged seed: E0 = %f   cond=|Emin|/|Emax| = %e'
                % (bastype, nB + 1, anzStreuBases, cand_ladder[-1][2],
                   cand_ladder[-1][0]))

            print('    best =  %2.3e , %2.3e , %2.3e' %
                  (cand_ladder[-1][0], cand_ladder[-1][1], cand_ladder[-1][2]))
            print('    worst=  %2.3e , %2.3e , %2.3e' %
                  (cand_ladder[0][0], cand_ladder[0][1], cand_ladder[0][2]))

            if ((removalGainFactor * np.abs(reff[pwpurge]) < np.abs(
                    cand_ladder[-1][pwpurge])) |
                ((np.abs(cand_ladder[-1][pwpurge]) < minCond) &
                 (reff[0] < minCond))):
                goPurge = True
                D0 = rectify_basis(cand_ladder[-1][4])
                print('1/%d ' % basisDim(D0), end='\n')

        initialCiv[3] = rectify_basis(cand_ladder[-1][4])
        # > nState > nBasis > end of stabilization

        print('\niniCivSeed:\n', initialCiv)

        initialCiv = essentialize_basis(initialCiv, MaxBVsPERcfg=bvma)

        print('\niniCivSeed (strat):\n', initialCiv)

        ma = blunt_ev(initialCiv[0],
                      initialCiv[1],
                      initialCiv[2],
                      initialCiv[3],
                      wrkdir='',
                      nzopt=zop,
                      costring=costr,
                      bin_path=BINBDGpath,
                      mpipath=MPIRUN,
                      einzel_file_path=wrkDir,
                      potNN=potnn,
                      potNNN=potnnn,
                      parall=-1,
                      anzcores=max(2, min(len(initialCiv[0]), MaxProc)),
                      tnni=10,
                      jay=Jay,
                      dia=False)

        ewN, ewH = NormHamDiag(ma)

        parCond = np.min(np.abs(ewN)) / np.max(np.abs(ewN))

        print(
            '>\n>\n>\n> stabilized initial basis: C-nbr = %4.4e E0 = %4.4e\n>\n>\n>'
            % (parCond, ewH[-1]))

        # > nState > nBasis > optimize each orb-ang, spin-iso cfg in a number of cycles
        for nCgfCycle in range(CgfCycles):

            # > nState > nBasis > nCfgCycle > optimize a single angular-momentum configuration, e.g. l1=1,l2=1,L=2,s12=....
            unis = []
            for ncfg in range(len(initialCiv[0])):
                if initialCiv[0][ncfg] in unis:
                    continue
                else:
                    unis.append(initialCiv[0][ncfg])

            for nUcfg in range(len(unis)):

                print(
                    '\n> basType %s > basSet %d/%d > nUcfg %d/%d > : Optimizing cfg: '
                    % (bastype, nB + 1, anzStreuBases, nUcfg + 1, len(unis)),
                    unis[nUcfg])
                # > nState > nBasis > nCfgCycle > nUcfg > allow each cfg to evolve over nRaces
                for nGen in range(nRaces):

                    print(
                        '\n> basType %s > basSet %d/%d > nCfg %d/%d > nGen %d/%d: breeding...'
                        % (bastype, nB + 1, anzStreuBases, nUcfg + 1,
                           len(unis), nGen + 1, nRaces))

                    Ais = copy.deepcopy(initialCiv)

                    bvPerCfg = [len(iws) for iws in Ais[1]]
                    maxBVinCfg = np.cumsum(bvPerCfg)
                    anzBV = sum(bvPerCfg)

                    # find all cfgs with angular-momentum structure = nUcfg
                    ret = [
                        n for n in range(len(initialCiv[0]))
                        if initialCiv[0][n] == unis[nUcfg]
                    ]

                    # parents are all vectors in those blocks

                    parentIWs = sum([Ais[1][nCfg] for nCfg in ret], [])
                    parentRWs = sum([Ais[2][nCfg] for nCfg in ret], [])

                    chiBVs = []
                    # produce the offspring cfg for nCfg

                    # from the iw's of the parent cfg, select mother/father pairs
                    if len(parentIWs) > 1:

                        iwpairs = [
                            ip for ip in list(
                                product(range(len(parentIWs)), repeat=2))
                            if ip[0] != ip[1]
                        ]

                        np.random.shuffle(iwpairs)
                        iwpairs = iwpairs[:int((len(parentIWs)) / 2)]

                    else:
                        iwpairs = [(0, np.nan)]

                    daughtersonI = []
                    Ais[2].append([])

                    nCh = 1

                    while nCh < cradleCapacity:

                        for iws in iwpairs:

                            daughtersonR = []
                            motherI = parentIWs[iws[0]]
                            if np.isnan(iws[1]) == False:
                                fatherI = parentIWs[iws[1]]
                            else:
                                fatherI = motherI * np.random.random()

                            daughtersonI.append(
                                intertwining(motherI,
                                             fatherI,
                                             mutation_rate=muta_initial))

                            for nrw in range(len(parentRWs[iws[0]])):
                                motherR = parentRWs[iws[0]][nrw]
                                if iws[1] >= 0:
                                    fatherR = parentRWs[iws[1]][nrw]
                                else:
                                    fatherR = motherR * np.random.random()

                                daughtersonR.append(
                                    intertwining(motherR,
                                                 fatherR,
                                                 mutation_rate=muta_initial))

                            rw1 = list(np.array(daughtersonR)[:, 1])
                            rw1.sort()
                            rw1 = rw1[::-1]
                            rw2 = list(np.array(daughtersonR)[:, 0])
                            rw2.sort()
                            rw2 = rw2[::-1]

                            Ais[2][-1].append(rw2)
                            Ais[2][-1].append(rw1)

                            Ais[3] = Ais[3] + [[
                                anzBV + nCh,
                                list(range(1, 1 + len(rw1)))
                            ]]
                            Ais[3] = Ais[3] + [[
                                anzBV + nCh + 1,
                                list(range(1, 1 + len(rw2)))
                            ]]

                            nCh += 2

                    dstmp = np.array(daughtersonI).flatten()
                    dstmp.sort()
                    daughtersonI = list(dstmp)[::-1]
                    Ais[1].append(daughtersonI)

                    Ais[0].append(unis[nUcfg])
                    Ais = essentialize_basis(Ais, MaxBVsPERcfg=bvma)

                    bvPerCfg = [len(iws) for iws in Ais[1]]
                    maxBVinCfg = np.cumsum(bvPerCfg)
                    anzBV = sum(bvPerCfg)

                    #bvLowerBnd = 0 if nCfg == 0 else maxBVinCfg[nCfg - 1]
                    #bvUpperBnd = maxBVinCfg[nCfg]

                    # now we have to rate all BV's for this cfg
                    # we do this by removing the BV, and test how well the cfg performs w/o it
                    # we remove the one which does make the least difference
                    #                    chID = 1
                    #                    tmpBas.sort()
                    #                    for bvTrail in [
                    #                            bv for bv in tmpBas
                    #                            if bvLowerBnd < bv[0] <= bvUpperBnd
                    #                    ]:
                    #                        childidentifier = [chID]
                    #                        cpy = copy.deepcopy(tmpBas)
                    #                        cpy.remove(bvTrail)
                    #                        cpy = rectify_basis(cpy)
                    #                        childishParaSets.append([
                    #                            cpy, Jay, costr, zop, 10, childidentifier,
                    #                            BINBDGpath, minCond, denseEVinterval
                    #                        ])
                    #                        chID += 1

                    parBVs = Ais[3][:anzBV - nCh + 1]
                    chBVs = Ais[3][anzBV - nCh + 1:]

                    #print('\nparBVs:\n', parBVs)
                    #print('\nchBVs:\n', chBVs)

                    childishParaSets = []
                    childishParaSets.append([
                        parBVs, Jay, costr, zop, 10, [-1], BINBDGpath, minCond,
                        denseEVinterval
                    ])

                    for bvTrail in chBVs:
                        childidentifier = [bvTrail[0]]
                        cpy = copy.deepcopy(parBVs)
                        cpy.append(bvTrail)
                        cpy = rectify_basis(cpy)
                        cpy.sort()
                        childishParaSets.append([
                            cpy, Jay, costr, zop, 10, childidentifier,
                            BINBDGpath, minCond, denseEVinterval
                        ])

                    # do this only if all children are in basis[3], otherwise their widths are removed!
                    # do it (!) to break frgs which became too long after the children were added

                    #print('\nAis:\n', Ais)
                    Ais = essentialize_basis(Ais, MaxBVsPERcfg=bvma)
                    #print('\nAis (strat):\n', Ais)
                    ma = blunt_ev(Ais[0],
                                  Ais[1],
                                  Ais[2],
                                  parBVs,
                                  wrkdir='',
                                  nzopt=zop,
                                  costring=costr,
                                  bin_path=BINBDGpath,
                                  mpipath=MPIRUN,
                                  einzel_file_path=wrkDir,
                                  potNN=potnn,
                                  potNNN=potnnn,
                                  parall=-1,
                                  anzcores=max(2, min(len(Ais[0]), MaxProc)),
                                  tnni=10,
                                  jay=Jay,
                                  dia=False)

                    ewN, ewH = NormHamDiag(ma)

                    parCond = np.min(np.abs(ewN)) / np.max(np.abs(ewN))

                    if ewH == []:
                        print('old generation already unstable:\n', Ais)
                        exit()

                    print(
                        'reference for the new gen: Dim(all children + parents) = %d) B(GS, parents) = %8.4f  C-number (parents): %4.3e\nfitness (parents)  = '
                        % (basisDim(Ais[3]), ewH[-1], parCond),
                        basQ(ewN, ewH, minCond, denseEVinterval))

                    tst = np.random.choice(np.arange(len(childishParaSets)),
                                           size=min(maxOnTrail,
                                                    len(childishParaSets)),
                                           replace=False)

                    if not 0 in tst:
                        tst = tst.tolist() + [0]

                    if maxOnTrail < len(childishParaSets):
                        tkkg = [childishParaSets[t] for t in tst]
                        childishParaSets = tkkg

                    # x) the parallel environment is set up in sets(chunks) of bases
                    #    in order to limit the number of files open simultaneously
                    split_points = [
                        n * maxParLen
                        for n in range(1 +
                                       int(len(childishParaSets) / maxParLen))
                    ] + [len(childishParaSets) + 1024]

                    Parchunks = [
                        childishParaSets[split_points[i]:split_points[i + 1]]
                        for i in range(len(split_points) - 1)
                    ]

                    #for childishParaSet in childishParaSets:
                    #    print(childishParaSet[0], childishParaSet[1],
                    #          childishParaSet[3])
                    cand_list = []
                    for chunk in Parchunks:
                        pool = ThreadPool(
                            max(min(MaxProc, len(childishParaSets)), 2))
                        jobs = []
                        for procnbr in range(len(chunk)):
                            recv_end, send_end = multiprocessing.Pipe(False)

                            pars = chunk[procnbr]

                            p = multiprocessing.Process(target=endmat,
                                                        args=(pars, send_end))
                            jobs.append(p)
                            cand_list.append(recv_end)
                            p.start()
                        for proc in jobs:
                            proc.join()

                    cand_ladder = [x.recv() for x in cand_list]

                    # ranking following condition-number (0) or quality (1)  or E(GS) (2)
                    cand_ladder.sort(key=lambda tup: np.abs(tup[pwopt]))

                    for cand in cand_ladder[-3:]:
                        print(cand[:4])
                    # optimum when deleting one

                    for cand in cand_ladder[::-1]:
                        if ((cand[3] != [-1]) & (cand[0] > minCond)):
                            parvenue = cand
                            break

                    if parvenue == []:
                        continue

                    print('\nparvenue:\n', parvenue)

                    Ais[3] = parvenue[-1]

                    #print(Ais[0], Ais[3])

                    initialCiv = essentialize_basis(Ais, MaxBVsPERcfg=bvma)

                    #print(initialCiv[0], initialCiv[3])

                    initialCiv = condense_Basis(initialCiv,
                                                unis,
                                                bvsPERcfg=bvma)

                    #print(initialCiv[0], initialCiv[3])

        exit()

        # after having evolved each configuration, stabilize the basis
        # and remove its least siginificant vectors
        ma = blunt_ev(initialCiv[0],
                      initialCiv[1],
                      initialCiv[2],
                      initialCiv[3],
                      wrkdir='',
                      nzopt=zop,
                      costring=costr,
                      bin_path=BINBDGpath,
                      mpipath=MPIRUN,
                      einzel_file_path=wrkDir,
                      potNN=potnn,
                      potNNN=potnnn,
                      parall=-1,
                      anzcores=max(2, min(len(initialCiv[0]), MaxProc)),
                      tnni=10,
                      jay=Jay,
                      dia=False)

        D0 = initialCiv[3]
        # purge just entire bv sets with identical internal width
        print('Stratifying generation (%d) -- criterion: %s' %
              (nGen, purgeStr))

        while goPurge:
            newpopList = []
            goPurge = False
            ParaSets = []
            ParaSets.append([
                D0, Jay, costr, zop, 10, [0, 0], BINBDGpath, minCond,
                denseEVinterval
            ])
            for bvTrail in D0:
                bvID = [int(bvTrail[0]), int(''.join(map(str, bvTrail[1])))]
                cpy = copy.deepcopy(D0)
                cpy.remove(bvTrail)
                ParaSets.append([
                    cpy, Jay, costr, zop, 10, bvID, BINBDGpath, minCond,
                    denseEVinterval
                ])

            tst = np.random.choice(np.arange(len(ParaSets)),
                                   size=min(maxOnPurge, len(ParaSets)),
                                   replace=False)
            if not 0 in tst:
                tst = tst.tolist() + [0]
            if maxOnPurge < len(ParaSets):
                tkkg = [ParaSets[t] for t in tst]
                ParaSets = tkkg
            # x) the parallel environment is set up in sets(chunks) of bases
            #    in order to limit the number of files open simultaneously
            split_points = [
                n * maxParLen
                for n in range(1 + int(len(ParaSets) / maxParLen))
            ] + [len(ParaSets) + 1024]
            Parchunks = [
                ParaSets[split_points[i]:split_points[i + 1]]
                for i in range(len(split_points) - 1)
            ]
            cand_list = []
            print(
                'checking each vector in %d-dim basis on its effect on the stability,'
                % len(ParaSets))
            for chunk in Parchunks:
                pool = ThreadPool(max(min(MaxProc, len(ParaSets)), 2))
                jobs = []
                for procnbr in range(len(chunk)):
                    recv_end, send_end = multiprocessing.Pipe(False)
                    pars = chunk[procnbr]
                    p = multiprocessing.Process(target=endmat,
                                                args=(pars, send_end))
                    jobs.append(p)
                    cand_list.append(recv_end)
                    p.start()
                for proc in jobs:
                    proc.join()
            cand_ladder = [x.recv() for x in cand_list]
            # ranking following condition-number (0) or quality (1)
            cand_ladder.sort(key=lambda tup: np.abs(tup[pwpurge]))
            reff = [[bvm[0], bvm[1], bvm[2]] for bvm in cand_ladder
                    if bvm[3] == [0, 0]][0]
            for cand in cand_ladder[-3:]:
                print(cand[:4])
            print(
                '\n> basType %s > basSet %d/%d: purged seed: E0 = %f   cond=|Emin|/|Emax| = %e'
                % (bastype, nB + 1, anzStreuBases, cand_ladder[-1][2],
                   cand_ladder[-1][0]))
            print('    best =  %2.3e , %2.3e , %2.3e' %
                  (cand_ladder[-1][0], cand_ladder[-1][1], cand_ladder[-1][2]))
            print('    worst=  %2.3e , %2.3e , %2.3e' %
                  (cand_ladder[0][0], cand_ladder[0][1], cand_ladder[0][2]))
            if ((removalGainFactor * np.abs(reff[pwpurge]) < np.abs(
                    cand_ladder[-1][pwpurge])) |
                ((np.abs(cand_ladder[-1][pwpurge]) < minCond) &
                 (reff[0] < minCond))):
                goPurge = True
                D0 = rectify_basis(cand_ladder[-1][4])
                print('1/%d ' % basisDim(D0), end='\n')
        initialCiv[3] = rectify_basis(cand_ladder[-1][4])
        # > nState > nBasis > end of stabilization
        initialCiv = essentialize_basis(initialCiv, MaxBVsPERcfg=bvma)
