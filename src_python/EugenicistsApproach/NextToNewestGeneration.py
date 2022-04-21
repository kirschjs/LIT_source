from bridgeA3 import *
from genetic_width_growth import *

import glob, copy
import shlex
import multiprocessing
from multiprocessing.pool import ThreadPool
from itertools import permutations, product

from PSI_parallel_M import span_initial_basis

if os.path.isdir(litpath3He) == False:
    subprocess.check_call(['mkdir -p', litpath3He])
    subprocess.check_call(['mkdir -p', respath])

with open(respath + 'dtype.dat', 'w') as outf:
    outf.write(dt)

os.chdir(litpath3He)

dbg = False

arglist = sys.argv

if arglist[1:] != []:
    StreuBases = np.arange(int(arglist[1]), int(arglist[2]) + 1)
    anzStreuBases = len(StreuBases)
else:
    anzStreuBases = 1
    StreuBases = np.arange(1, anzStreuBases + 1)

bastypes = [boundstatekanal] + streukas if StreuBases[0] == 1 else streukas

if 1 in StreuBases:
    if os.path.isdir(helionpath) != False:
        print('<ECCE> removing the existing helion folder\n%s.' % helionpath)
        os.system('rm -rf ' + helionpath)
    subprocess.check_call([' -p', helionpath])
    subprocess.check_call(['mkdir -p', helionpath + 'basis_struct/'])

finalStatePaths = [litpath3He[:-1] + '-%d/' % nB for nB in StreuBases]
for finalStatePath in finalStatePaths:
    if os.path.isdir(finalStatePath) == True:
        print('<ECCE> removing the existing final-state folder\n%s' %
              finalStatePath)
        os.system('rm -rf ' + finalStatePath)
    subprocess.check_call(['mkdir -p', finalStatePath])
    subprocess.check_call(['mkdir -p', finalStatePath + 'basis_struct/'])

# > optimize the various basis types, e.g., in case of the npp system:
# > helion ground state, final J=1/2- and J=3/2- states
for bastype in bastypes:

    Jay = float(bastype.split('^')[0][-3:])
    Jaystring = '%s' % str(Jay)[:3]

    # number of final-state bases which are grown with the above-set criteria

    costr = ''
    zop = 31 if tnni == 11 else 14
    for nn in range(1, zop):
        if bastype == boundstatekanal:
            cf = 1.0 if (1 <= nn <= 28) else 0.0
        else:
            cf = 1.0 if (1 <= nn <= 28) else 0.0
        costr += '%12.7f' % cf if (nn % 7 != 0) else '%12.7f\n' % cf

    # numerical stability
    minDiffwidthsINT = 10**-2
    minDiffwidthsREL = 10**-3
    maxParLen = 18

    # rating criterion in order to stratify the initial seed basis
    # default and reasonable for an unstable, random seed =0 (C-nbr)
    pwpurge = 0

    # rating criterion in order to stratify the *stabilized* initial seed basis
    # ensuing the one-time optimization of each of its configurations
    # here, we want to sort out idling vectors, i.e, such which do neither harm by causing
    # small C-nbrs nor increase the fitness significantly
    pwSig = 1

    # rating criterion for the offspring, i.e., based on what measure is a single basis-vector block
    # added to the parent population?
    # at present, any function of the norm/Hamiltonian matrices/spectra are conceivable
    # ECCE: criteria derived from the overlap with some input initial state, e.g., E1|A-body bound state>
    #       is not yet implemented!
    pwopt = 1
    purgeStr = ['Condition number', 'pulchritude',
                'Ground-state energy'][pwpurge]

    # evolution criteria
    minCond = 10**-6
    denseEVinterval = [-10., 150.0]
    removalGainFactor = 1.5
    maxOnPurge = 43
    maxOnTrail = 10**2
    muta_initial = 0.5

    # get the initial, random basis seed to yield thresholds close to the reuslts in a complete basis
    chThreshold = -6.5 if bastype == boundstatekanal else -1.5

    CgfCycles = 3
    # nRaces := |i|
    nRaces = 4 if bastype == boundstatekanal else 8

    cradleCapacity = 42

    # > nState > produce/optimize/grow multiple bases with pseudo-random initial seeds
    for nB in range(anzStreuBases):

        wrkDir = helionpath if bastype == boundstatekanal else finalStatePaths[
            nB]
        basisPath = wrkDir + 'basis_struct/'

        os.chdir(wrkDir)

        os.system('cp %s .' % potnn)
        os.system('cp %s .' % potnnn)

        gsEnergy = 42.0

        while gsEnergy >= chThreshold:

            t0 = time.perf_counter()

            seedMat = span_initial_basis(basisType=bastype,
                                         ini_grid_bounds=[
                                             0.006, 7.25, 0.004, 6.5, 0.005,
                                             5.25, 0.001, 4.5
                                         ],
                                         ini_dims=[8, 8, 8, 8],
                                         coefstr=costr,
                                         anzOp=zop)

            t1 = time.perf_counter()
            print(f"Seed basis generation in {np.abs(t0 - t1):0.4f} seconds.")

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
                bvv for bvv in ewH
                if denseEVinterval[0] < bvv < denseEVinterval[1]
            ])
            gsEnergy = ewH[-1]
            basCond = np.min(np.abs(ewN)) / np.max(np.abs(ewN))
            attractiveness = loveliness(gsEnergy, basCond, anzSigEV, minCond)

            print(
                '\n> basType %s > basSet %d/%d: seed basis: E0 = %f   cond=|Emin|/|Emax| = %e'
                % (bastype, nB + 1, anzStreuBases, gsEnergy, basCond))

            if gsEnergy >= chThreshold:
                print(
                    'ECCE! seed does not expand states with E<%f => new sowing attempt.'
                    % chThreshold)

                subprocess.call('rm -rf DMOUT.* && rm -rf DRDMOUT.*',
                                shell=True)
                subprocess.call('rm -rf TQUAOUT.* && rm -rf TDQUAOUT.*',
                                shell=True)
                continue

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

            # set of unique angular, spin, and isospin configurations
            # ecce: each of these cfg's might appear multiple times if the
            # number of radial widths associated with it exceeds <bvma>
            unisA = []
            for ncfg in range(len(initialCiv[0])):
                if initialCiv[0][ncfg] in unisA:
                    continue
                else:
                    unisA.append(initialCiv[0][ncfg])

            nbv = 0
            for cfg in range(len(initialCiv[0])):
                nbvc = 0
                for bv in initialCiv[1][cfg]:
                    nbv += 1
                    nbvc += 1
                    initialCiv[3] += [[
                        nbv,
                        np.array(
                            range(1, 1 +
                                  len(initialCiv[2][cfg][nbvc - 1]))).tolist()
                    ]]

            # > nState > nBasis > stabilize the seed basis
            goPurge = True

            #print('\n\nSeed Basis (naive):\n\n', initialCiv)

            initialCiv = condense_basis(initialCiv, MaxBVsPERcfg=10)

            #print('\n\nSeed Basis (condensed):\n\n', initialCiv, '\n\n')

            D0 = initialCiv[3]

            # purge just entire bv sets with identical internal width
            print(
                '\n> basType %s > basSet %d/%d: Stratifying the initial seed -- criterion: %s'
                % (bastype, nB + 1, anzStreuBases, purgeStr))
            t0 = time.perf_counter()
            while goPurge:
                newpopList = []
                goPurge = False
                ParaSets = []

                ParaSets.append([
                    D0, Jay, costr, zop, 10, [0, 0], BINBDGpath, minCond,
                    denseEVinterval
                ])

                for bvTrail in D0:

                    bvID = [
                        int(bvTrail[0]),
                        int(''.join(map(str, bvTrail[1])))
                    ]
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

                if dbg:
                    print(
                        '   rating basis vectors in %d-dim basis on their effect on the stability,'
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

                if dbg:
                    for cand in cand_ladder[-3:]:
                        print(cand[:4])

                print(
                    '\n> basType %s > basSet %d/%d: purged seed: E0 = %f   C-nbr=|Emin|/|Emax| = %e'
                    % (bastype, nB + 1, anzStreuBases, cand_ladder[-1][2],
                       cand_ladder[-1][0]))

                if dbg:
                    print('    best =  %2.3e , %2.3e , %2.3e' %
                          (cand_ladder[-1][0], cand_ladder[-1][1],
                           cand_ladder[-1][2]))
                    print('    worst=  %2.3e , %2.3e , %2.3e' %
                          (cand_ladder[0][0], cand_ladder[0][1],
                           cand_ladder[0][2]))

                # only if the removal of a basis-vector block (see [Kir Dipl, p.38])
                # 1) increases the cond. number significantly, or,
                # 2) in case of a prohibitively unstable reference, any removal which
                #    stabilizes the set above a preset threshold (minCond) is acceptable
                if ((removalGainFactor * np.abs(reff[pwpurge]) < np.abs(
                        cand_ladder[-1][pwpurge])) |
                    ((np.abs(cand_ladder[-1][pwpurge]) < minCond) &
                     (reff[0] < minCond))):
                    goPurge = True
                    D0 = rectify_basis(cand_ladder[-1][4])
                    print(
                        '                   removal of 1/%d basis-vector blocks is advantageous.'
                        % len(D0),
                        end='')

            t1 = time.perf_counter()
            print(
                f"Seed basis generation stabilized in {np.abs(t0 - t1):0.4f} seconds."
            )

            subprocess.call('rm -rf DMOUT.* && rm -rf DRDMOUT.*', shell=True)
            subprocess.call('rm -rf TQUAOUT.* && rm -rf TDQUAOUT.*',
                            shell=True)

            initialCiv[3] = rectify_basis(cand_ladder[-1][4])
            # > nState > nBasis > end of stabilization

            initialCiv = essentialize_basis(initialCiv, MaxBVsPERcfg=bvma)

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
                          potNN='./%s' % nnStr,
                          potNNN='./%s' % nnnStr,
                          parall=-1,
                          anzcores=max(2, min(len(initialCiv[0]), MaxProc)),
                          tnnii=tnni,
                          jay=Jay,
                          dia=False)

            ewN, ewH = NormHamDiag(ma)

            parCond = np.min(np.abs(ewN)) / np.max(np.abs(ewN))

            gsEnergy = ewH[-1]

            print(
                '\n> basType %s > basSet %d/%d: stabilized initial basis: C-nbr = %4.4e E0 = %4.4e\n\n>>> COMMENCING OPTIMIZATION <<<\n'
                % (bastype, nB + 1, anzStreuBases, parCond, gsEnergy))

            # count unique cfg's after purge and exit() if the purge removed one of them
            # entirely
            unis = []
            for ncfg in range(len(initialCiv[0])):
                if initialCiv[0][ncfg] in unis:
                    continue
                else:
                    unis.append(initialCiv[0][ncfg])
            if len(unis) != len(unisA):
                print(
                    'Elemental cfg of the seed was removed entirely during purge.\n new round of sowing.'
                )
                gsEnergy = 42.0

        # > nState > nBasis > optimize each orb-ang, spin-iso cfg in a number of cycles
        for nCgfCycle in range(CgfCycles):

            # > nState > nBasis > nCfgCycle > optimize a single angular-momentum configuration, e.g. l1=1,l2=1,L=2,s12=....

            for nUcfg in range(len(unis)):

                print(
                    '> basType %s > basSet %d/%d > cfgCycle %d/%d > nUcfg %d/%d > : Optimizing cfg: '
                    % (bastype, nB + 1, anzStreuBases, nCgfCycle + 1,
                       CgfCycles, nUcfg + 1, len(unis)), unis[nUcfg])
                # > nState > nBasis > nCfgCycle > nUcfg > allow each cfg to evolve over nRaces
                for nGen in range(nRaces):

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

                    print(
                        '\n> basType %s > basSet %d/%d > cfgCycle %d/%d > nUcfg %d/%d > nGen %d/%d: breeding %d basis-vector blocks from %d parent blocks.'
                        % (bastype, nB + 1, anzStreuBases, nCgfCycle + 1,
                           CgfCycles, nUcfg + 1, len(unis), nGen + 1, nRaces,
                           cradleCapacity, len(parentIWs)))

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
                                  potNN='./%s' % nnStr,
                                  potNNN='./%s' % nnnStr,
                                  parall=-1,
                                  anzcores=max(2, min(len(Ais[0]), MaxProc)),
                                  tnnii=tnni,
                                  jay=Jay,
                                  dia=False)

                    ewN, ewH = NormHamDiag(ma)

                    parLove, parCond = basQ(ewN, ewH, minCond, denseEVinterval)

                    if ewH == []:
                        print('old generation already unstable:\n', Ais)
                        exit()

                    print(
                        'reference for the new gen: Dim(parents+offspring) = %d; parents: B(GS) = %8.4f  C-nbr = %4.3e  fitness = %4.3e'
                        % (basisDim(Ais[3]), ewH[-1], parCond, parLove))

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

                    subprocess.call('rm -rf DMOUT.* && rm -rf DRDMOUT.*',
                                    shell=True)
                    subprocess.call('rm -rf TQUAOUT.* && rm -rf TDQUAOUT.*',
                                    shell=True)

                    cand_ladder = [x.recv() for x in cand_list]

                    # ranking following condition-number (0) or quality (1)  or E(GS) (2)
                    cand_ladder.sort(key=lambda tup: np.abs(tup[pwopt]))

                    if dbg:
                        for cand in cand_ladder[-3:]:
                            print(cand[:4])
                    # optimum when deleting one

                    cand = cand_ladder[-1]

                    if ((cand[3] != [-1]) & (cand[0] > minCond)):
                        parvenue = cand

                        if dbg: print('\nparents:\n', parLove, parCond)
                        print('\nparvenue:\n', parvenue[:4])
                        Ais[3] = parvenue[-1]

                        initialCivL = essentialize_basis(Ais,
                                                         MaxBVsPERcfg=bvma)

                        #print('\n\niniciV+optChild (naive):\n\n', initialCivL)

                        initialCiv = condense_basis(initialCivL,
                                                    MaxBVsPERcfg=bvma)

                        #print('\n\niniciV+optChild (strat):\n\n', initialCiv)
                    else:
                        Ais[3] = parBVs
                        print(
                            'All children are spoiled. Starting anew with the initial civilization.\n'
                        )

            initialCivL = essentialize_basis(initialCiv, MaxBVsPERcfg=bvma)

            initialCiv = condense_basis(initialCivL, MaxBVsPERcfg=bvma)

            # after having evolved each configuration, stabilize the basis
            # and remove its least siginificant vectors before repeating
            # the optimization of the individual configurations
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
                          potNN='./%s' % nnStr,
                          potNNN='./%s' % nnnStr,
                          parall=-1,
                          anzcores=max(2, min(len(initialCiv[0]), MaxProc)),
                          tnnii=tnni,
                          jay=Jay,
                          dia=False)

            D0 = initialCiv[3]
            # purge just entire bv sets with identical internal width

            print(
                '\n> basType %s > basSet %d/%d > cfgCycle %d/%d: Final stabilization, i.e, removal of insignificant(%s) basis-vector blocks.'
                % (bastype, nB + 1, anzStreuBases, nCgfCycle + 1, CgfCycles,
                   purgeStr))

            goPurge = True
            while goPurge:
                newpopList = []
                goPurge = False
                ParaSets = []
                ParaSets.append([
                    D0, Jay, costr, zop, 10, [0, 0], BINBDGpath, minCond,
                    denseEVinterval
                ])
                for bvTrail in D0:
                    bvID = [
                        int(bvTrail[0]),
                        int(''.join(map(str, bvTrail[1])))
                    ]
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
                cand_ladder.sort(key=lambda tup: np.abs(tup[pwSig]))
                reff = [[bvm[0], bvm[1], bvm[2]] for bvm in cand_ladder
                        if bvm[3] == [0, 0]][0]
                for cand in cand_ladder[-3:]:
                    print(cand[:4])
                print(
                    '\n> basType %s > basSet %d/%d: purged seed: E0 = %f   cond=|Emin|/|Emax| = %e'
                    % (bastype, nB + 1, anzStreuBases, cand_ladder[-1][2],
                       cand_ladder[-1][0]))
                print('    best =  %2.3e , %2.3e , %2.3e' %
                      (cand_ladder[-1][0], cand_ladder[-1][1],
                       cand_ladder[-1][2]))
                print(
                    '    worst=  %2.3e , %2.3e , %2.3e' %
                    (cand_ladder[0][0], cand_ladder[0][1], cand_ladder[0][2]))

                # more fastidious significance threshold because we subjected the vectors
                # added to the initial seed already to a vetting process
                removalGainFactor2 = 2 * removalGainFactor

                if ((removalGainFactor2 * np.abs(reff[pwSig]) < np.abs(
                        cand_ladder[-1][pwSig])) |
                    ((np.abs(cand_ladder[-1][pwSig]) < minCond) &
                     (reff[0] < minCond))):
                    goPurge = True
                    D0 = rectify_basis(cand_ladder[-1][4])
                    print('removing 1/%d basis-vector blocks.' % len(D0),
                          end='\n')

            subprocess.call('rm -rf DMOUT.* && rm -rf DRDMOUT.*', shell=True)
            subprocess.call('rm -rf TQUAOUT.* && rm -rf TDQUAOUT.*',
                            shell=True)
            initialCiv[3] = rectify_basis(cand_ladder[-1][4])
            # > nState > nBasis > end of stabilization
            initialCivL = essentialize_basis(initialCiv, MaxBVsPERcfg=bvma)

            initialCiv = condense_basis(initialCivL, MaxBVsPERcfg=bvma)

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
                      potNN='./%s' % nnStr,
                      potNNN='./%s' % nnnStr,
                      parall=-1,
                      anzcores=max(2, min(len(initialCiv[0]), MaxProc)),
                      tnnii=tnni,
                      jay=Jay,
                      dia=False)

        ewN, ewH = NormHamDiag(ma)

        optLove, optCond = basQ(ewN, ewH, minCond, denseEVinterval)

        print(
            '\n> basType %s > basSet %d/%d: optimized basis: C-nbr = %4.4e E0 = %4.4e fitness = %4.4e\n\n'
            % (bastype, nB + 1, anzStreuBases, optCond, ewH[-1], optLove))

        # Output on tape; further processing via A3...py
        suf = 'ref' if bastype == boundstatekanal else 'fin'

        lfrags = np.array(initialCiv[0])[:, 1].tolist()
        sfrags = np.array(initialCiv[0])[:, 0].tolist()
        n3_inlu(8, fn=basisPath + 'INLU_%s' % suf, fr=lfrags, indep=-1)
        n3_inlu(8, fn=basisPath + 'INLUCN_%s' % suf, fr=lfrags, indep=-1)
        n3_inob(sfrags, 8, fn=basisPath + 'INOB_%s' % suf, indep=-1)
        n3_inob(sfrags, 15, fn=basisPath + 'DRINOB_%s' % suf, indep=-1)

        subprocess.call('cp INQUA_M ' + basisPath + 'INQUA_V18_%s' % suf,
                        shell=True)
        subprocess.call('cp INEN ' + basisPath + 'INEN_%s' % suf, shell=True)
        subprocess.call('cp INSAM ' + basisPath, shell=True)

        subprocess.call('rm -rf ./inen_*', shell=True)
        subprocess.call('rm -rf ./endout_*', shell=True)
        subprocess.call('rm -rf ./MATOUTB_*', shell=True)
        subprocess.call('rm -rf DMOUT.* && rm -rf DRDMOUT.*', shell=True)
        subprocess.call('rm -rf TQUAOUT.* && rm -rf TDQUAOUT.*', shell=True)

        fullBasfile, actBasfile = write_basis_on_tape(initialCiv,
                                                      Jay,
                                                      bastype,
                                                      baspath=basisPath)

        if bastype != boundstatekanal:
            AbasOutStr = respath + 'Ssigbasv3heLIT_%s_BasNR-%d.dat' % (
                bastype, StreuBases[nB])
            FbasOutStr = respath + 'SLITbas_full_%s_BasNR-%d.dat' % (
                bastype, StreuBases[nB])
            subprocess.call('cp %s %s' % (fullBasfile, FbasOutStr), shell=True)
            subprocess.call('cp %s %s' % (actBasfile, AbasOutStr), shell=True)
        else:
            AbasOutStr = respath + 'Ssigbasv3heLIT_%s.dat' % bastype
            FbasOutStr = respath + 'SLITbas_full_%s.dat' % bastype
            subprocess.call('cp %s %s' % (fullBasfile, FbasOutStr), shell=True)
            subprocess.call('cp %s %s' % (actBasfile, AbasOutStr), shell=True)

        matoutstr = '%smat_%s' % (
            respath, bastype + '_BasNR-%d' % StreuBases[nB]
        ) if bastype != boundstatekanal else '%smat_%s' % (respath, bastype)
        subprocess.call('cp MATOUTB %s' % matoutstr, shell=True)
        print(
            'channel %s: Basis structure, Norm, and Hamiltonian written into %s'
            % (bastype, respath + 'mat_' + bastype))

        srcDir = litpath3He[:-1]
        os.system('rsync -r -u %s %s' % (srcDir, bkpdir))
        if bastype != boundstatekanal:
            os.system('rsync -r -u %s %s' % (wrkDir[:-1], bkpdir))
        # for the bound-state/initial-state channel, consider only one basis set
        if bastype == boundstatekanal:
            break