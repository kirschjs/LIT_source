from bridgeA3 import *
from genetic_width_growth import *

import glob, copy
import shlex
import multiprocessing
from multiprocessing.pool import ThreadPool
from itertools import permutations, product

from PSI_parallel_M import span_initial_basis

if os.path.isdir(litpath3He) == False:
    subprocess.check_call(['mkdir', '-p', litpath3He])
    subprocess.check_call(['mkdir', '-p', respath])

with open(respath + 'dtype.dat', 'w') as outf:
    outf.write(dt)

os.chdir(litpath3He)

dbg = False

# maximum number of attempts to furnish a random initial basis which
# satisfies the defined minimal quality and stability criteria
max_iter = 12

# call with arg1<0 : boundsate
#           a b    : streubases from a to b
#           a a    : a single basis is optimized
arglist = sys.argv

if arglist[1:] != []:
    # for par_run.py operation
    StreuBases = np.arange(int(arglist[1]), int(arglist[2]) + 1)
    anzStreuBases = len(StreuBases)
    bastypes = [boundstatekanal] if int(arglist[1]) < 0 else streukas
else:
    # for manual operation
    anzStreuBases = 14
    StreuBases = np.arange(1, anzStreuBases + 1)
    bastypes = [boundstatekanal] + streukas  #

if boundstatekanal in bastypes:
    if os.path.isdir(helionpath) != False:
        print('<ECCE> removing the existing helion folder\n%s.' % helionpath)
        os.system('rm -rf ' + helionpath)
    subprocess.check_call(['mkdir', '-p', helionpath])
    subprocess.check_call(['mkdir', '-p', helionpath + 'basis_struct/'])

for streuka in streukas:
    if streuka in bastypes:
        finalStatePaths = [litpath3He[:-1] + '-%d/' % nB for nB in StreuBases]
        for finalStatePath in finalStatePaths:
            if os.path.isdir(finalStatePath) == True:
                print('<ECCE> removing the existing final-state folder\n%s' %
                      finalStatePath)
                os.system('rm -rf ' + finalStatePath)
            subprocess.check_call(['mkdir', '-p', finalStatePath])
            subprocess.check_call(
                ['mkdir', '-p', finalStatePath + 'basis_struct/'])
        break

# > optimize the various basis types, e.g., in case of the npp system:
# > helion ground state, final J=1/2- and J=3/2- states
for bastype in bastypes:

    # ini_dims = [BS(int),BS(rel),SCATT(int),SCATT(rel)]
    init_dims = [8, 20, 8, 24]

    # lower and upper bounds for the grids from which the initial seed state is taken
    # 1-4: initial state, 1-2(jacobi1), 3-4(jacobi2)
    # 5-8: final   states,5-6(jacobi1), 7-8(jacobi2)
    ini_grid_bnds = [0.2, 19.25, 0.001, 18.5, 0.01, 19.25, 0.001, 18.5]

    Jay = float(bastype.split('^')[0][-3:])
    Jaystring = '%s' % str(Jay)[:3]

    # number of final-state bases which are grown with the above-set criteria

    costr = ''
    zop = 31 if tnni == 11 else 14
    for nn in range(1, zop):
        #if bastype == boundstatekanal:
        #    cf = 1.0 if (1 <= nn <= 28) else 0.0
        #else:
        #    cf = 1.0 if (1 <= nn <= 28) else 0.0
        cf = 1.0 if (nn != 8) else 1.
        # for contact interaction
        #cf = 1.0 if ((nn == 14) | (nn == 2) | (nn == 1)) else 0.
        costr += '%12.7f' % cf if (nn % 7 != 0) else '%12.7f\n' % cf

    # numerical stability
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
    minCond = 10**-14
    denseEVinterval = [-10., 180.0]
    removalGainFactor = 1.05
    maxOnPurge = 113
    maxOnTrail = 43
    muta_initial = 0.075

    BDGdeu = 2.224
    BDG3h = 8.482
    BDG3he = 7.72
    # get the initial, random basis seed to yield thresholds close to the reuslts in a complete basis
    chThreshold = -6.51 if bastype == boundstatekanal else -1.6

    CgfCycles = 2
    # nRaces := |i|
    nRaces = 1 if bastype == boundstatekanal else 3

    cradleCapacity = 6

    # > nState > produce/optimize/grow multiple bases with pseudo-random initial seeds
    for nB in range(anzStreuBases):

        wrkDir = helionpath if bastype == boundstatekanal else finalStatePaths[
            nB]
        basisPath = wrkDir + 'basis_struct/'

        os.chdir(wrkDir)

        os.system('cp %s .' % potnn)
        os.system('cp %s .' % potnnn)

        gsEnergy = 42.0

        seed_attempts = 0
        while ((gsEnergy >= chThreshold) | (gsEnergy < -1.2 * BDG3he) &
               (seed_attempts < max_iter)):

            seed_attempts += 1
            t0 = time.perf_counter()

            try:
                wrkVol = du(pathbase)
                while int(wrkVol) > homeQuota:
                    print('wrkDir holds %d bytes. Waiting for 60s to shrink.' %
                          int(wrkVol))
                    time.sleep(60)
                    wrkVol = du(pathbase)
            except:
                print(
                    '(ecce) disk-usage-assessment failure. I will continue, aware of the increased crash risk!'
                )

            # ini_grid_bnds = [bs_int_low,bs_int_up,bs_rel_low,bs_rel_up,SC_int_low,SC_int_up,SC_rel_low,SC_rel_up]

            seedMat = span_initial_basis(basisType=bastype,
                                         ini_grid_bounds=ini_grid_bnds,
                                         ini_dims=init_dims,
                                         coefstr=costr,
                                         anzOp=zop)

            t1 = time.perf_counter()
            print(
                f"%d-Seed basis generation in {np.abs(t0 - t1):0.4f} seconds."
                % seed_attempts)

            smartEV, basCond = smart_ev(seedMat, threshold=10**-9)
            # > nState > nBasis > stabilize the seed basis
            goPurge = True if (basCond < minCond) else False

            anzSigEV = len([
                bvv for bvv in smartEV
                if denseEVinterval[0] < bvv < denseEVinterval[1]
            ])
            gsEnergy = smartEV[-1]
            attractiveness = loveliness(gsEnergy, basCond, anzSigEV, minCond)

            print(
                '\n> basType %s > basSet %d/%d: seed basis: E0 = %f   cond=|Emin|/|Emax| = %e'
                % (bastype, nB + 1, anzStreuBases, gsEnergy, basCond))

            if ((gsEnergy >= chThreshold) | (gsEnergy < -1.2 * BDG3he)):
                print(
                    'ECCE! seed does not expand states with E<%f => new sowing attempt.'
                    % chThreshold)

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

            #print('\n\nSeed Basis (naive):\n\n', initialCiv)

            initialCiv = condense_basis(initialCiv, MaxBVsPERcfg=bvma)

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
                    D0, Jay, costr, zop, tnni, [0, 0], BINBDGpath, minCond,
                    denseEVinterval
                ])

                for bvTrail in D0:

                    if len(bvTrail[1]) > 1:
                        bvID = [
                            int(bvTrail[0]),
                            int(''.join(map(str, bvTrail[1])))
                        ]
                        cpy = copy.deepcopy(D0)
                        cpy.remove(bvTrail)

                        ParaSets.append([
                            cpy, Jay, costr, zop, tnni, bvID, BINBDGpath,
                            minCond, denseEVinterval
                        ])

                #for ca in ParaSets:
                #    print(ca[0], '\n\n')

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
                reff = [[bvm[0], bvm[1], bvm[2]] for bvm in cand_ladder
                        if bvm[3] == [0, 0]][0]

                # ranking following condition-number (0) or quality (1)
                condTh = minCond
                gsDiff = 0.005
                empt = True
                while empt:
                    testlist = [
                        cand for cand in cand_ladder
                        if ((cand[0] > condTh)
                            & (np.abs(cand[2] - gsEnergy) < gsDiff))
                    ]
                    if len(testlist) > 3:
                        stab_ladder = testlist
                        empt = False
                    condTh = condTh * 0.5
                    gsDiff += 0.001

                #for cand in stab_ladder:  #[-3:]:
                #    print(cand[:4])

                stab_ladder.sort(key=lambda tup: np.abs(tup[pwpurge]))

                print(
                    '\n> basType %s > basSet %d/%d: purged seed: E0 = %f   C-nbr=|Emin|/|Emax| = %e'
                    % (bastype, nB + 1, anzStreuBases, stab_ladder[-1][2],
                       stab_ladder[-1][0]))

                if dbg:
                    print('    best =  %2.3e , %2.3e , %2.3e' %
                          (stab_ladder[-1][0], stab_ladder[-1][1],
                           stab_ladder[-1][2]))
                    print('    worst=  %2.3e , %2.3e , %2.3e' %
                          (stab_ladder[0][0], stab_ladder[0][1],
                           stab_ladder[0][2]))

                print(
                    'removed of 1/%d basis-vector blocks to gain stability.' %
                    len(D0),
                    end='')

                if ((np.abs(stab_ladder[-1][0]) < minCond)):
                    goPurge = True
                    newBas = stab_ladder[-1][4] if stab_ladder[-1][3] != [
                        0, 0
                    ] else stab_ladder[-2][4]
                    D0 = rectify_basis(newBas)
                else:
                    goPurge = False
                    D0 = rectify_basis(stab_ladder[-1][4])

            t1 = time.perf_counter()
            print(
                f"\n\nSeed basis generation stabilized in {np.abs(t0 - t1):0.4f} seconds."
            )

            try:
                initialCiv[3] = rectify_basis(stab_ladder[-1][4])
                # > nState > nBasis > end of stabilization
                initialCiv = essentialize_basis(initialCiv, MaxBVsPERcfg=bvma)
            except:
                print(
                    'empty candidate ladder. I will use the un-purged civilization.'
                )

            try:
                wrkVol = du(pathbase)
                while int(wrkVol) > homeQuota:
                    print('wrkDir holds %d bytes. Waiting for 60s to shrink.' %
                          int(wrkVol))
                    time.sleep(60)
                    wrkVol = du(pathbase)
            except:
                print(
                    '(ecce) disk-usage-assessment failure. I will continue, aware of the increased crash risk!'
                )

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
                          jay=Jay)

            smartEV, parCond = smart_ev(ma, threshold=10**-9)

            gsEnergy = smartEV[-1]

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

            print(unis, unisA)
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

                            for nrw in range(
                                    np.min([len(pa) for pa in parentRWs])):
                                #for nrw in range(len(parentRWs[iws[0]])):
                                motherR = parentRWs[iws[0]][nrw]
                                if iws[1] >= 0:
                                    #print(iws[1], nrw, '\n', ret, unis, '\n',
                                    #      initialCiv[0])
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
                        parBVs, Jay, costr, zop, tnni, [-1], BINBDGpath,
                        minCond, denseEVinterval
                    ])

                    for bvTrail in chBVs:
                        childidentifier = [bvTrail[0]]
                        cpy = copy.deepcopy(parBVs)
                        cpy.append(bvTrail)
                        cpy = rectify_basis(cpy)
                        cpy.sort()
                        childishParaSets.append([
                            cpy, Jay, costr, zop, tnni, childidentifier,
                            BINBDGpath, minCond, denseEVinterval
                        ])

                    # do this only if all children are in basis[3], otherwise their widths are removed!
                    # do it (!) to break frgs which became too long after the children were added

                    #print('\nAis:\n', Ais)
                    Ais = essentialize_basis(Ais, MaxBVsPERcfg=bvma)
                    #print('\nAis (strat):\n', Ais)

                    subprocess.call('rm -rf TQUAOUT.*', shell=True)
                    subprocess.call('rm -rf TDQUAOUT.*', shell=True)
                    try:
                        wrkVol = du(pathbase)
                        while int(wrkVol) > homeQuota:
                            print(
                                'wrkDir holds %d bytes. Waiting for 60s to shrink.'
                                % int(wrkVol))
                            time.sleep(60)
                            wrkVol = du(pathbase)
                    except:
                        print(
                            '(ecce) disk-usage-assessment failure. I will continue, aware of the increased crash risk!'
                        )

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
                                  jay=Jay)

                    smartEV, parCond = smart_ev(ma, threshold=10**-9)
                    anzSigEV = len([
                        bvv for bvv in smartEV
                        if denseEVinterval[0] < bvv < denseEVinterval[1]
                    ])
                    gsEnergy = smartEV[-1]
                    parLove = loveliness(gsEnergy, parCond, anzSigEV, minCond)

                    print(
                        'reference for the new gen: Dim(parents+offspring) = %d; parents: B(GS) = %8.4f  C-nbr = %4.3e  fitness = %4.3e'
                        % (basisDim(Ais[3]), gsEnergy, parCond, parLove))

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

                    if dbg:
                        for cand in cand_ladder[-3:]:
                            print(cand[:4])
                    # optimum when deleting one

                    cand = cand_ladder[-1]

                    if ((cand[3] != [-1]) & (cand[0] > minCond * 10**-2)):
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

            try:
                wrkVol = du(pathbase)
                while int(wrkVol) > homeQuota:
                    print('wrkDir holds %d bytes. Waiting for 60s to shrink.' %
                          int(wrkVol))
                    time.sleep(60)
                    wrkVol = du(pathbase)
            except:
                print(
                    '(ecce) disk-usage-assessment failure. I will continue, aware of the increased crash risk!'
                )

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
                          jay=Jay)

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
                    D0, Jay, costr, zop, tnni, [0, 0], BINBDGpath, minCond,
                    denseEVinterval
                ])
                for bvTrail in D0:

                    if len(bvTrail[1]) > 1:
                        bvID = [
                            int(bvTrail[0]),
                            int(''.join(map(str, bvTrail[1])))
                        ]
                        cpy = copy.deepcopy(D0)
                        cpy.remove(bvTrail)
                        ParaSets.append([
                            cpy, Jay, costr, zop, tnni, bvID, BINBDGpath,
                            minCond, denseEVinterval
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
                        cand_ladder[-1][pwSig])) &
                    ((np.abs(cand_ladder[-1][pwSig]) < minCond)
                     #&(reff[0] < minCond)
                     )):
                    goPurge = True

                    if np.min([len(bv[1]) for bv in cand_ladder[-1][4]]) < 1:
                        print('%$**&!@#:  ',
                              [len(bv[1]) for bv in cand_ladder[-1][4]])
                        #exit()
                    D0 = rectify_basis(cand_ladder[-1][4])
                    print('removing 1/%d basis-vector blocks.' % len(D0),
                          end='\n')

            initialCiv[3] = rectify_basis(cand_ladder[-1][4])
            # > nState > nBasis > end of stabilization
            initialCivL = essentialize_basis(initialCiv, MaxBVsPERcfg=bvma)

            initialCiv = condense_basis(initialCivL, MaxBVsPERcfg=bvma)

        try:
            wrkVol = du(pathbase)
            while int(wrkVol) > homeQuota:
                print('wrkDir holds %d bytes. Waiting for 60s to shrink.' %
                      int(wrkVol))
                time.sleep(60)
                wrkVol = du(pathbase)
        except:
            print(
                '(ecce) disk-usage-assessment failure. I will continue, aware of the increased crash risk!'
            )

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
                      jay=Jay)

        smartEV, optCond = smart_ev(ma, threshold=10**-9)
        anzSigEV = len([
            bvv for bvv in smartEV
            if denseEVinterval[0] < bvv < denseEVinterval[1]
        ])
        gsEnergy = smartEV[-1]
        optLove = loveliness(gsEnergy, parCond, anzSigEV, minCond)

        print(
            '\n> basType %s > basSet %d/%d: optimized basis: C-nbr = %4.4e E0 = %4.4e fitness = %4.4e\n\n'
            % (bastype, nB + 1, anzStreuBases, optCond, gsEnergy, optLove))

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

        subprocess.call('rm -rf TQUAOUT.*', shell=True)
        subprocess.call('rm -rf TDQUAOUT.*', shell=True)

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