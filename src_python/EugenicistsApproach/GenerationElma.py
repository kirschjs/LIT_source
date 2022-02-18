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

bastypes = streukas

HelBasDim = sum([
    len(ln.split()[1:])
    for ln in open(helionpath +
                   'basis_struct/LITbas_full_%s.dat' % boundstatekanal)
])
HeFrags = []
HeFrags = []
fragfile = [
    ln for ln in open(helionpath +
                      'basis_struct/frags_LIT_%s.dat' % boundstatekanal)
]

HeFrags = [fr.split(' ')[1].strip() for fr in fragfile]
HeFrags = [fr.split(' ')[0] for fr in fragfile]

for bastype in bastypes:
    # number of final-state bases which are grown with the above-set criteria
    anzStreuBases = 1

    costr = ''
    zop = 31 if tnni == 11 else 14
    for nn in range(1, zop):
        cf = 1.0 if (1 <= nn <= 28) else 0.0
        costr += '%12.7f' % cf if (nn % 7 != 0) else '%12.7f\n' % cf

    finalStatePaths = [
        litpath3He[:-1] + '-%d/' % nB for nB in range(anzStreuBases)
    ]
    for finalStatePath in finalStatePaths:
        if os.path.isdir(finalStatePath) != True:
            os.system('rm -rf ' + finalStatePath)
            print('No hamilton-optimized basis found in%s .' % finalStatePath)
            exit()

    # numerical stability
    minDiffwidthsINT = 10**-2
    minDiffwidthsREL = 10**-3
    maxParLen = 18

    # evolution criteria
    minCond = 10**-9
    denseEVinterval = [-10., 150.0]
    removalGainFactor = 1.5
    maxOnPurge = 43
    muta_initial = 0.75

    # nRaces := |i|
    nRaces = 5

    nbrOff = 16
    MaxOff = 22

    for nB in range(anzStreuBases):

        wrkDir = finalStatePaths[nB]
        basisPath = wrkDir + 'basis_struct/'

        os.chdir(wrkDir)

        print('\n>>> Basistype: %s\n >> Basis Set number: %d/%d ' %
              (bastype, nB + 1, anzStreuBases))

        Jay = float(bastype.split('^')[0][-3:])
        Jaystring = '%s' % str(Jay)[:3]

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

        # 3) breed the first offspring generation
        # each pair of basis vectors procreates => |parents + offspring| = 2*|parents|
        # P_0 -> {P_0^i=parents+offspring_i,i=1,...,nRaces}

        Civilizations = []
        Civilizations.append([cfgs, intwLIT, rws, []])
        nbv = 0
        for cfg in range(len(Civilizations[0][0])):
            nbvc = 0
            for bv in Civilizations[0][1][cfg]:
                nbv += 1
                nbvc += 1
                Civilizations[0][3] += [[
                    nbv,
                    np.array(
                        range(
                            1, 1 +
                            len(Civilizations[0][2][cfg][nbvc - 1]))).tolist()
                ]]

        nCivi = 0
        if nbv > 1:

            for nCivi in range(nRaces):

                cfgbnds = np.add.accumulate(
                    [len(iws) for iws in Civilizations[nCivi][1]])
                cfgbnds = np.insert(cfgbnds, 0, 0)

                lfragTNG = np.array(Civilizations[nCivi][0])[:, 1].tolist()
                sfragTNG = np.array(Civilizations[nCivi][0])[:, 0].tolist()
                insam(len(lfragTNG))

                # i) each initial parent-offspring basis is likely unstable if comprised of all 'families'(Dfull)

                ma = blunt_ev(Civilizations[-1][0],
                              Civilizations[-1][1],
                              Civilizations[-1][2],
                              Civilizations[-1][3],
                              wrkdir='',
                              nzopt=zop,
                              costring=costr,
                              bin_path=BINBDGpath,
                              mpipath=MPIRUN,
                              einzel_file_path=wrkDir,
                              potNN=potnn,
                              potNNN=potnnn,
                              parall=-1,
                              anzcores=max(2, min(len(lfragTNG), MaxProc)),
                              tnni=10,
                              jay=Jay,
                              dia=False)

                ewN, ewH = NormHamDiag(ma)

                parCond = np.min(np.abs(ewN)) / np.max(np.abs(ewN))
                print(
                    ' >> hamiltonian-optimized basis: E0 = %f   cond=|Emin|/|Emax| = %e'
                    % (ewH[-1], parCond))

                # sift through the parents and purge it of 'ideling' individuals
                # which reduce the fitness

                D0 = Civilizations[-1][3]

                parBV = D0[-1][0]

                Ais = copy.deepcopy(Civilizations[-1])
                print(Ais[-2])
                exit()
                # generate a pair of offspring from a randomly selected couple
                newBV = True
                growthType = 'add an additional CFG'

                print('\nciv-%d)  growthType: %s' % (nCivi, growthType))

                # bv with new int an relw. => append new offspring cfg's

                if newBV:

                    anzBV = sum([len(iws) for iws in Ais[1]])

                    childishParaSets = []
                    childishParaSets.append([
                        Ais[3], Jay, costr, zop, 10, [-1], BINBDGpath, minCond,
                        denseEVinterval
                    ])
                    chiBV = nbrOff

                    # produce an offspring cfg for each parent cfg
                    for cfg in range(len(Ais[0])):

                        childidentifier = [cfg]

                        tmpBas = copy.deepcopy(Ais[3])

                        # from the iw's of the parent cfg, select mother/father pairs
                        iwpairs = [
                            ip for ip in list(
                                product(range(len(Ais[1][cfg])), repeat=2))
                            if ip[0] != ip[1]
                        ]
                        np.random.shuffle(iwpairs)
                        iwpairs = iwpairs[:int((len(Ais[1][cfg])) / 2)]

                        #Ais[0].append(Ais[0][cfg])
                        Ais[2].append([])

                        daughterson = []
                        for iws in iwpairs:
                            mother = Ais[1][cfg][iws[0]]
                            father = Ais[1][cfg][iws[1]]
                            daughterson.append(
                                intertwining(mother,
                                             father,
                                             mutation_rate=muta_initial))
                            rwshake1 = 0.9 + 0.2 * np.random.random()
                            rwshake2 = 0.9 + 0.2 * np.random.random()
                            rwa = list(
                                np.random.choice(rwshake1 *
                                                 np.array(Ais[2][cfg][iws[0]]),
                                                 2,
                                                 replace=False))
                            rwa.sort()
                            rw1 = rwa[::-1]
                            rwa = list(
                                np.random.choice(rwshake2 *
                                                 np.array(Ais[2][cfg][iws[1]]),
                                                 2,
                                                 replace=False))
                            rwa.sort()
                            rw2 = rwa[::-1]
                            Ais[2][-1].append(rw1)
                            anzBV += 1
                            tmpBas.append(
                                [anzBV, list(range(1, 1 + len(rw1)))])
                            Ais[2][-1].append(rw2)
                            anzBV += 1
                            tmpBas.append(
                                [anzBV, list(range(1, 1 + len(rw2)))])

                        dstmp = np.array(daughterson).flatten()
                        dstmp.sort()
                        daughterson = list(dstmp)[::-1]

                        Ais[1].append(daughterson)

                        childishParaSets.append([
                            tmpBas, Jay, costr, zop, 10, childidentifier,
                            BINBDGpath, minCond, denseEVinterval
                        ])

                # 2) calc. matrices including all children
                print(Ais[0])
                print(childishParaSets[0])
                exit()

                # x) the parallel environment is set up in sets(chunks) of bases
                #    in order to limit the number of files open simultaneously
                split_points = [
                    n * maxParLen
                    for n in range(1 + int(len(childishParaSets) / maxParLen))
                ] + [len(childishParaSets) + 42]

                Parchunks = [
                    childishParaSets[split_points[i]:split_points[i + 1]]
                    for i in range(len(split_points) - 1)
                ]

                pipe_list = []

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
                        pipe_list.append(recv_end)
                        p.start()
                    for proc in jobs:
                        proc.join()

                child_ladder = [x.recv() for x in pipe_list]
                # ranking (0=cond; 1=fit; 2=bind)
                print('\n>>> Offspring order based on %s <<<' % pur[pw])
                child_ladder.sort(key=lambda tup: np.abs(tup[pw]))
                child_ladder.reverse()

                for ch in child_ladder[:3]:
                    print(ch[:4])

                if child_ladder[0][0] > 0.0:

                    if newBV:

                        nbrParBV = sum(
                            [len(iw) for iw in Civilizations[-1][1]])
                        parvenue = child_ladder[0]

                        # check if sexiest basis includes an offspring; if [3][0]==-1
                        # the most glamorous basis is comprised of the parents, only!
                        if int(parvenue[3][0]) >= 0:
                            Aisopt = copy.deepcopy(Civilizations[-1])

                            nbrParentCFGs = len(Civilizations[-1][1])
                            optCFG = int(parvenue[3][0])
                            Aisopt[0].append(Civilizations[-1][0][optCFG])

                            optCHiw = nbrParentCFGs + int(parvenue[3][0])

                            Aisopt[1].append(Ais[1][optCHiw])
                            Aisopt[2].append(Ais[2][optCHiw])

                            newBVs = [[
                                nbrParBV + 1 + i,
                                list(range(1, 1 + len(Ais[2][optCHiw][i])))
                            ] for i in range(len(Ais[1][optCHiw]))]
                            for nbv in newBVs:
                                Aisopt[3].append(nbv)

                            stratifiedOptCivilization = Aisopt
                        else:
                            print('no parvenue! the old generation prevails.')
                            stratifiedOptCivilization = Civilizations[-1]

                    elif newRW:
                        Aisopt = copy.deepcopy(Civilizations[-1])

                        parvenue = child_ladder[0][3]
                        # Aisopt[0,1] do not change
                        # the basis is improved by a set of relative widths for a specific configuration

                        for neww in range(len(parvenue[:-1])):

                            Aisopt[2][parvenue[-1]][neww].append(
                                Ais[2][parvenue[-1]][neww][parvenue[neww] - 1])

                            for nbv in range(len(Aisopt[3])):
                                if (Aisopt[3][nbv][0] -
                                        1 == (cfgbnds[parvenue[-1]] + neww)):
                                    Aisopt[3][nbv][1].append(
                                        len(Aisopt[2][parvenue[-1]][neww]))

                        stratifiedOptCivilization = Aisopt

                    Civilizations.append(stratifiedOptCivilization)

                    print(
                        'civ-%d) a new, %d-dimensional race has evolved!\n -----------------------------'
                        % (nCivi, basisDim(stratifiedOptCivilization[3])))

                else:
                    print(
                        'civ-%d) only spoiled offspring was raised!\n       >>> increased mutation rate.\n -----------------------------'
                        % nCivi)
                    muta_initial *= 1.1
                    Civilizations.append(Civilizations[-1])
                    if muta_initial >= 1.0:
                        print(
                            '... although wild mutations were allowed. Aborting...'
                        )
                        break

        ma = blunt_ev(Civilizations[-1][0],
                      Civilizations[-1][1],
                      Civilizations[-1][2],
                      Civilizations[-1][3],
                      wrkdir='',
                      nzopt=zop,
                      costring=costr,
                      bin_path=BINBDGpath,
                      mpipath=MPIRUN,
                      einzel_file_path=wrkDir,
                      potNN=potnn,
                      potNNN=potnnn,
                      parall=-1,
                      anzcores=max(2, min(len(Civilizations[-1][0]), MaxProc)),
                      tnni=10,
                      jay=Jay,
                      dia=True)

        ewN, ewH = NormHamDiag(ma)

        print(
            'After %d generations,\n-----\n' % nCivi, Civilizations[-1],
            '\n-----\nemerged as the dominant culture.\nDim = %d) B(GS) = %8.4f  fit = '
            % (basisDim(Civilizations[-1][3]), ewH[-1]),
            basQ(ewN, ewH, minCond, denseEVinterval))

        suf = 'ref' if bastype == boundstatekanal else 'fin'

        lfrags = np.array(Civilizations[-1][0])[:, 1].tolist()
        sfrags = np.array(Civilizations[-1][0])[:, 0].tolist()
        n3_inlu(8, fn=basisPath + 'INLU_%s' % suf, fr=lfrags, indep=-1)
        n3_inlu(8, fn=basisPath + 'INLUCN_%s' % suf, fr=lfrags, indep=-1)
        n3_inob(sfrags, 8, fn=basisPath + 'INOB_%s' % suf, indep=-1)
        n3_inob(sfrags, 15, fn=basisPath + 'DRINOB_%s' % suf, indep=-1)

        os.system('cp INQUA_M ' + basisPath + 'INQUA_V18_%s' % suf)
        os.system('cp INEN ' + basisPath + 'INEN_%s' % suf)
        os.system('cp INSAM ' + basisPath)

        os.system('rm -rf ./inen_*')
        os.system('rm -rf ./endout_*')
        os.system('rm -rf ./MATOUTB_*')
        os.system('rm -rf ./T*OUT.*')
        os.system('rm -rf ./D*OUT.*')

        fullBasfile, actBasfile = write_basis_on_tape(Civilizations[-1],
                                                      Jay,
                                                      bastype,
                                                      baspath=basisPath)

        if bastype != boundstatekanal:
            AbasOutStr = respath + 'Ssigbasv3heLIT_%s_BasNR-%d.dat' % (bastype,
                                                                       nB)
            FbasOutStr = respath + 'SLITbas_full_%s_BasNR-%d.dat' % (bastype,
                                                                     nB)
            subprocess.call('cp %s %s' % (fullBasfile, FbasOutStr), shell=True)
            subprocess.call('cp %s %s' % (actBasfile, AbasOutStr), shell=True)
        else:
            AbasOutStr = respath + 'Ssigbasv3heLIT_%s.dat' % bastype
            FbasOutStr = respath + 'SLITbas_full_%s.dat' % bastype
            subprocess.call('cp %s %s' % (fullBasfile, FbasOutStr), shell=True)
            subprocess.call('cp %s %s' % (actBasfile, AbasOutStr), shell=True)

        matoutstr = '%smat_%s' % (
            respath, bastype + '_BasNR-%d' % nB
        ) if bastype != boundstatekanal else '%smat_%s' % (respath, bastype)
        subprocess.call('cp MATOUTB %s' % matoutstr, shell=True)
        print(
            'channel %s: Basis structure, Norm, and Hamiltonian written into %s'
            % (bastype, respath + 'mat_' + bastype))

        # for the bound-state/initial-state channel, consider only one basis set
        if bastype == boundstatekanal:
            break