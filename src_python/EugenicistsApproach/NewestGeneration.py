from bridgeA3 import *
from genetic_width_growth import *

import glob, copy
import shlex
import multiprocessing
from multiprocessing.pool import ThreadPool
from itertools import permutations, product

os.chdir(litpath3He)

#bastypes = streukas  #+[boundstatekanal]

for bastype in bastypes:

    print('\n#---------   Basistype: %s ---------------#\n' % bastype)
    # 1) calculation for ONE trail channel, only.
    angu = channels[bastype]
    Jstreu = float(bastype.split('^')[0][-3:])
    Jstreustring = '%s' % str(Jstreu)[:3]

    # 2) read the initial population of parents, e.g.,
    # - geom. width sets         (intwLIT,rws)
    # - (iso)spin configurations (cfgs)
    # - unstable, "full" basis   (Dfull)

    # Read first generation of parents == (B0)
    cfgs = [
        con.split() for con in open(basisPath + 'frags_LIT_J%s_%s.dat' %
                                    (Jstreustring, bastype))
    ]
    origCFGs = copy.deepcopy(cfgs)

    intwLIT = [
        np.array(ln.split()).astype(float).tolist()
        for ln in open(basisPath + 'intw3heLIT_J%s_%s.dat' %
                       (Jstreustring, bastype))
    ]
    relwLIT = [
        np.array(ln.split()).astype(float).tolist()
        for ln in open(basisPath + 'relw3heLIT_J%s_%s.dat' %
                       (Jstreustring, bastype))
    ]
    rws = []
    rw0 = 0
    for cfg in range(len(intwLIT)):
        rws.append([])
        for bv in range(len(intwLIT[cfg])):
            rws[-1].append(relwLIT[bv + rw0])
        rw0 += len(intwLIT[cfg])

    minDiffwidthsINT = 10**-2
    minDiffwidthsREL = 10**-3

    minCond = 10**-11
    maxE = 100.0

    maxParLen = 18
    removalGainFactor = 1.5
    muta_initial = 0.92
    # nRaces := |i|
    nRaces = 3 if bastype == boundstatekanal else 20
    nbrOff = 2
    MaxOff = 12

    dbg = False

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
                np.array(range(
                    1, 1 + len(Civilizations[0][2][cfg][nbvc - 1]))).tolist()
            ]]

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
                      einzel_file_path=v18uixpath,
                      potNN=potnn,
                      potNNN=potnnn,
                      parall=-1,
                      tnni=10,
                      jay=Jstreu,
                      dia=False)

        ewN, ewH = NormHamDiag(ma)

        #print('parents: \n', Civilizations[-1])

        if ewH == []:
            print('parent basis unstable:\n', Civilizations[-1][3])
            exit()

        print(
            'civ-%d)          parents: Dim = %d) B(GS) = %8.4f  fit = ' %
            (nCivi, basisDim(Civilizations[-1][3]), ewH[-1]), basQ(ewN, ewH))

        # sift through the parents and purge it of 'ideling' individuals
        # which yield reduce the fitness

        D0 = Civilizations[-1][3]

        # according to the quality criterion <pw> the offspring is ordered
        pw = 1 + nCivi % 2
        # for each civilization, the initial population of parents is purged, i.e.,
        # parents which destabilize the basis significantly are removed
        pwpurge = 0

        pur = [
            'condition number', 'quality = f(#EV<100, B(GS), cond. nbr.)',
            'B(GS)'
        ]

        print('              commencing purges (%s) ' % (pur[pwpurge]),
              end='\n')

        goPurge = True

        while goPurge:

            newpopList = []
            goPurge = False
            ParaSets = []

            D0flat = flatten_basis(D0)

            cpy = copy.deepcopy(D0flat)
            ParaSets.append([
                cpy, Jstreu, costr, zop, 10, [0, 0], BINBDGpath, minCond, maxE
            ])

            for bvTrail in D0flat:

                bvID = [int(bvTrail[0]), int(''.join(map(str, bvTrail[1])))]

                cpy = copy.deepcopy(D0flat)

                cpy.remove(bvTrail)
                cpy = rectify_basis(cpy)

                ParaSets.append([
                    cpy, Jstreu, costr, zop, 10, bvID, BINBDGpath, minCond,
                    maxE
                ])

            # x) the parallel environment is set up in sets(chunks) of bases
            #    in order to limit the number of files open simultaneously
            split_points = [
                n * maxParLen
                for n in range(1 + int(len(ParaSets) / maxParLen))
            ] + [len(ParaSets) + 42]

            Parchunks = [
                ParaSets[split_points[i]:split_points[i + 1]]
                for i in range(len(split_points) - 1)
            ]

            cand_list = []
            for chunk in Parchunks:
                pool = ThreadPool(anzproc)
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

            print('reference (cond,qual) = %4.4e,%4.4e' % (reff[0], reff[1]))

            print('best/worst=  %2.3e , %2.3e /  %2.3e , %2.3e' %
                  (cand_ladder[-1][0], cand_ladder[-1][1], cand_ladder[1][0],
                   cand_ladder[1][1]))

            if ((removalGainFactor * np.abs(reff[pwpurge]) < np.abs(
                    cand_ladder[-1][pwpurge])) |
                ((np.abs(cand_ladder[-1][pwpurge]) < minCond) &
                 (reff[0] < minCond))):
                goPurge = True
                D0 = rectify_basis(cand_ladder[-1][4])
                print('1/%d ' % basisDim(D0), end='\n')

        removalGainFactor *= 1.1

        #        D0 = rectify_basis(D0flat)

        # -- end of purges; the removal of any BV will now reduce the population's fitness

        n3_inen_bdg(D0, Jstreu, costr, fn='INEN', pari=0, nzop=zop, tni=tnni)

        subprocess.run([BINBDGpath + 'TDR2END_NORMAL.exe'],
                       capture_output=True,
                       text=True)

        matout = np.core.records.fromfile('MATOUTB', formats='f8', offset=4)

        ewN, ewH = NormHamDiag(matout)

        Civilizations[-1][3] = D0

        parBV = D0[-1][0]

        Ais = copy.deepcopy(Civilizations[-1])

        # generate a pair of offspring from a randomly selected couple
        newRW = True if nCivi % 2 != 0 else False
        newBV = True if nCivi % 2 == 0 else False
        growthType = 'relative widths to a CFG' if newRW else 'an additional CFG'

        print(
            '\nciv-%d) (purged) parents: Dim = %d) B(GS) = %8.4f  fit = %8.4f\n        >>> adding %s <<<'
            %
            (nCivi, basisDim(D0), ewH[-1],
             basQ(ewN, ewH, denseEnergyInterval=[-1000, maxE])[0], growthType))

        bvsPerCfg = 12
        # bv with new int an relw. => append new offspring cfg's

        if newBV:

            anzBV = sum([len(iws) for iws in Ais[1]])

            childishParaSets = []
            childishParaSets.append([
                Ais[3], Jstreu, costr, zop, 10, [-1], BINBDGpath, minCond, maxE
            ])
            chiBV = nbrOff

            # produce an offspring cfg for each parent cfg
            for cfg in range(len(Ais[0])):

                childidentifier = [cfg]

                tmpBas = copy.deepcopy(Ais[3])

                # from the iw's of the parent cfg, select mother/father pairs
                iwpairs = [
                    ip
                    for ip in list(product(range(len(Ais[1][cfg])), repeat=2))
                    if ip[0] != ip[1]
                ]
                np.random.shuffle(iwpairs)
                iwpairs = iwpairs[:int((len(Ais[1][cfg])) / 2)]

                Ais[0].append(Ais[0][cfg])
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
                    tmpBas.append([anzBV, list(range(1, 1 + len(rw1)))])
                    Ais[2][-1].append(rw2)
                    anzBV += 1
                    tmpBas.append([anzBV, list(range(1, 1 + len(rw2)))])

                dstmp = np.array(daughterson).flatten()
                dstmp.sort()
                daughterson = list(dstmp)[::-1]

                Ais[1].append(daughterson)

                childishParaSets.append([
                    tmpBas, Jstreu, costr, zop, 10, childidentifier,
                    BINBDGpath, minCond, maxE
                ])

        # the internal width is kept fixed, and the offspring expands the relative-width set of an existing BV
        if newRW:

            childishParaSets = []
            chiBV = nbrOff
            for cfg in range(len(Ais[0])):

                chcombos = product(range(
                    len(Ais[2][cfg][0]) + 1,
                    len(Ais[2][cfg][0]) + nbrOff + 1),
                                   repeat=len(Ais[1][cfg]))

                tmp = list(chcombos)
                np.random.shuffle(tmp)
                offgenMax = min(len(tmp), MaxOff)
                chcombos = tmp[:offgenMax]

                for combo in chcombos:

                    childidentifier = list(combo)
                    childidentifier.append(cfg)

                    tmpBas = copy.deepcopy(Ais[3])
                    for bvn in range(len(Ais[1][cfg])):

                        try:
                            for bvm in range(len(tmpBas)):
                                if (tmpBas[bvm][0] == (cfgbnds[cfg] + bvn +
                                                       1)):
                                    tmpBas[bvm] = [
                                        tmpBas[bvm][0],
                                        tmpBas[bvm][1] + [list(combo)[bvn]]
                                    ]
                        except:
                            print(cfg, bvn, list(combo)[bvn])
                            print(Ais[3])
                            print(cfgbnds[cfg])
                            print(tmpBas)

                    childishParaSets.append([
                        tmpBas, Jstreu, costr, zop, 10, childidentifier,
                        BINBDGpath, minCond, maxE
                    ])
                    chiBV += 1

                for bvn in range(len(Ais[1][cfg])):
                    l0 = len(Ais[2][cfg][bvn])
                    while len(Ais[2][cfg][bvn]) < l0 + nbrOff:
                        motherfather = np.random.choice(l0, 2, replace=False)
                        childRW = intertwining(
                            Ais[2][cfg][bvn][motherfather[0]],
                            Ais[2][cfg][bvn][motherfather[1]],
                            mutation_rate=0.2)
                        for ca in childRW:
                            if min([
                                    np.abs(ca - rw) for rw in Ais[2][cfg][bvn]
                            ]) > minDiffwidthsREL:
                                Ais[2][cfg][bvn].append(ca)
                                break
                    tmp = Ais[2][cfg][bvn]
                    Ais[2][cfg][bvn] = tmp

        # 2) calc. matrices including all children
        ma = blunt_ev(Ais[0],
                      Ais[1],
                      Ais[2],
                      Ais[3],
                      wrkdir='',
                      nzopt=zop,
                      costring=costr,
                      bin_path=BINBDGpath,
                      einzel_file_path=v18uixpath,
                      potNN=potnn,
                      potNNN=potnnn,
                      parall=-1,
                      tnni=10,
                      jay=Jstreu,
                      dia=False)

        # 3) rate certain partitions according to their fitness
        if glob.glob('MATOUTB*') != []:
            os.system('rm -rf MATOUTB*')

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
            pool = ThreadPool(anzproc)
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

                nbrParBV = sum([len(iw) for iw in Civilizations[-1][1]])
                parvenue = child_ladder[0]

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

            if newRW:
                Aisopt = copy.deepcopy(Civilizations[-1])

                parvenue = child_ladder[0][3]
                # Aisopt[0,1] do not change
                # the basis is improved by a set of relative widths for a specific configuration

                for neww in range(len(parvenue[:-1])):

                    Aisopt[2][parvenue[-1]][neww].append(
                        Ais[2][parvenue[-1]][neww][parvenue[neww] - 1])

                    for nbv in range(len(Aisopt[3])):
                        if (Aisopt[3][nbv][0] - 1 == (cfgbnds[parvenue[-1]] +
                                                      neww)):
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
                print('... although wild mutations were allowed. Aborting...')
                break

    ma = blunt_ev(Civilizations[-1][0],
                  Civilizations[-1][1],
                  Civilizations[-1][2],
                  Civilizations[-1][3],
                  wrkdir='',
                  nzopt=zop,
                  costring=costr,
                  bin_path=BINBDGpath,
                  einzel_file_path=v18uixpath,
                  potNN=potnn,
                  potNNN=potnnn,
                  parall=-1,
                  tnni=10,
                  jay=Jstreu,
                  dia=True)
    ewN, ewH = NormHamDiag(ma)

    print(
        'After %d generations,\n-----\n' % nCivi, Civilizations[-1],
        '\n-----\nemerged as the dominant culture.\nDim = %d) B(GS) = %8.4f  fit = '
        % (basisDim(Civilizations[-1][3]), ewH[-1]), basQ(ewN, ewH))

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

    write_basis_on_tape(Civilizations[-1], Jstreu, bastype, baspath=basisPath)
    subprocess.call('cp MATOUTB %smat_%s' % (respath, bastype), shell=True)
    print(
        'channel %s: Basis structure, Norm, and Hamiltonian written into %s' %
        (bastype, respath + 'mat_' + bastype))