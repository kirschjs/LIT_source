from bridgeA3 import *
from genetic_width_growth import *

os.chdir(litpath3He)

basisPath = litpath3He + 'basis_struct/'
bastypes = [boundstatekanal]  # + streukas

# 1) calculation for ONE trail channel, only.
angu = channels[bastypes[0]]
Jstreu = float(bastypes[0].split('^')[0][-3:])
Jstreustring = '%s' % str(Jstreu)[:3]

# 2) read the initial population of parents, e.g.,
# - geom. width sets         (intwLIT,rws)
# - (iso)spin configurations (cfgs)
# - unstable, "full" basis   (Dfull)

nGen = 0
# Read first generation of parents == (B0)
cfgs = [
    con.split() for con in open(basisPath + 'frags_LIT_J%s_%s.dat' %
                                (Jstreustring, bastypes[0]))
]
intwLIT = [
    np.array(ln.split()).astype(float).tolist()
    for ln in open(basisPath + 'intw3heLIT_J%s_%s.dat' %
                   (Jstreustring, bastypes[0]))
]
relwLIT = [
    np.array(ln.split()).astype(float).tolist()
    for ln in open(basisPath + 'relw3heLIT_J%s_%s.dat' %
                   (Jstreustring, bastypes[0]))
]
rws = []
rw0 = 0
for cfg in range(len(intwLIT)):
    rws.append([])
    for bv in range(len(intwLIT[cfg])):
        rws[-1].append(relwLIT[bv + rw0])
    rw0 += len(intwLIT[cfg])
dbg = False
ewMax = 1000
minDiffwidthsINT = 10**-2
minDiffwidthsREL = 10**-3
mDiff = 10**-5
minCond = 10**-6
nGenMax = 2
# nRaces := |i|
nRaces = 12
# nBasisSamples := |j|
nBasisSamples = 4
nbrOff = 2
targetDimfac0 = 0.95

Civilizations = []

# 3) breed the first offspring generation
# each pair of basis vectors procreates => |parents + offspring| = 2*|parents|
# P_0 -> {P_0^i=parents+offspring_i,i=1,...,nRaces}

P0 = [cfgs, intwLIT, rws, []]

for nCivi in range(nRaces):

    if nCivi == 0:
        # i) for each (iso)spin configuration, select pairs to procreate
        #    and breed offsping generation == (C0)
        iTNG = []
        rTNG = []
        for cfg in range(len(P0[1])):
            TNGt = []
            # in case of a one-person population, carve an Eve from the rip of Adam
            if len(P0[1][cfg]) == 1:
                P0[1][cfg] += [P0[1][cfg][0] * np.random.random()]
                P0[2][cfg] += [
                    (np.array(P0[2][cfg][0]) * np.random.random()).tolist()
                ]
            inxs = np.arange(len(P0[1][cfg]))
            np.random.shuffle(inxs)
            parentpairs = np.reshape(
                inxs, (-1, 2)) if len(inxs) % 2 == 0 else np.reshape(
                    np.append(inxs, np.random.choice(inxs[:-1])), (-1, 2))
            for coupl in parentpairs:
                muta_rate = 0.1
                father = P0[1][cfg][coupl[0]]
                mother = P0[1][cfg][coupl[1]]
                son, daughter = intertwining(father,
                                             mother,
                                             mutation_rate=muta_rate)
                gogo = True
                while ((np.min(
                        np.array([
                            np.abs(float(father) - float(daughter)),
                            np.abs(float(father) - float(son)),
                            np.abs(float(mother) - float(son)),
                            np.abs(float(mother) - float(daughter)),
                            np.abs(float(son) - float(daughter))
                        ])) < minDiffwidthsINT) & (gogo)):
                    muta_rate += 0.01
                    if muta_rate >= 1:
                        muta_rate = 0.1
                        gogo = False
                    son, daughter = intertwining(father,
                                                 mother,
                                                 mutation_rate=muta_rate)
                    assert ((np.isnan(son) == False) &
                            (np.isnan(daughter) == False))
                TNGt.append([son, daughter])
            if TNGt == []:
                print('no offpring...')
                continue
            iTNG.append(
                np.sort(
                    np.random.choice(
                        np.ravel(TNGt), size=len(P0[1][cfg]),
                        replace=False))[::-1].astype(float).tolist())
            rTNG.append([])
            for relwS in range(len(P0[1][cfg])):
                # in case of a one-person population, carve an Eve from the rip of Adam
                if len(P0[2][cfg][relwS]) == 1:
                    P0[2][cfg][relwS] += [
                        P0[2][cfg][relwS][0] * np.random.random()
                    ]
                inxs = np.arange(len(P0[2][cfg][relwS]))
                np.random.shuffle(inxs)
                parentpairs = np.reshape(
                    inxs, (-1, 2)) if len(inxs) % 2 == 0 else np.reshape(
                        np.append(inxs, np.random.choice(inxs[:-1])), (-1, 2))
                TNGt = []
                for coupl in parentpairs:
                    muta_rate = 0.1
                    father = P0[2][cfg][relwS][coupl[0]]
                    mother = P0[2][cfg][relwS][coupl[1]]
                    son, daughter = intertwining(father,
                                                 mother,
                                                 mutation_rate=muta_rate)
                    gogo = True
                    while ((np.min(
                            np.array([
                                np.abs(float(father) - float(daughter)),
                                np.abs(float(father) - float(son)),
                                np.abs(float(mother) - float(son)),
                                np.abs(float(mother) - float(daughter)),
                                np.abs(float(son) - float(daughter))
                            ])) < minDiffwidthsREL) & (gogo == True)):
                        muta_rate += 0.01
                        if muta_rate >= 1:
                            muta_rate = 0.1
                            gogo = False
                        son, daughter = intertwining(father,
                                                     mother,
                                                     mutation_rate=muta_rate)
                        assert ((np.isnan(son) == False) &
                                (np.isnan(daughter) == False))
                    TNGt.append([son, daughter])
                rTNG[-1].append(
                    np.sort(np.ravel(TNGt))[::-1].astype(float).tolist())
            if iTNG == []:
                print('no offspring...')
                exit()
            iwTNG = P0[1] + iTNG
            rwTNG = P0[2] + rTNG

        for cfg in range(len(rwTNG)):
            anzrelmin = np.max([len(relS) for relS in rwTNG[cfg]])
            for relwS in range(len(rwTNG[cfg])):
                if len(rwTNG[cfg][relwS]) < anzrelmin:
                    rwTNG[cfg][relwS] = np.linspace(np.min(rwTNG[cfg][relwS]),
                                                    np.max(rwTNG[cfg][relwS]),
                                                    anzrelmin)[::-1]
        if dbg:
            print('IW(parents)  :\n', P0[1])
            print('IW(offspring):\n', iTNG, '\n')
            print('IW(family)   :\n', iwTNG, '\nRW(family)   :', rwTNG)

        P0 = [2 * P0[0], iwTNG, rwTNG, []]
        Civilizations.append(P0)

    Dfull = []
    anzPro = 0
    for cfg in range(len(Civilizations[nCivi][1])):
        for bv in range(len(Civilizations[nCivi][1][cfg])):
            anzPro += 1
            tnp = np.arange(1,
                            len(Civilizations[nCivi][2][cfg][bv]) +
                            1).tolist()
            Dfull.append([anzPro, tnp])
    if Civilizations[nCivi][3] == []:
        Civilizations[nCivi][3] = Dfull

    cfgbnds = np.add.accumulate([len(iws) for iws in Civilizations[nCivi][1]])
    cfgbnds = np.insert(cfgbnds, 0, 0)

    lfragTNG = np.array(Civilizations[nCivi][0])[:, 1].tolist()
    sfragTNG = np.array(Civilizations[nCivi][0])[:, 0].tolist()
    insam(len(lfragTNG))

    # i) each initial parent-offspring basis is likely unstable if comprised of all 'families'(Dfull)

    n3_inlu(8, fn='INLU', fr=lfragTNG, indep=parall)
    os.system(BINBDGpath + 'DRLUD.exe')
    n3_inlu(8, fn='INLUCN', fr=lfragTNG, indep=parall)
    os.system(BINBDGpath + 'LUDW_CN.exe')
    n3_inob(sfragTNG, 8, fn='INOB', indep=parall)
    os.system(BINBDGpath + 'KOBER.exe')
    n3_inob(sfragTNG, 15, fn='INOB', indep=parall)
    os.system(BINBDGpath + 'DROBER.exe')

    he3inquaBS(intwi=Civilizations[-1][1],
               relwi=Civilizations[-1][2],
               potf=potnn,
               inquaout='INQUA_M_0')
    parallel_mod_of_3inqua(lfragTNG,
                           sfragTNG,
                           infile='INQUA_M_0',
                           outfile='INQUA_M',
                           einzel_path=v18uixpath)
    subprocess.call('cp INQUA_M inq_tmp', shell=True)
    if parall == -1:
        subprocess.run([
            'mpirun', '-np',
            '%d' % anzproc, BINBDGpath + 'V18_PAR/mpi_quaf_v7'
        ])
        subprocess.run([BINBDGpath + 'V18_PAR/sammel'])

    if tnni == 11:
        he3inquaBS(intwi=Civilizations[-1][1],
                   relwi=Civilizations[-1][2],
                   potf=potnnn)
        parallel_mod_of_3inqua(lfragTNG,
                               sfragTNG,
                               infile='INQUA_M',
                               outfile='INQUA_M',
                               tni=1,
                               einzel_path=v18uixpath)
        if parall == -1:
            subprocess.run([
                'mpirun', '-np',
                '%d' % anzproc, BINBDGpath + 'UIX_PAR/mpi_drqua_v7'
            ])
            subprocess.run([BINBDGpath + 'UIX_PAR/SAMMEL-uix'])
        else:
            subprocess.run([BINBDGpath + 'DRQUA_AK_M.exe'])

    if dbg:
        print('all qua-MEs calculated. Entering individual BV assessment.')

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

    if ewH == []:
        print('parent basis unstable:\n', Civilizations[-1][3])

        exit()

    basCond = np.min(np.abs(ewN)) / np.max(np.abs(ewN)) if ewH != [] else 0.

    print('civ-%d) initial stability (condition number): ' % nCivi, basCond)
    exit()

    # generate stable-basis candidates for the offspring set nCivi
    # < A_i -> A^j_i >
    # initial j-basis dimension which might still result in an entirely unstable
    # set of bases; if so, the targetDim is reduced until a stable set is found
    targetDim = int(basisDim(Civilizations[-1][3]) * targetDimfac0)

    go = True if basCond < minCond else False

    while (go):
        # select the candidates for stable bases randomly from the full, but unstable basis
        competingCivilizations = [
            select_random_basis(Civilizations[-1][3], targetDim)
            for n in range(nBasisSamples)
        ]

        oAijList = []

        print(
            'stabilizing %d-dimensional basis with %d random %d-dim samples' %
            (basisDim(Civilizations[-1][3]), nBasisSamples, targetDim))

        for oAij in competingCivilizations:
            # assess oAij's resilience
            n3_inen_bdg(oAij,
                        Jstreu,
                        costr,
                        fn='INEN',
                        pari=0,
                        nzop=zop,
                        tni=tnni)
            if parall == -1:
                subprocess.run([BINBDGpath + 'TDR2END_NORMAL.exe'],
                               capture_output=True,
                               text=True)
            else:
                subprocess.run([BINBDGpath + 'DR2END_NORMAL.exe'])
            suche_fehler()
            matout = np.core.records.fromfile('MATOUTB',
                                              formats='f8',
                                              offset=4)
            specN, specH = NormHamDiag(matout)
            if specH != []:
                Quala = specH[-1]
                Qualb = len([bvv for bvv in specN if bvv < ewMax]) + len(
                    [bvv for bvv in specH if bvv < ewMax])
                Cond = np.min(np.abs(specN)) / np.max(np.abs(specN))
            else:
                print(
                    'eigensystem calculation failed for rnd basis sample of people set %d'
                    % nCivi)
                Quala = 10**4
                Qualb = 0
                Cond = 0.
            # admit the random basis of the norm's condition number exceeds a lower bound
            if Cond > minCond:
                oAijList.append([oAij, Quala, Qualb, Cond])
        if oAijList != []:
            # ECCE: the list of oAij's comprises only stable elements (condNbr > minCond)
            # hence, from these, we pick the ``strongest''!
            nQual = 1
            idx = np.array([elem[nQual] for elem in oAijList]).argsort()[::-1]
            oAijList = [eww for eww in np.array(oAijList)[idx]]
            # reset the candidate population to the strongest and most stable Doi element
            Dopt = oAijList[::-1][-1]
            go = False
        else:
            # if none of the random bases is stable, take a new, lower-dimensional sample
            if dbg:
                print('all bases unstable. reducing dimension!')
            targetDim = int(targetDim * targetDimfac0)
            assert targetDim > 3

    if dbg & (basCond < minCond):
        print(
            'people set (i) = %d chose the fittest of |j| = %d samples:\n' %
            (nCivi, len(oAijList)), 'quality:', Dopt[1:],
            '\n dimensionality: %d/%d' %
            (basisDim(Dopt[0]), basisDim(Civilizations[-1][3])))

    # sift through the candidate population and purge it of 'ideling' individuals
    # which are those without whom the ground state's energy drops by less than some epsilon > 0
    D0 = Dopt[0] if basCond < minCond else Civilizations[-1][3]
    go = True
    nPW = 0

    while go:

        if dbg:
            print('testing removal of 1/%d basis vectors' % basisDim(D0))

        nPW += 1
        newpopList = []
        go = False
        D0flat = flatten_basis(D0)
        for bvTrail in ['Raeuber'] + D0flat:
            cpy = D0flat.copy()
            if bvTrail in cpy:
                cpy.remove(bvTrail)
            n3_inen_bdg(cpy,
                        Jstreu,
                        costr,
                        fn='INEN',
                        pari=0,
                        nzop=zop,
                        tni=tnni)
            subprocess.run([BINBDGpath + 'TDR2END_NORMAL.exe'],
                           capture_output=True,
                           text=True)
            matout = np.core.records.fromfile('MATOUTB',
                                              formats='f8',
                                              offset=4)
            specN, specH = NormHamDiag(matout)
            # number of norm AND ham ev < TH
            if bvTrail == 'Raeuber':
                Quala = specH[-1]
                continue

            if specH != []:
                basQuala = specH[-1]
                basQualb = len([bvv for bvv in specN if bvv < ewMax]) + len(
                    [bvv for bvv in specH if bvv < ewMax])
            else:
                basQuala = 10**4
                basQualb = 0

            if np.abs(basQuala - Quala) < mDiff:
                newpopList.append([
                    cpy,
                    np.abs(basQuala - Quala), basQualb,
                    np.min(np.abs(specN)) / np.max(np.abs(specN))
                ])
        if newpopList != []:
            go = True
            nQual = 1
            idx = np.array([elem[nQual]
                            for elem in newpopList]).argsort()[::-1]
            newpopList = [eww for eww in np.array(newpopList)[idx]]
            D0 = rectify_basis(newpopList[0][0])

    D0 = rectify_basis(D0flat)
    n3_inen_bdg(D0, Jstreu, costr, fn='INEN', pari=0, nzop=zop, tni=tnni)
    if parall == -1:
        subprocess.run([BINBDGpath + 'TDR2END_NORMAL.exe'],
                       capture_output=True,
                       text=True)
    else:
        subprocess.run([BINBDGpath + 'DR2END_NORMAL.exe'])

    matout = np.core.records.fromfile('MATOUTB', formats='f8', offset=4)

    ewN, ewH = NormHamDiag(matout)

    if ewH != []:
        basQuala_i = ewH[-1]
        basQualb_i = len([bvv for bvv in ewN if bvv < ewMax]) + len(
            [bvv for bvv in ewH if bvv < ewMax])
        basCond_i = np.min(np.abs(ewN)) / np.max(np.abs(ewN))
    else:
        basQuala_i = 10**3
        basQualb_i = 0
        basCond_i = 10**-19

    if dbg:
        print('D0 bare:', D0, '\n')

    print(
        'B0(gen %d/%d, civ %d/%d, basDim = %d) = ' %
        (nGen + 1, nGenMax, nCivi + 1, nRaces, basisDim(D0)), ewH[-1])

    if nCivi > 2:
        exit()

    # write a stable and fit population (i) on TAPE
    bvcfgs = [[
        ncfg for ncfg in range(len(cfgbnds) - 1)
        if ((bv[0] > cfgbnds[ncfg]) & (bv[0] <= cfgbnds[ncfg + 1]))
    ] for bv in D0]
    bvcfgnames = (np.array(2 * Civilizations[-1][0])[sum(bvcfgs, [])].tolist())
    D0 = np.concatenate((np.array(D0), bvcfgnames, bvcfgs), axis=-1).tolist()
    D1 = {}
    for cfg in np.unique(bvcfgnames, axis=0).tolist():
        D1[str(cfg)] = []
        for bv in D0:
            if bv[2:4] == cfg:
                D1[str(cfg)].append(bv[:2] + [bv[-1]])
    rwTNGt = []
    for nn in range(len(Civilizations[-1][2])):
        for mm in range(len(Civilizations[-1][2][nn])):
            rwTNGt.append(np.sort(Civilizations[-1][2][nn][mm])[::-1].tolist())

    if dbg:
        print('(iso)spin configurations:')
        print(Civilizations[-1][0], '\n')
        print('internal widths:')
        print(Civilizations[-1][1], '\n')
        print('relative widths:')
        print(Civilizations[-1][2], '\n')
        print('cfg boundaries:')
        print(cfgbnds, '\n')
        print(bvcfgs)
        print('D0: ', D0, '\n')
        print('D1: ', D1, '\n')

    cfgsNext = []
    iwTNGNext = []
    rwTNGNext = []

    for stru in D1:
        for rwlen in np.unique([len(cff[1]) for cff in D1[stru]]):

            cfgsNext.append(Civilizations[-1][0][D1[stru][0][2]])

            iwTNGNext.append([
                sum(Civilizations[-1][1], [])[bvs[0] - 1] for bvs in D1[stru]
                if len(bvs[1]) == rwlen
            ])
            rwTNGNext.append([
                np.array(rwTNGt[bvs[0] - 1])[(np.array(bvs[1]) -
                                              1).tolist()].tolist()
                for bvs in D1[stru] if len(bvs[1]) == rwlen
            ])

    if dbg:
        print('NEXT (iso)spin configurations:')
        print(cfgsNext, '\n')
        print('NEXT internal widths:')
        print(iwTNGNext, '\n')
        print('NEXT relative widths:')
        print(rwTNGNext, '\n')

    Civilizations.append([cfgsNext, iwTNGNext, rwTNGNext, []])

    if dbg:
        print('fit and stable basis %d/%d populated.' % (nCivi, nRaces))

    Dp = []
    anzPro = 0
    for cfg in range(len(Civilizations[-1][1])):
        for bv in range(len(Civilizations[-1][1][cfg])):
            anzPro += 1
            tnp = np.arange(1, len(Civilizations[-1][2][cfg][bv]) + 1).tolist()
            Dp.append([anzPro, tnp])

    Ais = [
        Civilizations[-1][0], Civilizations[-1][1], Civilizations[-1][2], Dp
    ]

    for nGen in range(nGenMax):

        print('\n Parent generation %d: |Ais| = %d' % (nGen, basisDim(Ais[3])))

        ma = blunt_ev(Ais[0],
                      Ais[1],
                      Ais[2],
                      Ais[3],
                      dia=False,
                      wrkdir='',
                      nzopt=zop,
                      costring=costr,
                      bin_path=BINBDGpath,
                      einzel_file_path=v18uixpath,
                      potNN=potnn,
                      potNNN=potnnn,
                      parall=-1,
                      tnni=10,
                      jay=Jstreu)

        # generate a pair of offspring from a randomly selected couple
        offs_set = []
        for n in range(nbrOff):
            son, daughter = breed_offspring(Ais[1], Ais[2])
            offs_set.append(daughter)
            offs_set.append(son)

        # add <offs_set> to the parent population
        assert len(offs_set) == (2 * nbrOff)
        Ai = []
        for child in range(len(offs_set)):

            ais0 = Ais[0].copy()
            ais1 = Ais[1].copy()
            ais2 = Ais[2].copy()
            ais3 = Ais[3].copy()
            ais0.append(Ais[0][offs_set[child][2]])
            ais1.append([offs_set[child][0]])
            ais2.append([[offs_set[child][1]]])

            aitmp0 = ais0.copy()
            aitmp1 = ais1.copy()
            aitmp2 = ais2.copy()
            aitmp3 = ais3.copy()

            Dp1 = []
            anzPro = 0
            for cfg in range(len(ais1)):
                for bv in range(len(ais1[cfg])):
                    anzPro += 1
                    tnp = np.arange(1, len(ais2[cfg][bv]) + 1).tolist()
                    Dp1.append([anzPro, tnp])

            aitmp3 = Dp1

            Ai.append([aitmp0, aitmp1, aitmp2, aitmp3])

        for m in range(len(Ai)):
            ma = blunt_ev(Ai[m][0],
                          Ai[m][1],
                          Ai[m][2],
                          Ai[m][3],
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

            basQuala_m = ewH[-1]
            basQualb_m = len([bvv for bvv in ewN if bvv < ewMax]) + len(
                [bvv for bvv in ewH if bvv < ewMax])
            basCond_m = np.min(np.abs(ewN)) / np.max(np.abs(ewN))
            if basCond_m > minCond:
                Ai[m].append([basQuala_m, basQualb_m, basCond_m])
                if dbg:
                    print('%d) ' % m, basQuala_m, basQualb_m, basCond_m)
            else:
                if dbg:
                    print('%d) ' % m, basQuala_m, basQualb_m, basCond_m,
                          '  (NOT added!)')

        qualCrit = 2
        criterium = [gs[4][qualCrit] for gs in Ai if len(gs) > 4]
        if qualCrit == 0:
            idx = np.argmin(criterium) if criterium != [] else -1
        elif qualCrit == 2:
            idx = np.argmax(criterium) if criterium != [] else -1

        if dbg:
            print('quality of generation %d:\n' % nGen)
            print('fittest element: %d' % idx)

        Ais = Ai[idx][:4] if idx >= 0 else []
        print(Ai[idx][-1])

    if Ais != []:
        Civilizations[-1][1].append(Ais[1][0])
        Civilizations[-1][2].append(Ais[2][0])
        Civilizations[-1][0].append(Ais[0][0])
        print('a new, %d-dimensional race has evolved!' %
              len(sum(sum(Ais[2], []), [])))
        for Civilization in Civilizations:
            print('\n', Civilization)

    else:
        print('no fruitful new generation emerged to replace the elder...')
