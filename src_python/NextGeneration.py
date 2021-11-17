from bridgeA3 import *
from genetic_width_growth import *

import glob, copy
import shlex
import multiprocessing
from multiprocessing.pool import ThreadPool

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

# Read first generation of parents == (B0)
cfgs = [
    con.split() for con in open(basisPath + 'frags_LIT_J%s_%s.dat' %
                                (Jstreustring, bastypes[0]))
]
origCFGs = copy.deepcopy(cfgs)

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

muta_initial = 0.2
# nRaces := |i|
nRaces = 14
nbrOff = 42
targetDimfac0 = 0.95

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
            np.array(range(1, 1 +
                           len(Civilizations[0][2][cfg][nbvc - 1]))).tolist()
        ]]

for nCivi in range(nRaces):

    #    Dfull = []
    #    anzPro = 0
    #    for cfg in range(len(Civilizations[nCivi][1])):
    #        for bv in range(len(Civilizations[nCivi][1][cfg])):
    #            anzPro += 1
    #            tnp = np.arange(1,
    #                            len(Civilizations[nCivi][2][cfg][bv]) +
    #                            1).tolist()
    #            Dfull.append([anzPro, tnp])
    #    if Civilizations[nCivi][3] == []:
    #        Civilizations[nCivi][3] = Dfull

    cfgbnds = np.add.accumulate([len(iws) for iws in Civilizations[nCivi][1]])
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

    if ewH == []:
        print('parent basis unstable:\n', Civilizations[-1][3])

        exit()

    basCond = np.min(np.abs(ewN)) / np.max(np.abs(ewN))

    print(
        'civ-%d) stability (condition number) the old civilization: ' % nCivi,
        basCond)

    # sift through the candidate population and purge it of 'ideling' individuals
    # which are those without whom the ground state's energy drops by less than some epsilon > 0
    D0 = Civilizations[-1][3]

    go = True
    nPW = 0

    while go:

        nPW += 1
        newpopList = []
        go = False

        D0flat = flatten_basis(D0)

        for bvTrail in ['Raeuber'] + D0flat:

            cpy = D0flat.copy()
            if bvTrail in cpy:
                #print('assessing BV: ', bvTrail)
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
                #print('reference E=', Quala)

            else:
                basQuala = specH[-1]
                basQualb = len([bvv for bvv in specN if bvv < ewMax]) + len(
                    [bvv for bvv in specH if bvv < ewMax])

                bvTdiff = np.abs(basQuala - Quala)
                if bvTdiff < mDiff:
                    print('irrelevant BV found: ', bvTrail)
                    print('dE = ', bvTdiff)
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
            print('removed 1/%d basis vectors' % basisDim(D0))

    D0 = rectify_basis(D0flat)
    n3_inen_bdg(D0, Jstreu, costr, fn='INEN', pari=0, nzop=zop, tni=tnni)

    subprocess.run([BINBDGpath + 'TDR2END_NORMAL.exe'],
                   capture_output=True,
                   text=True)

    matout = np.core.records.fromfile('MATOUTB', formats='f8', offset=4)

    ewN, ewH = NormHamDiag(matout)

    basQuala_i = ewH[-1]
    basQualb_i = len([bvv for bvv in ewN if bvv < ewMax]) + len(
        [bvv for bvv in ewH if bvv < ewMax])
    basCond_i = np.min(np.abs(ewN)) / np.max(np.abs(ewN))

    print('civ-%d) basDim = %d) B(GS) = ' % (nCivi, basisDim(D0)), ewH[-1])

    Civilizations[-1][3] = D0

    parBV = D0[-1][0]

    Ais = copy.deepcopy(Civilizations[-1])

    print('\n Parent generation: |Ais| = %d' % (basisDim(Ais[3])))
    # generate a pair of offspring from a randomly selected couple
    offs_set = []
    for n in range(nbrOff):
        son, daughter = breed_offspring(Ais[1], Ais[2])
        offs_set.append(daughter)
        offs_set.append(son)
    # 1) assemble offspring in rectified set
    offs_set.sort(key=lambda tup: tup[2])
    bvsPerCfg = 18

    for childcfgs in np.unique([chf[2] for chf in offs_set]):
        cset = [ch for ch in offs_set if ch[2] == childcfgs]
        tmp = [
            cset[i * bvsPerCfg:(i + 1) * bvsPerCfg]
            for i in range(int(len(cset) / bvsPerCfg))
        ]
        if (len(cset) % bvsPerCfg) != 0:
            tmp += [cset[-(len(cset) % bvsPerCfg):]]
        iwtmp = [[tt[0] for tt in tmpp] for tmpp in tmp]
        rwtmp = [[[tt[1]] for tt in tmpp] for tmpp in tmp]
        for iws in range(len(iwtmp)):
            Ais[0].append(Civilizations[-1][0][childcfgs])
            Ais[1].append(iwtmp[iws])
            Ais[2].append(rwtmp[iws])

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
    chiBV = len(offs_set)
    childishParaSets = [[
        Ais[3] + [[ch, [1]]], Jstreu, costr, zop, 10, ch, BINBDGpath
    ] for ch in range(parBV + 1, 1 + parBV + chiBV)]
    # 3) rate certain partitions according to their fitness
    print('commencing child vetting...')
    if glob.glob('MATOUTB*') != []:
        os.system('rm -rf MATOUTB*')

    pool = ThreadPool(anzproc)
    jobs = []
    pipe_list = []
    for procnbr in range(len(childishParaSets)):
        recv_end, send_end = multiprocessing.Pipe(False)
        pars = childishParaSets[procnbr]
        p = multiprocessing.Process(target=endmat, args=(pars, send_end))
        jobs.append(p)
        pipe_list.append(recv_end)
        p.start()
    for proc in jobs:
        proc.join()
    child_ladder = [x.recv() for x in pipe_list]
    child_ladder_un = copy.deepcopy(child_ladder)

    # bdg ranking
    child_ladder.sort(key=lambda tup: tup[1])

    # attr ranking
    child_ladder.sort(key=lambda tup: tup[0])
    child_ladder.reverse()

    if child_ladder[0][0] > 0.0:

        chopt = child_ladder[0][3] - parBV - 1

        # 3-1) re-order the offspring list
        parvenue = np.array(offs_set)[chopt]
        Aisopt = Civilizations[-1]
        Aisopt[0].append(Civilizations[-1][0][int(parvenue[2])])
        Aisopt[1].append([parvenue[0]])
        Aisopt[2].append([[parvenue[1]]])
        Aisopt[3].append([parBV + 1, [1]])

        if dbg:
            for off in range(len(child_ladder)):
                if off == chopt:
                    print(child_ladder_un[off], offs_set[off], '  (*)')
                else:
                    print(child_ladder_un[off])

        if idx == -1:
            print(
                'civ-%d) all offspring *and* parents are unfit; INCONSISTENCY!'
                % (nCivi))
            exit()

        stratifiedOptCivilization = condense_Basis(Aisopt,
                                                   origCFGs,
                                                   bvsPERcfg=bvsPerCfg)

        Civilizations.append(stratifiedOptCivilization)

        print(
            'civ-%d) a new, %d-dimensional race has evolved!\n -----------------------------'
            % (nCivi, basisDim(stratifiedOptCivilization[3])))

    else:
        print(
            'civ-%d) spoiled offspring, only, was raised!\n -----------------------------'
            % nCivi)
        muta_initial *= 1.1
        if muta_initial >= 1.0:
            print('... although wild mutations were allowed. Aborting...')
            break

print('After %d generations,\n-----\n' % nCivi, Civilizations[-1],
      '\n-----\nemerged as the dominant culture.')

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
              jay=0.5,
              dia=True)

os.system('rm -rf ./inen_*')
os.system('rm -rf ./endout_*')
os.system('rm -rf ./MATOUTB_*')
os.system('rm -rf ./T*OUT.*')
os.system('rm -rf ./D*OUT.*')

exit()