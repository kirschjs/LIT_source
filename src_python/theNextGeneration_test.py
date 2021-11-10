from bridgeA3 import *
from genetic_width_growth import *

os.chdir(litpath3He)

basisPath = litpath3He + 'basis_struct/'
bastypes = [boundstatekanal]  # + streukas

for bastype in bastypes:

    angu = channels[bastype]

    Jstreu = float(bastype.split('^')[0][-3:])
    Jstreustring = '%s' % str(Jstreu)[:3]

    # Read a (hopefully stable) parent basis == (D')
    Dprime = [[
        int(bv.strip().split()[0]),
        np.array(bv.strip().split()[1:]).astype(int).tolist()
    ] for bv in open(basisPath + 'LITbas_full_J%s_%s.dat' %
                     (Jstreustring, bastype))][:2]

    nGen = 0

    # Read first generation of parents == (B0)
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

    nGenMax = 10
    while nGen < nGenMax:

        # for each (iso)spin configuration, select pairs to procreate
        # and breed offsping generation == (C0)
        iTNG = []
        rTNG = []
        for cfg in range(len(intwLIT)):

            TNGt = []

            inxs = np.arange(len(intwLIT[cfg]))
            np.random.shuffle(inxs)
            parentpairs = np.reshape(
                inxs, (-1, 2)) if len(inxs) % 2 == 0 else np.reshape(
                    np.append(inxs, np.random.choice(inxs[:-1])), (-1, 2))

            for coupl in parentpairs:
                father = intwLIT[cfg][coupl[0]]
                mother = intwLIT[cfg][coupl[1]]
                son, daughter = intertwining(father, mother)
                TNGt.append([son, daughter])
            iTNG.append(
                np.sort(
                    np.random.choice(
                        np.ravel(TNGt), size=len(intwLIT[cfg]),
                        replace=False))[::-1].astype(float).tolist())

            rTNG.append([])
            for relwS in range(len(intwLIT[cfg])):
                inxs = np.arange(len(rws[cfg][relwS]))
                np.random.shuffle(inxs)
                parentpairs = np.reshape(
                    inxs, (-1, 2)) if len(inxs) % 2 == 0 else np.reshape(
                        np.append(inxs, np.random.choice(inxs[:-1])), (-1, 2))
                TNGt = []
                for coupl in parentpairs:
                    father = rws[cfg][relwS][coupl[0]]
                    mother = rws[cfg][relwS][coupl[1]]
                    son, daughter = intertwining(father, mother)
                    TNGt.append([son, daughter])
                rTNG[-1].append(
                    np.sort(np.ravel(TNGt))[::-1].astype(float).tolist())

        iwTNG = intwLIT + iTNG
        rwTNG = rws + rTNG

        tngBasis = []
        anzOld = len(Dprime)

        anzNew = 0
        dimDfull = 0
        for cfg in range(len(iTNG)):
            for bv in range(len(iTNG[cfg])):
                dimDfull += len(rwTNG[cfg][bv])
                anzNew += 1
                tnp = np.arange(1, len(rTNG[cfg][bv]) + 1).tolist()
                tngBasis.append([anzOld + anzNew, tnp])

        cfgs = [
            con.split() for con in open(basisPath + 'frags_LIT_J%s_%s.dat' %
                                        (Jstreustring, bastype))
        ]

        lfragTNG = 2 * [cc[1] for cc in cfgs]
        sfragTNG = 2 * [cc[0] for cc in cfgs]
        cfgs = 2 * cfgs
        cfgbnds = np.add.accumulate([len(iws) for iws in iwTNG])
        cfgbnds = np.insert(cfgbnds, 0, 0)

        # calculate Norm/Hamiltonian system in parent+offspring basis == (Dfull)
        Dfull = Dprime + tngBasis

        # test every c(andidate) vector in Dfull which is not in the initial D0 = Dprime
        D0 = Dprime.copy()

        # update initial parent population:
        # relwLIT, rws, intwLIT
        rectbas = []
        for bv in D0:
            if bv[0] in [b[0] for b in rectbas]:
                rectbas[[b[0]
                         for b in rectbas].index(bv[0])][1].append(bv[1][0])
            else:
                rectbas.append(bv)
        idx = np.array([b[0] for b in rectbas]).argsort()[::-1]
        rectbas = [[bb[0], np.sort(bb[1]).tolist()] for bb in rectbas]

        tmpBasis = rectbas
        print('relevant basis vectors:')
        print(tmpBasis, '\n')
        print('(iso)spin configurations:')
        print(cfgs, '\n')
        print('internal widths:')
        print(iwTNG, '\n')
        print('cfg boundaries:')
        print(cfgbnds, '\n')
        print('relative widths:')
        print(rwTNG, '\n')

        bvcfgs = [[
            ncfg for ncfg in range(len(cfgbnds) - 1)
            if ((bv[0] > cfgbnds[ncfg]) & (bv[0] <= cfgbnds[ncfg + 1]))
        ] for bv in tmpBasis]

        print(bvcfgs)

        exit()
