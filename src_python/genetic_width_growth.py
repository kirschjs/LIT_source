import struct
import numpy as np


def breed_offspring(iws, rws, parentBVs=[]):

    anzBV = len(sum(sum(rws, []), []))

    cand = np.random.choice(range(1, 1 + anzBV),
                            (1,
                             2)).tolist()[0] if parentBVs == [] else parentBVs

    cf = []
    for bv in cand:
        ws = 0
        for cfg in range(len(iws)):
            for nbv in range(len(iws[cfg])):
                for rw in range(len(rws[cfg][nbv])):
                    ws += 1
                    if ws == bv:
                        cf.append(cfg)

    mw = retr_widths(cand[0], iws, rws)
    fw = retr_widths(cand[1], iws, rws)

    try:
        chint = intertwining(mw[0], fw[0], mutation_rate=0.1)
        chrel = intertwining(mw[1], fw[1], mutation_rate=0.1)
    except:
        print(cand[0])
        print(cand[1])
        exit()

    c0 = [chint[0], chrel[0], cf[0]]
    c1 = [chint[1], chrel[1], cf[1]]

    return [c0, c1]


def retr_widths(bvNbr, iws, rws):

    assert bvNbr <= len(sum(sum(rws, []), []))

    ws = 0
    for cfg in range(len(iws)):
        for nbv in range(len(iws[cfg])):
            for rw in range(len(rws[cfg][nbv])):
                ws += 1
                if ws == bvNbr:
                    return [iws[cfg][nbv], rws[cfg][nbv][rw]]


def flatten_basis(basis):

    fb = []
    for bv in basis:
        for rw in bv[1]:
            fb.append([bv[0], [rw]])
    return fb


def rectify_basis(basis):

    rectbas = []
    for bv in basis:
        if bv[0] in [b[0] for b in rectbas]:
            rectbas[[b[0] for b in rectbas].index(bv[0])][1].append(bv[1][0])
        else:
            rectbas.append(bv)
    idx = np.array([b[0] for b in rectbas]).argsort()[::-1]
    rectbas = [[bb[0], np.sort(bb[1]).tolist()] for bb in rectbas]

    return rectbas


def basisDim(bas=[]):
    dim = 0
    for bv in bas:
        for rw in bv[1]:
            dim += 1
    return dim


def select_random_basis(basis_set, target_dim):

    assert target_dim < basisDim(basis_set)

    basv = []
    for bv in basis_set:
        for rw in bv[1]:
            basv.append([bv[0], [rw]])

    # split in 2 in order to sample first (parents) and second (children) part of the basis
    basv1 = basv[:int(len(basv) / 2)]
    basv2 = basv[int(len(basv) / 2):]

    np.random.shuffle(basv1)
    np.random.shuffle(basv2)
    tmp = basv1[:int(target_dim / 2)] + basv2[:int(target_dim / 2)]

    basv = (np.array(tmp)[np.array([ve[0] for ve in tmp]).argsort()]).tolist()

    rndBas = rectify_basis(basv)

    return rndBas


def float_to_bin(num):
    return format(struct.unpack('!I', struct.pack('!f', num))[0], '032b')


def bin_to_float(binary):
    return struct.unpack('!f', struct.pack('!I', int(binary, 2)))[0]


def intertwining(p1, p2, mutation_rate=0.0, wMin=0.0001, wMax=120., dbg=False):

    Bp1 = float_to_bin(p1)
    Bp2 = float_to_bin(p2)

    assert len(Bp1) == len(Bp2)
    assert mutation_rate < 1

    mutationMask = np.random.choice(2,
                                    p=[1 - mutation_rate, mutation_rate],
                                    size=len(Bp1))

    pivot = np.random.randint(0, len(Bp1))

    Bchild1 = Bp1[:pivot] + Bp2[pivot:]
    Bchild2 = Bp2[:pivot] + Bp1[pivot:]

    Bchild2mutated = ''.join(
        (mutationMask | np.array(list(Bchild2)).astype(int)).astype(str))
    Bchild1mutated = ''.join(
        (mutationMask | np.array(list(Bchild1)).astype(int)).astype(str))

    Fc1 = np.abs(bin_to_float(Bchild1mutated))
    Fc2 = np.abs(bin_to_float(Bchild2mutated))

    if (np.isnan(Fc1) | np.isnan(Fc2) | (Fc1 < wMin) | (Fc1 > wMax) |
        (Fc2 < wMin) | (Fc2 > wMax)):
        Fc1 = Fc2 = 42.1

    if (dbg | np.isnan(Fc1) | np.isnan(Fc2)):
        print('parents (binary)        :%12.4f%12.4f' % (p1, p2))
        print('parents (decimal)       :', Bp1, ';;', Bp2)
        print('children (binary)       :', Bchild1, ';;', Bchild2)
        print('children (decimal)      :%12.4f%12.4f' % (Fc1, Fc2))

    return Fc1, Fc2


#intertwining(0.1, 0.01, 0.1, True)