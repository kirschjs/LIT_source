from bridgeA2 import *
from two_particle_functions import *
from parameters_and_constants import *
import operator

bastypes = streukas
bastypes = [boundstatekanal]
bastypes = [boundstatekanal] + streukas

for bastype in bastypes:
    angu = channels[bastype]
    Jstreu = float(bastype.split('^')[0][-1])
    Jstreustring = '%s' % str(Jstreu)[:3]

    if os.path.isdir(v18uixpath) == False:
        os.mkdir(v18uixpath)
    if os.path.isdir(litpathD) == False:
        os.mkdir(litpathD)

    os.chdir(litpathD)

    if os.path.isdir(litpathD + 'basis_struct/') == False:
        os.mkdir(litpathD + 'basis_struct/')
    os.chdir(litpathD + 'basis_struct/')

    lit_w = {}
    lit_rw = {}
    lfrags = []
    sfrags = []
    lfrags2 = []
    sfrags2 = []
    for lcfg in range(len(angu)):
        sfrags = sfrags + angu[lcfg][1]
        for scfg in angu[lcfg][1]:
            lfrags = lfrags + [angu[lcfg][0]]

    wli = 'lin'

    he_iw = he_rw = he_frgs = ob_stru = lu_stru = sbas = []

    ew_threshold = 10**(-10)

    bvma = 6

    if bastype == boundstatekanal:

        ngeosets = 1
        nw = 20
        wi, wf = 0.0002, 8.5

    else:
        ngeosets = 1
        nw = 36
        wi, wf = 0.0002, 5.15

    # to include all reference basis states in the augmented basis
    #lit_rw_sparse = np.empty(max(len(sfrags), len(ob_stru)), dtype=list)
    # to use a small subset comprising all internal, but not all relative parameters sets
    lit_rw_sparse = np.empty(len(sfrags), dtype=list)

    #  -- internal widths --------------------------------------------------
    lit_w_tmp = [
        np.abs(
            np.geomspace(
                start=wi,  # * (1 + 1 * n / (ngeosets)),
                stop=wf,  # * (1 + n / ngeosets),
                num=nw,
                endpoint=True,
                dtype=None)) for n in range(ngeosets)
    ]
    lit_w_tmp = [np.array([vv for vv in lit_w_tmp[0] if vv < 100])]

    lit_w[0] = np.reshape(lit_w_tmp, (1, -1))[0]

    #    lit_w_tmp_add = np.abs(
    #        np.geomspace(
    #            start=wf, stop=2 * wf, num=nwadd, endpoint=True, dtype=None))
    #    lit_w[0] = [
    #        float(x)
    #        for x in np.sort(np.concatenate((lit_w_tmp, lit_w_tmp_add), axis=0))
    #    ] if nwadd != 0 else lit_w_tmp

    #        lit_w[frg] = sparse(lit_w[frg], mindist=mindist_int)
    for frg in range(1, len(lfrags)):
        lit_w_tmp = [
            frg * (lit_w[0][n + 1] - lit_w[0][n]) / len(lfrags)
            for n in range(ngeosets * nw - 1)
        ] + [lit_w[0][-1]]
        lit_w[frg] = lit_w_tmp + lit_w[0]

    widi = []
    for n in range(len(lit_w)):

        tmp = np.sort(lit_w[n])[::-1]
        #tmp = sparse(tmp, mindist_int)
        zer_per_ws = int(np.ceil(len(tmp) / bvma))

        bins = [0 for nmmm in range(zer_per_ws + 1)]
        bins[0] = 0
        for mn in range(len(tmp)):
            bins[1 + mn % zer_per_ws] += 1
        bnds = np.cumsum(bins)

        tmp2 = [list(tmp[bnds[nn]:bnds[nn + 1]]) for nn in range(zer_per_ws)]
        sfrags2 += len(tmp2) * [sfrags[n]]
        lfrags2 += len(tmp2) * [lfrags[n]]

        widi += tmp2

    anzBV = sum([len(zer) for zer in widi])
    print('# Basisvektoren = %d' % anzBV)
    print(lfrags2)

    if ((bastype != boundstatekanal) | (new_deuteron)):
        sbas = []
        bv = 1
        for n in range(len(lfrags2)):
            sbas += [[bv, range(1, 1 + len(widi[n]))]]
            bv += 1
    else:
        sfrags2 = ob_stru
        lfrags2 = lu_stru
        widi = he_iw

    path_bas_int_rel_pairs = litpathD + 'basis_struct/LITbas_full_J%s_%s.dat' % (
        Jstreustring, bastype)
    if os.path.exists(path_bas_int_rel_pairs):
        os.remove(path_bas_int_rel_pairs)

    with open(path_bas_int_rel_pairs, 'wb') as f:
        np.savetxt(f, [[jj[0], kk] for jj in sbas for kk in jj[1]], fmt='%d')
    f.close()

    path_frag_stru = litpathD + 'basis_struct/frags_LIT_J%s_%s.dat' % (
        Jstreustring, bastype)
    if os.path.exists(path_frag_stru): os.remove(path_frag_stru)
    with open(path_frag_stru, 'wb') as f:
        np.savetxt(
            f,
            np.column_stack([sfrags2, lfrags2]),
            fmt='%s',
            delimiter=' ',
            newline=os.linesep)
    f.close()

    path_intw = litpathD + 'basis_struct/intwDLIT_J%s_%s.dat' % (Jstreustring,
                                                                 bastype)
    if os.path.exists(path_intw): os.remove(path_intw)
    with open(path_intw, 'wb') as f:
        for ws in widi:
            np.savetxt(f, [ws], fmt='%12.6f', delimiter=' ; ')

    f.close()

    os.chdir(v18uixpath)
    print('Calculating in %s' % v18uixpath)
    n2_inlu(8, fn='INLUCN', fr=lfrags2, npol=npoli)
    os.system(BINBDGpath + 'LUDW_CN.exe')
    n2_inob(sfrags2, 8, fn='INOB')
    os.system(BINBDGpath + 'KOBER.exe')

    DinquaBS(intwi=widi, potf=potnn, npol=npoli)
    #n2_inen_bdg(sbas, Jstreu, costr, fn='INEN', pari=0)
    n2_inen_pol(sbas, Jstreu, costr, fn='INEN', pari=0, npol=npoli)

    if bastype == boundstatekanal:
        n2_inlu(
            8,
            fn=deuteronpath + 'INLUCN_ref',
            fr=lfrags2,
            indep=-1,
            npol=npoli)
        n2_inob(sfrags2, 8, fn=deuteronpath + 'INOB_ref', indep=-1)
        os.system('cp INQUA_M ' + deuteronpath + 'INQUA_V18_ref')
        os.system('cp INEN ' + deuteronpath + 'INEN_ref')

    subprocess.run([BINBDGpath + 'QUAFL_M.exe'])

    subprocess.run([BINBDGpath + 'DR2END_NORMAL.exe'])

    suche_fehler()

    subprocess.call('cp MATOUTB %smat_%s' % (respath, bastype), shell=True)
    #matout = [line for line in open('MATOUTB')]
    matout = np.core.records.fromfile('MATOUTB', formats='f8', offset=4)

    goodEVs = smart_ev(matout, ew_threshold)

    print('the good guys: ', np.real(goodEVs[:4]), ' ... ',
          np.real(goodEVs[-4:]))