import operator

from bridgeA3 import *
from parameters_and_constants import *

with open(respath + 'dtype.dat', 'w') as outf:
    outf.write(dt)
outf.close()

for bastype in bastypes:
    angu = channels[bastype]

    Jstreu = float(bastype.split('^')[0][-3:])
    Jstreustring = '%s' % str(Jstreu)[:3]

    if os.path.isdir(v18uixpath) == False:
        os.mkdir(v18uixpath)
    if os.path.isdir(litpath3He) == False:
        os.mkdir(litpath3He)

    os.chdir(litpath3He)

    if os.path.isdir(litpath3He + 'basis_struct/') == False:
        os.mkdir(litpath3He + 'basis_struct/')
    os.chdir(litpath3He + 'basis_struct/')

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

    he_iw = he_rw = he_frgs = ob_stru = lu_stru = sbas = []

    inv_scale_i = 0.25
    inv_scale_r = 5.1

    nGrd = 14

    nwint = 14
    nwrel = 4

    grd_type = 'geo'  #'cum'  #'poly'  #

    cal += ['diag']

    mindist_int = 0.001

    iLcutoff = [12., 4., 3.]
    rLcutoff = [12., 4., 3.]

    if bastype == boundstatekanal:

        rel_scale = 1.
        wi, wf, nw = 0.001, 11.5, [nwint
                                   for n in lfrags]  # for lit-state continuum

    else:

        rel_scale = 0.1
        wi, wf, nw = 0.0001, 3.1, [nwint
                                   for n in lfrags]  # for lit-state continuum

    # to include all reference basis states in the augmented basis
    #lit_rw_sparse = np.empty(max(len(sfrags), len(ob_stru)), dtype=list)
    # to use a small subset comprising all internal, but not all relative parameters sets
    lit_rw_sparse = np.empty(len(sfrags), dtype=list)

    for frg in range(len(lfrags)):

        Lsum = np.sum([int(ie) for ie in lfrags[frg]])

        #  -- internal widths --------------------------------------------------
        offset = 1.
        #if (sfrags[frg][-1] == 'y'):
        offset += 0.1 * frg / (1 + len(lfrags[frg]))

        wii = wi * offset
        wff = wf * offset

        if grd_type == 'geo':
            lit_w_tmp = np.abs(
                np.geomspace(
                    #lit_w_tmp = np.abs(np.linspace(
                    start=wii,
                    stop=wff,
                    num=nw[frg],
                    endpoint=True,
                    dtype=None))

        lit_w[frg] = lit_w_tmp
        lit_w[frg] = [
            ww for ww in sparse(lit_w[frg], mindist=mindist_int)
            if ww < iLcutoff[int(
                np.max([float(lfrags[frg][0]),
                        float(lfrags[frg][1])]))]
        ]

        #  -- relative widths --------------------------------------------------
        #geomspace
        wir, wfr, nwr = rel_scale * wi, rel_scale * wf, nwrel * len(lit_w[frg])

        offset = 0.95 + 0.1 * np.random.random()

        wiir = wir
        wffr = wfr

        if grd_type == 'geo':
            lit_w_tmp = np.geomspace(
                #lit_w_tmp = np.linspace(
                start=wiir,
                stop=wffr,
                num=nwr,
                endpoint=True,
                dtype=None)

        lit_w_tmp = offset * lit_w_tmp

        lit_rw_tmp = [
            ww for ww in np.abs(
                np.sort(np.array(lit_w_tmp).flatten())[::-1].tolist())
            if ww < rLcutoff[int(
                np.max([float(lfrags[frg][0]),
                        float(lfrags[frg][1])]))]
        ]

        lit_rw[frg] = []
        for bv in range(len(lit_w[frg])):
            lit_rw[frg].append(lit_rw_tmp[bv::len(lit_w[frg])])

        minL = np.min([len(ws) for ws in lit_rw[frg]])
        lit_rw[frg] = [wd[-minL:] for wd in lit_rw[frg]]

    widi = []
    widr = []

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
        tmp3 = [lit_rw[n][bnds[nn]:bnds[nn + 1]] for nn in range(zer_per_ws)]
        sfrags2 += len(tmp2) * [sfrags[n]]
        lfrags2 += len(tmp2) * [lfrags[n]]

        widi += tmp2
        widr += tmp3

    anzBV = sum([len(zer) for zer in widi])
    print('# Basisvektoren = %d' % anzBV)
    print(lfrags2, sfrags2)

    if ((bastype != boundstatekanal) | (new_helion)):
        sbas = []
        bv = 1
        for n in range(len(lfrags2)):
            bvv = 0
            for m in range(len(widi[n])):
                bvv += 1
                #sbas += [[bv, [(bvv) % (1 + len(widr[n][m]))]]]
                sbas += [[
                    bv,
                    [
                        x
                        for x in range(1, 1 +
                                       max([len(wid) for wid in widr[n]]), 1)
                    ]
                ]]
                bv += 1

    else:
        sfrags2 = ob_stru
        lfrags2 = lu_stru
        widi = he_iw
        widr = he_rw
    path_bas_dims = litpath3He + 'basis_struct/LITbas_dims_J%s_%s.dat' % (
        Jstreustring, bastype)
    with open(path_bas_dims, 'wb') as f:
        np.savetxt(f, [np.size(wid) for wid in widr], fmt='%d')
        f.seek(NEWLINE_SIZE_IN_BYTES, 2)
        f.truncate()
    f.close()

    path_bas_int_rel_pairs = litpath3He + 'basis_struct/LITbas_full_J%s_%s.dat' % (
        Jstreustring, bastype)
    if os.path.exists(path_bas_int_rel_pairs):
        os.remove(path_bas_int_rel_pairs)

    with open(path_bas_int_rel_pairs, 'w') as oof:
        #np.savetxt(f, [[jj[0], kk] for jj in sbas for kk in jj[1]], fmt='%d')
        #f.seek(NEWLINE_SIZE_IN_BYTES, 2)
        #f.truncate()
        so = ''
        for bv in sbas:
            so += '%4s' % str(bv[0])
            for rww in bv[1]:
                so += '%4s' % str(rww)
            so += '\n'
        oof.write(so)

    oof.close()

    path_frag_stru = litpath3He + 'basis_struct/frags_LIT_J%s_%s.dat' % (
        Jstreustring, bastype)
    if os.path.exists(path_frag_stru): os.remove(path_frag_stru)
    with open(path_frag_stru, 'wb') as f:
        np.savetxt(f,
                   np.column_stack([sfrags2, lfrags2]),
                   fmt='%s',
                   delimiter=' ',
                   newline=os.linesep)
        f.seek(NEWLINE_SIZE_IN_BYTES, 2)
        f.truncate()
    f.close()

    path_intw = litpath3He + 'basis_struct/intw3heLIT_J%s_%s.dat' % (
        Jstreustring, bastype)
    if os.path.exists(path_intw): os.remove(path_intw)
    with open(path_intw, 'wb') as f:
        for ws in widi:
            np.savetxt(f, [ws], fmt='%12.6f', delimiter=' ')
        f.seek(NEWLINE_SIZE_IN_BYTES, 2)
        f.truncate()
    f.close()

    path_relw = litpath3He + 'basis_struct/relw3heLIT_J%s_%s.dat' % (
        Jstreustring, bastype)
    if os.path.exists(path_relw): os.remove(path_relw)
    with open(path_relw, 'wb') as f:
        for wss in widr:
            for ws in wss:
                np.savetxt(f, [ws], fmt='%12.6f', delimiter=' ')
        f.seek(NEWLINE_SIZE_IN_BYTES, 2)
        f.truncate()
    f.close()

    if (('einzel' in cal) & (parall == -1)):
        if os.path.isdir(v18uixpath + 'eob/') == False:
            os.mkdir(v18uixpath + 'eob/')
        os.chdir(v18uixpath + 'eob/')
        n3_inob([
            'he_no1',
            'he_no1y',
            'he_no2',
            'he_no2y',
            'he_no3',
            'he_no3y',
            'he_no5',
            'he_no5y',
            'he_no6',
            'he_no6y',
        ],
                8,
                fn='INOB',
                indep=+1)
        os.system(BINBDGpath + 'KOBER.exe')
        if os.path.isdir(v18uixpath + 'eob-tni/') == False:
            os.mkdir(v18uixpath + 'eob-tni/')
        os.chdir(v18uixpath + 'eob-tni/')
        n3_inob([
            'he_no1',
            'he_no1y',
            'he_no2',
            'he_no2y',
            'he_no3',
            'he_no3y',
            'he_no5',
            'he_no5y',
            'he_no6',
            'he_no6y',
        ],
                15,
                fn='INOB',
                indep=+1)
        os.system(BINBDGpath + 'DROBER.exe')
        if os.path.isdir(v18uixpath + 'elu/') == False:
            os.mkdir(v18uixpath + 'elu/')
        os.chdir(v18uixpath + 'elu/')
        n3_inlu(8,
                fn='INLUCN',
                fr=[
                    '000',
                    '202',
                    '022',
                    '110',
                    '101',
                    '011',
                    '111',
                    '112',
                    '211',
                    '212',
                    '213',
                    '121',
                    '122',
                    '212',
                    '222',
                    '221',
                    '220',
                ],
                indep=+1)
        os.system(BINBDGpath + 'LUDW_CN.exe')
        if os.path.isdir(v18uixpath + 'elu-tni/') == False:
            os.mkdir(v18uixpath + 'elu-tni/')
        os.chdir(v18uixpath + 'elu-tni/')
        n3_inlu(8,
                fn='INLU',
                fr=[
                    '000',
                    '202',
                    '022',
                    '110',
                    '101',
                    '011',
                    '111',
                    '112',
                    '211',
                    '212',
                    '213',
                    '121',
                    '122',
                    '212',
                    '222',
                    '221',
                    '220',
                ],
                indep=+1)
        os.system(BINBDGpath + 'DRLUD.exe')

    os.chdir(v18uixpath)
    print('Calculating in %s' % v18uixpath)
    print(lfrags2)
    n3_inlu(8, fn='INLU', fr=lfrags2, indep=parall)
    os.system(BINBDGpath + 'DRLUD.exe')
    n3_inlu(8, fn='INLUCN', fr=lfrags2, indep=parall)
    os.system(BINBDGpath + 'LUDW_CN.exe')
    n3_inob(sfrags2, 8, fn='INOB', indep=parall)
    os.system(BINBDGpath + 'KOBER.exe')
    n3_inob(sfrags2, 15, fn='INOB', indep=parall)
    os.system(BINBDGpath + 'DROBER.exe')

    he3inquaBS(intwi=widi, relwi=widr, potf=potnn)

    parallel_mod_of_3inqua(lfrags2,
                           sfrags2,
                           infile='INQUA_M',
                           outfile='INQUA_M',
                           einzel_path=v18uixpath)

    insam(len(lfrags2))

    n3_inen_bdg(sbas, Jstreu, costr, fn='INEN', pari=0, nzop=zop, tni=tnni)

    if bastype == boundstatekanal:
        n3_inlu(8, fn=helionpath + 'INLU_ref', fr=lfrags2, indep=-1)
        n3_inlu(8, fn=helionpath + 'INLUCN_ref', fr=lfrags2, indep=-1)
        n3_inob(sfrags2, 8, fn=helionpath + 'INOB_ref', indep=-1)
        n3_inob(sfrags2, 15, fn=helionpath + 'DRINOB_ref', indep=-1)
        os.system('cp INQUA_M ' + helionpath + 'INQUA_V18_ref')
        os.system('cp INEN ' + helionpath + 'INEN_ref')
        os.system('cp INSAM ' + helionpath)

    if parall == -1:
        subprocess.run([
            'mpirun', '-np',
            '%d' % anzproc, BINBDGpath + 'V18_PAR/mpi_quaf_v7'
        ])
        subprocess.run([BINBDGpath + 'V18_PAR/sammel'])
    else:
        subprocess.run([BINBDGpath + 'QUAFL_M.exe'])

    if tnni == 11:
        he3inquaBS(intwi=widi, relwi=widr, potf=potnnn)

        parallel_mod_of_3inqua(lfrags2,
                               sfrags2,
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
            subprocess.run([BINBDGpath + 'TDR2END_NORMAL.exe'])
            subprocess.call('cp OUTPUT out_normal', shell=True)
        else:
            subprocess.run([BINBDGpath + 'DRQUA_AK_M.exe'])
            subprocess.run([BINBDGpath + 'DR2END_AK.exe'])
    elif tnni == 10:
        if parall == -1:
            subprocess.run([BINBDGpath + 'TDR2END_NORMAL.exe'])
            subprocess.call('cp OUTPUT out_normal', shell=True)
        else:
            subprocess.run([BINBDGpath + 'DR2END_NORMAL.exe'])

    suche_fehler()

    if bastype == boundstatekanal:
        os.system('cp INQUA_M ' + helionpath + 'INQUA_UIX_ref')

    ew_threshold = 10**(-8)

    subprocess.call('cp MATOUTB %smat_%s' % (respath, bastype), shell=True)
    #matout = [line for line in open('MATOUTB')]
    matout = np.core.records.fromfile('MATOUTB', formats='f8', offset=4)

    if 'diag' in cal:
        goodEVs = smart_ev(matout, ew_threshold)

        print('the good, big, and the smallest: ', np.real(goodEVs[:4]),
              ' ... ', np.real(goodEVs[-3:]))
        print('#E < 0                         : ',
              len([dg for dg in np.real(goodEVs) if float(dg) <= 0.0]))