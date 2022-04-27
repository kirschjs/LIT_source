import operator

from bridgeA3 import *
from parameters_and_constants import *


def span_initial_basis(
    basisType,
    coefstr,
    anzOp=14,
    ini_grid_bounds=[0.1, 9.5, 0.01, 4.5, 0.2, 10.5, 0.02, 5.5],
    ini_dims=[4, 8, 4, 8],
):

    wrkDir = os.getcwd()

    angu = channels[basisType]
    Jstreu = float(basisType.split('^')[0][-3:])
    Jstreustring = '%s' % str(Jstreu)[:3]

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

    # minimal distance allowed for between width parameters
    mindist_int = 0.001

    # lower bound for width parameters '=' IR cutoff (broadest state)
    rWmin = 0.0001

    # orbital-angular-momentum dependent upper bound '=' UV cutoff (narrowest state)
    iLcutoff = [12., 4., 3.]
    rLcutoff = [12., 4., 3.]
    if basisType == boundstatekanal:
        nwint = ini_dims[0]
        nwrel = ini_dims[1]
        rel_scale = 1.
        wi, wf, nw = ini_grid_bounds[0], ini_grid_bounds[1], [
            nwint for n in lfrags
        ]  # initial helion bound state
    else:
        nwint = ini_dims[2]
        nwrel = ini_dims[3]
        rel_scale = 1.
        wi, wf, nw = ini_grid_bounds[4], ini_grid_bounds[5], [
            nwint for n in lfrags
        ]  # final-state continuum

    if nwrel >= rwma:
        print(
            'The set number for relative width parameters per basis vector > max!'
        )
        exit()

    lit_rw_sparse = np.empty(len(sfrags), dtype=list)
    for frg in range(len(lfrags)):
        Lsum = np.sum([int(ie) for ie in lfrags[frg]])
        #  -- internal widths --------------------------------------------------
        offset = 1.
        #if (sfrags[frg][-1] == 'y'):
        if nw[frg] != 1:
            offset += 0.1 * frg / (1 + len(lfrags[frg]))

        wii = wi * offset
        wff = wf * offset

        lit_w_tmp = np.abs(
            np.geomspace(start=wii,
                         stop=wff,
                         num=nw[frg],
                         endpoint=True,
                         dtype=None))

        if nw[frg] != 1:
            lit_w_tmp = np.sort([wd * np.random.random()
                                 for wd in lit_w_tmp])[::-1]

        lit_w[frg] = lit_w_tmp

        lit_w[frg] = [
            ww for ww in sparse(lit_w[frg], mindist=mindist_int)
            if rWmin < ww < iLcutoff[int(
                np.max([float(lfrags[frg][0]),
                        float(lfrags[frg][1])]))]
        ]
        #  -- relative widths --------------------------------------------------

        if basisType == boundstatekanal:
            wir, wfr, nwr = rel_scale * ini_grid_bounds[
                2], rel_scale * ini_grid_bounds[3], nwrel * len(lit_w[frg])
        else:
            wir, wfr, nwr = rel_scale * ini_grid_bounds[
                6], rel_scale * ini_grid_bounds[7], nwrel * len(lit_w[frg])

        offset = 0.1 + 0.2 * np.random.random() if nwr != 1 else 1.0

        wiir = wir
        wffr = wfr

        lit_w_tmp = np.geomspace(start=wiir,
                                 stop=wffr,
                                 num=nwr,
                                 endpoint=True,
                                 dtype=None)

        lit_w_tmp = offset * lit_w_tmp

        if nwr != 1:
            lit_w_tmp = np.sort([wd * np.random.random()
                                 for wd in lit_w_tmp])[::-1]

        lit_rw_tmp = [
            ww for ww in np.abs(
                np.sort(np.array(lit_w_tmp).flatten())[::-1].tolist())
            if rWmin < ww < rLcutoff[int(
                np.max([float(lfrags[frg][0]),
                        float(lfrags[frg][1])]))]
        ]
        if lit_rw_tmp == []:
            lit_rw_tmp = [np.random.random()]

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
    print(
        'seed state with (%d) basis-vector blocks with [orbital][(iso)spin] configurations:'
        % anzBV)
    print(lfrags2, sfrags2, '\n')

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
                    x for x in range(1, 1 + max([len(wid)
                                                 for wid in widr[n]]), 1)
                ]
            ]]
            bv += 1

    path_bas_dims = wrkDir + '/basis_struct/LITbas_dims_%s.dat' % basisType
    with open(path_bas_dims, 'wb') as f:
        np.savetxt(f, [np.size(wid) for wid in widr], fmt='%d')
        f.seek(NEWLINE_SIZE_IN_BYTES, 2)
        f.truncate()
    f.close()
    path_bas_int_rel_pairs = wrkDir + '/basis_struct/LITbas_full_%s.dat' % basisType
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
    path_frag_stru = wrkDir + '/basis_struct/frags_LIT_%s.dat' % basisType
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
    path_intw = wrkDir + '/basis_struct/intw3heLIT_%s.dat' % basisType
    if os.path.exists(path_intw): os.remove(path_intw)
    with open(path_intw, 'wb') as f:
        for ws in widi:
            np.savetxt(f, [ws], fmt='%12.6f', delimiter=' ')
        f.seek(NEWLINE_SIZE_IN_BYTES, 2)
        f.truncate()
    f.close()
    path_relw = wrkDir + '/basis_struct/relw3heLIT_%s.dat' % basisType
    if os.path.exists(path_relw): os.remove(path_relw)
    with open(path_relw, 'wb') as f:
        for wss in widr:
            for ws in wss:
                np.savetxt(f, [ws], fmt='%12.6f', delimiter=' ')
        f.seek(NEWLINE_SIZE_IN_BYTES, 2)
        f.truncate()
    f.close()
    if os.path.isdir(wrkDir + '/eob/') == False:
        subprocess.check_call(['mkdir', '-p', wrkDir + '/eob/'])
    os.chdir(wrkDir + '/eob/')
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
    if os.path.isdir(wrkDir + '/eob-tni/') == False:
        subprocess.check_call(['mkdir', '-p', wrkDir + '/eob-tni/'])
    os.chdir(wrkDir + '/eob-tni/')
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
    if os.path.isdir(wrkDir + '/elu/') == False:
        subprocess.check_call(['mkdir', '-p', wrkDir + '/elu/'])
    os.chdir(wrkDir + '/elu/')
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
                '123',
                '121',
                '122',
                '212',
                '222',
                '221',
                '220',
            ],
            indep=+1)
    os.system(BINBDGpath + 'LUDW_CN.exe')
    if os.path.isdir(wrkDir + '/elu-tni/') == False:
        subprocess.check_call(['mkdir', '-p', wrkDir + '/elu-tni/'])
    os.chdir(wrkDir + '/elu-tni/')
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
                '123',
                '121',
                '122',
                '212',
                '222',
                '221',
                '220',
            ],
            indep=+1)
    os.system(BINBDGpath + 'DRLUD.exe')

    os.chdir(wrkDir)

    n3_inlu(8, fn='INLU', fr=lfrags2, indep=parall)
    os.system(BINBDGpath + 'DRLUD.exe')
    n3_inlu(8, fn='INLUCN', fr=lfrags2, indep=parall)
    os.system(BINBDGpath + 'LUDW_CN.exe')
    n3_inob(sfrags2, 8, fn='INOB', indep=parall)
    os.system(BINBDGpath + 'KOBER.exe')
    n3_inob(sfrags2, 15, fn='INOB', indep=parall)
    os.system(BINBDGpath + 'DROBER.exe')
    he3inquaBS(intwi=widi, relwi=widr, potf='./%s' % nnStr)
    parallel_mod_of_3inqua(lfrags2,
                           sfrags2,
                           infile='INQUA_M',
                           outfile='INQUA_M',
                           einzel_path=wrkDir + '/')
    insam(len(lfrags2))

    anzproc = max(2, min(len(lfrags2), MaxProc))
    #print('Anzahl der Sklaven + 1: %d' % anzproc)
    #exit()

    n3_inen_bdg(sbas, Jstreu, coefstr, fn='INEN', pari=0, nzop=anzOp, tni=tnni)

    if parall == -1:
        wrkVol = du(pathbase)
        while int(wrkVol) > homeQuota:
            print('wrkDir holds %d bytes. Waiting for 60s to shrink.' %
                  int(wrkVol))
            time.sleep(60)
            wrkVol = du(pathbase)
        subprocess.run([
            MPIRUN, '-np',
            '%d' % anzproc, BINBDGpath + 'V18_PAR/mpi_quaf_v7'
        ])
        subprocess.run([BINBDGpath + 'V18_PAR/sammel'])
        subprocess.call('rm -rf DMOUT.*', shell=True)
    else:
        subprocess.run([BINBDGpath + 'QUAFL_M.exe'])
    if tnni == 11:
        he3inquaBS(intwi=widi, relwi=widr, potf='./%s' % nnnStr)
        parallel_mod_of_3inqua(lfrags2,
                               sfrags2,
                               infile='INQUA_M',
                               outfile='INQUA_M',
                               tni=1,
                               einzel_path=wrkDir + '/')
        if parall == -1:
            while int(wrkVol) > homeQuota:
                print('wrkDir holds %d bytes. Waiting for 60s to shrink.' %
                      int(wrkVol))
                time.sleep(60)
                wrkVol = du(pathbase)
            subprocess.run([
                MPIRUN, '-np',
                '%d' % anzproc, BINBDGpath + 'UIX_PAR/mpi_drqua_v7'
            ])
            subprocess.run([BINBDGpath + 'UIX_PAR/SAMMEL-uix'])
            subprocess.call('rm -rf DRDMOUT.*', shell=True)
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

    ew_threshold = 10**(-8)

    matout = np.core.records.fromfile('MATOUTB', formats='f8', offset=4)

    return matout


#    goodEVs = smart_ev(matout, ew_threshold)
#    print('the good, big, and the smallest: ', np.real(goodEVs[:4]), ' ... ',
#          np.real(goodEVs[-3:]))
#    print('#E < 0                         : ',
#          len([dg for dg in np.real(goodEVs) if float(dg) <= 0.0]))