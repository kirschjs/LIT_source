from bridge import *
from three_particle_functions import *
from triton_width_gen import *
from parameters_and_constants import *
import operator
from BasisVisualization import visbas
from plot_basis_3bdy import plotbasis3bdy
from plot_spectrum import plotHspec
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs
from sklearn.neighbors import BallTree
from diversipy import *

bastypes = streukas
bastypes = [boundstatekanal]
bastypes = streukas + [boundstatekanal]

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

    wli = 'lin'

    bvma = 8
    rvma = 20

    if bastype == boundstatekanal:

        inv_scale_i = 2.0
        inv_scale_r = 15.0
        nwinv = 4

        rel_scale = 0.01

        nwint = 8
        nwrel = 8
        siffux = '_uix'
        he_iw, he_rw, he_frgs = retrieve_he3_widths(
            pathbase + '/systems/he3/dim70/INQUA_N%s' % siffux)
        ob_stru = read_inob(pathbase + '/systems/he3/dim70/INOB%s' % siffux)
        lu_stru = read_inlu(pathbase + '/systems/he3/dim70/INLU%s' % siffux)

        sbas = get_bsv_rw_idx(
            inen=pathbase + '/systems/he3/dim70/INEN%s' % siffux)
        if len(lu_stru) != len(ob_stru):
            print('reference he3 input inconsistent.')
            exit()

        he_iw = he_rw = he_frgs = ob_stru = lu_stru = sbas = []

        wi, wf, nw = 0.012, 1.02, [
            np.max([2, int(nwint - np.max([1, int(n[0]) + int(n[1])]))])
            for n in lfrags
        ]  # for lit-state continuum
        mindist_int = 0.01
        #mindist_rel = 0.04
        #print('#iw / Z = ', nw, '\n#rw     =  ', nwrel, lfrags)
        #print('3-helium reference structure:\n', he_frgs, '\n', ob_stru, '\n',
        #      lu_stru)

    else:

        inv_scale_i = 2.0
        inv_scale_r = 15.0
        nwinv = 4

        rel_scale = 0.01

        nwint = 8
        nwrel = 8

        he_iw = he_rw = he_frgs = ob_stru = lu_stru = sbas = []

        wi, wf, nw = 0.01, 0.5, [
            np.max([2, int(nwint - np.max([1, int(n[0]) + int(n[1])]))])
            for n in lfrags
        ]  # for lit-state continuum
        mindist_int = 0.01
        #mindist_rel = 0.1

    # to include all reference basis states in the augmented basis
    #lit_rw_sparse = np.empty(max(len(sfrags), len(ob_stru)), dtype=list)
    # to use a small subset comprising all internal, but not all relative parameters sets
    lit_rw_sparse = np.empty(len(sfrags), dtype=list)

    for frg in range(len(lfrags)):

        if wli == 'lin':
            #  -- internal widths --------------------------------------------------
            offset = frg / len(lfrags) * (wf - wi) / nw[frg]

            lit_w_tmp = np.abs(
                np.geomspace(
                    start=wi + offset,
                    stop=wf + offset,
                    num=nw[frg],
                    endpoint=True,
                    dtype=None))

            lit_w_tmp_add = np.abs(
                np.geomspace(
                    start=wf + offset,
                    stop=inv_scale_i * (2 * wf - wi) + offset,
                    num=nwinv,
                    endpoint=True,
                    dtype=None))

            lit_w[frg] = [
                float(x) for x in np.sort(
                    np.concatenate((lit_w_tmp, lit_w_tmp_add), axis=0))
            ] if nwinv != 0 else lit_w_tmp

            lit_w[frg] = sparse(lit_w[frg], mindist=mindist_int)
            #  -- relative widths --------------------------------------------------
            #linspace
            wir, wfr, nwr = rel_scale * wi, rel_scale * wf, nwrel
            offset = (1 + frg) / len(lfrags) * (wfr - wir) / nwr

            lit_w_tmp = np.geomspace(
                start=wir + offset,
                stop=wfr + offset,
                num=nwr,
                endpoint=True,
                dtype=None)

            lit_w_add = np.geomspace(
                start=wfr + offset,
                stop=inv_scale_r * (2 * wfr - wir) + offset,
                num=nwinv,
                endpoint=True,
                dtype=None)

            lit_w_tmp = [
                float(x) for x in np.sort(
                    np.concatenate((lit_w_tmp, lit_w_add), axis=0))
            ] if nwinv != 0 else lit_w_tmp

            lit_rw[frg] = np.abs(
                np.sort(np.array(lit_w_tmp).flatten())[::-1].tolist())

            for nf in range(len(ob_stru)):
                if ((ob_stru[nf][:-1] == sfrags[frg][:-1]) &
                    (lu_stru[nf] == lfrags[frg])):
                    lit_w[frg] = np.concatenate(
                        (lit_w[frg], np.array(he_iw[nf])), axis=0)
                    #           lit_rw_sparse[frg] = np.array(he_rw[nf])
                    lit_rw[frg] = np.concatenate(
                        (lit_rw[frg], np.array(he_rw[nf])), axis=0)
                    lit_rw[frg] = sparse(lit_rw[frg], mindist=mindist_int)
                    #lit_rw[frg] = sparse2(
                    #    lit_w[frg], lit_rw[frg], mindist=mindist_rel)
                lit_rw[frg] = sparse(lit_rw[frg], mindist=mindist_int)
                #lit_rw[frg] *= np.sort(np.random.uniform(0.4, 1.1, 1))[::-1]

    #print(np.sort(lit_w[0])[::-1])
    #print(np.sort(lit_rw[0])[::-1])
    #exit()

    widi = []
    widr = []

    for n in range(len(lit_w)):

        tmp = np.sort(lit_w[n])[::-1]
        tmp = sparse(tmp, mindist_int)
        zer_per_ws = int(np.ceil(len(tmp) / bvma))

        bins = [0 for nmmm in range(zer_per_ws + 1)]
        bins[0] = 0
        for mn in range(len(tmp)):
            bins[1 + mn % zer_per_ws] += 1
        bnds = np.cumsum(bins)

        tmp2 = [list(tmp[bnds[nn]:bnds[nn + 1]]) for nn in range(zer_per_ws)]

        tmpr = np.sort(lit_rw[n])[::-1]
        tmpr = sparse(tmpr, mindist_int)
        r_per_ws = int(np.ceil(len(tmpr) / rvma))

        binsr = [0 for nmmm in range(r_per_ws + 1)]
        binsr[0] = 0
        for mn in range(len(tmpr)):
            binsr[1 + mn % r_per_ws] += 1
        bndsr = np.cumsum(binsr)

        tmp2r = [list(tmpr[bndsr[nn]:bndsr[nn + 1]]) for nn in range(r_per_ws)]

        tmp2i = []
        for ti in tmp2:
            tmp2i += [ti] * (len(bndsr) - 1)

        tmp2r = tmp2r * (len(bnds) - 1)

        sfrags2 += len(tmp2r) * [sfrags[n]]
        lfrags2 += len(tmp2r) * [lfrags[n]]

        widi += tmp2i
        widr += tmp2r

    anzBV = sum([len(zer) for zer in widi])
    print('# Basisvektoren = %d' % anzBV)

    if ((bastype != boundstatekanal) | (new_helion)):
        sbas = []
        bv = 1
        for n in range(len(lfrags2)):

            rwind = nearest(widi[n], widr[n])

            for m in range(len(widi[n])):
                if 43 == 43:
                    #sbas += [[
                    #    bv,
                    #    [x * ((x + m % 2) % 2) for x in range(len(widr[n]))]
                    #]]
                    sbas += [[bv, [x for x in range(1, 1 + len(widr[n]))]]]
                else:
                    sbas += [[bv, rwind[1][m][::-1][:-7]]]
                bv += 1
    else:
        sfrags2 = ob_stru
        lfrags2 = lu_stru
        widi = he_iw
        widr = he_rw

    path_bas_int_rel_pairs = litpath3He + 'basis_struct/LITbas_full_J%s_%s.dat' % (
        Jstreustring, bastype)
    if os.path.exists(path_bas_int_rel_pairs):
        os.remove(path_bas_int_rel_pairs)

    with open(path_bas_int_rel_pairs, 'wb') as f:
        np.savetxt(f, [[jj[0], kk] for jj in sbas for kk in jj[1]], fmt='%d')
    f.close()

    path_frag_stru = litpath3He + 'basis_struct/frags_LIT_J%s_%s.dat' % (
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

    path_intw = litpath3He + 'basis_struct/intw3heLIT_J%s_%s.dat' % (
        Jstreustring, bastype)
    if os.path.exists(path_intw): os.remove(path_intw)
    with open(path_intw, 'wb') as f:
        for ws in widi:
            np.savetxt(f, [ws], fmt='%12.6f', delimiter=' ; ')

    f.close()

    path_relw = litpath3He + 'basis_struct/relw3heLIT_J%s_%s.dat' % (
        Jstreustring, bastype)
    if os.path.exists(path_relw): os.remove(path_relw)
    with open(path_relw, 'wb') as f:
        for ws in widr:
            np.savetxt(f, [ws], fmt='%12.6f', delimiter=' ; ')
    f.close()

    if 'einzel' in cal:
        if os.path.isdir(v18uixpath + 'eob/') == False:
            os.mkdir(v18uixpath + 'eob/')
        os.chdir(v18uixpath + 'eob/')
        n3_inob(
            [
                'he_no1', 'he_no1i', 'he_no2', 'he_no2i', 'he_no6', 'he_no6i',
                'he_no3', 'he_no3i', 'he_no4i', 'he_no5', 'he_no5i'
            ],
            8,
            fn='INOB',
            indep=+1)
        os.system(BINBDGpath + 'KOBER.exe')
        if os.path.isdir(v18uixpath + 'eob-tni/') == False:
            os.mkdir(v18uixpath + 'eob-tni/')
        os.chdir(v18uixpath + 'eob-tni/')
        n3_inob(
            [
                'he_no1', 'he_no1i', 'he_no2', 'he_no2i', 'he_no6', 'he_no6i',
                'he_no3', 'he_no3i', 'he_no4i', 'he_no5', 'he_no5i'
            ],
            15,
            fn='INOB',
            indep=+1)
        os.system(BINBDGpath + 'DROBER.exe')
        if os.path.isdir(v18uixpath + 'elu/') == False:
            os.mkdir(v18uixpath + 'elu/')
        os.chdir(v18uixpath + 'elu/')
        n3_inlu(
            8,
            fn='INLUCN',
            fr=[
                '000', '202', '022', '110', '101', '011', '111', '112', '211',
                '212', '213', '121', '122', '212', '222', '221', '220'
            ],
            indep=+1)
        os.system(BINBDGpath + 'LUDW_CN.exe')
        if os.path.isdir(v18uixpath + 'elu-tni/') == False:
            os.mkdir(v18uixpath + 'elu-tni/')
        os.chdir(v18uixpath + 'elu-tni/')
        n3_inlu(
            8,
            fn='INLU',
            fr=[
                '000', '202', '022', '110', '101', '011', '111', '112', '211',
                '212', '213', '121', '122', '212', '222', '221', '220'
            ],
            indep=+1)
        os.system(BINBDGpath + 'DRLUD.exe')

    os.chdir(v18uixpath)
    print('Calculating in %s' % v18uixpath)
    n3_inlu(8, fn='INLU', fr=lfrags2, indep=-1)
    os.system(BINBDGpath + 'DRLUD.exe')
    n3_inlu(8, fn='INLUCN', fr=lfrags2, indep=-1)
    os.system(BINBDGpath + 'LUDW_CN.exe')
    n3_inob(sfrags2, 8, fn='INOB', indep=-1)
    os.system(BINBDGpath + 'KOBER.exe')
    n3_inob(sfrags2, 15, fn='INOB', indep=-1)
    os.system(BINBDGpath + 'DROBER.exe')

    he3inqua(intwi=widi, relwi=widr, potf=potnn)

    parallel_mod_of_3inqua(
        lfrags2,
        sfrags2,
        infile='INQUA_N',
        outfile='INQUA_N',
        einzel_path=v18uixpath)

    insam(len(lfrags2))

    n3_inen_bdg(sbas, Jstreu, costr, fn='INEN', pari=0, nzop=31, tni=11)

    if bastype == boundstatekanal:
        n3_inlu(8, fn=helionpath + 'INLU_ref', fr=lfrags2, indep=-1)
        n3_inlu(8, fn=helionpath + 'INLUCN_ref', fr=lfrags2, indep=-1)
        n3_inob(sfrags2, 8, fn=helionpath + 'INOB_ref', indep=-1)
        n3_inob(sfrags2, 15, fn=helionpath + 'DRINOB_ref', indep=-1)
        os.system('cp INQUA_N ' + helionpath + 'INQUA_V18_ref')
        os.system('cp INEN ' + helionpath + 'INEN_ref')
        os.system('cp INSAM ' + helionpath)

    subprocess.run(
        ['mpirun', '-np',
         '%d' % anzproc, BINBDGpath + 'V18_PAR/mpi_quaf_v6'])
    subprocess.run([BINBDGpath + 'V18_PAR/sammel'])

    he3inqua(intwi=widi, relwi=widr, potf=potnnn)

    parallel_mod_of_3inqua(
        lfrags2,
        sfrags2,
        infile='INQUA_N',
        outfile='INQUA_N',
        tni=1,
        einzel_path=v18uixpath)

    subprocess.run([
        'mpirun', '-np',
        '%d' % anzproc, BINBDGpath + 'UIX_PAR/mpi_drqua_uix'
    ])
    subprocess.run([BINBDGpath + 'UIX_PAR/SAMMEL-uix'])
    subprocess.run([BINBDGpath + 'TDR2END_AK.exe'])

    suche_fehler()

    if bastype == boundstatekanal:
        os.system('cp INQUA_N ' + helionpath + 'INQUA_UIX_ref')

    ew_threshold = 10**(-9)

    subprocess.call('cp MATOUT mat_%s' % bastype, shell=True)
    matout = [line for line in open('MATOUT')]

    goodEVs = smart_ev(matout, ew_threshold)

    print('the good guys: ', np.real(goodEVs[:4]), ' ... ',
          np.real(goodEVs[-4:]))

    subprocess.call('rm ' + v18uixpath + 'TQUAOUT.*', shell=True)
    subprocess.call('rm ' + v18uixpath + 'TDQUAOUT.*', shell=True)
    subprocess.call('rm ' + v18uixpath + 'DMOUT.*', shell=True)
    subprocess.call('rm ' + v18uixpath + 'DRDMOUT.*', shell=True)