from bridge import *
from three_particle_functions import *
from triton_width_gen import *
from parameters_and_constants import *
import operator
from BasisVisualization import visbas
from plot_basis_3bdy import plotbasis3bdy
from plot_spectrum import plotHspec
from scipy.cluster.vq import vq, kmeans, whiten
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.cluster import KMeans

use_old_widths = 0

if os.path.isdir(litpath3He) == False:
    os.mkdir(litpath3He)
    os.mkdir(litpath3He + 'results/')

if os.path.isdir(litpath3He + 'lit_bas_rhs/') == False:
    os.mkdir(litpath3He + 'lit_bas_rhs/')
os.chdir(litpath3He + 'lit_bas_rhs/')

for streukanal in streukas:

    Jstreu = float(streukanal.split('^')[0])
    Jstreustring = '%s' % str(Jstreu)[:3]

    ite = 0
    EBDG = 666
    #while ((abs(EBDG + 2.22) > 0.05) & (ite < 10)):

    while ((EBDG == 666) & (ite < 1)):
        ite += 1

        lfrags = []
        sfrags = []
        lfrags2 = []
        sfrags2 = []

        for lcfg in range(len(streukanaele3He[streukanal])):
            sfrags = sfrags + streukanaele3He[streukanal][lcfg][1]
            for scfg in streukanaele3He[streukanal][lcfg][1]:
                lfrags = lfrags + [streukanaele3He[streukanal][lcfg][0]]

        print('\n3-body LIT-basis components (bare):%5d :' % len(lfrags))

        # read width sets of the ground-state (3-helium) basis
        he_iw, he_rw, frgs = retrieve_he3_widths(v18uixpath + 'INQUA_N_UIX')

        if use_old_widths:

            litbas_full = np.loadtxt(
                'LITbas_full_J%s.dat' % Jstreu).astype(int)
            try:
                litbas_red = np.loadtxt(
                    'LITbas_red_J%s.dat' % Jstreu).astype(int)
            except:
                litbas_red = litbas_full

            lit_iw = [
                np.array(ln.split(';')).astype(float).tolist()
                for ln in open(litpath3He + 'intw3heLIT_J%s.dat' % Jstreu)
            ]

            lit_rw = [
                np.array(ln.split(';')).astype(float).tolist()
                for ln in open(litpath3He + 'relw3heLIT_J%s.dat' % Jstreu)
            ]

            frg = 0
            for intwn in range(len(lit_iw)):
                if lit_iw[intwn][0] <= lit_iw[np.max([0, intwn - 1])][0]:
                    sfrags2 += [sfrags[frg]]
                    lfrags2 += [lfrags[frg]]
                else:
                    frg += 1
                    sfrags2 += [sfrags[frg]]
                    lfrags2 += [lfrags[frg]]

        else:

            # generate widths for the two coordinates for every LIT-basis spin-/angular-momentum comfiguration

            lit_iw = []

            lit_rw = []

            for frg in range(len(lfrags)):

                contiwidths1 = []
                contiwidths2 = []
                bachanwidths1 = []
                bachanwidths2 = []

                nbrwidthsets = 2
                for gaussset in range(nbrwidthsets):

                    shaker = (1. - 2. * np.random.random()
                              ) * 0.5 * np.array(gauss_center)
                    bachanwidths1 = [
                        w
                        for w in np.random.normal(gauss_center[0] + shaker[0],
                                                  gauss_sigma[0], dimC[0])
                        if ((np.abs(w) > eps_low[0]) & (np.abs(w) < eps_up[0]))
                    ]
                    bachanwidths2 = [
                        w
                        for w in np.random.normal(gauss_center[1] + shaker[1],
                                                  gauss_sigma[1], dimC[1])
                        if ((np.abs(w) > eps_low[1]) & (np.abs(w) < eps_up[1]))
                    ]

                    contiwidths1 = np.sort(
                        np.abs(np.concatenate((contiwidths1,
                                               bachanwidths1))))[::-1]
                    contiwidths2 = np.sort(
                        np.abs(np.concatenate((contiwidths2,
                                               bachanwidths2))))[::-1]

                iw = contiwidths1  #sparsify(rho1widths, min_single_width_dist_int)  #
                rw = contiwidths2  #sparsify(rho2widths, min_single_width_dist_rel)  #

                #widthdata = list(zip(iw, rw))
                widthdata = np.array([[a, b] for b in iw for a in rw])
                #print(iw, len(rw), np.shape(widthdata))

                km = KMeans(n_clusters=maxdim).fit(widthdata)
                kmlabels = km.labels_

                fig = plt.figure()
                ax = fig.add_axes([0, 0, 1, 1])
                ax.scatter(
                    widthdata[:, 0],
                    widthdata[:, 1],
                    c=kmlabels.astype(np.float),
                    edgecolor='k')
                ax.scatter(
                    km.cluster_centers_[:, 0],
                    km.cluster_centers_[:, 1],
                    c='red',
                    edgecolor='k')
                #    plt.show()

                #ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)

                #bandwidth = estimate_bandwidth(
                #    widthdata, quantile=cluster_quanil, n_samples=15000)
                #ms.fit(widthdata)
                #labels = ms.labels_

                cluster_centers = km.cluster_centers_

                labels_unique = np.unique(kmlabels)
                n_clusters_ = len(labels_unique)

                #print("number of estimated clusters : %d" % n_clusters_)

                #print('number of formed clusters : %d' % len(cluster_centers))
                if len(cluster_centers) < 2:
                    print('no clusters found.')
                    exit()

                iw = np.sort(cluster_centers[:, 0])[::-1]
                rw = np.sort(cluster_centers[:, 1])[::-1]

                lit_iw += np.array(
                    np.array_split(iw, int(np.ceil(len(iw) / 12.)))).tolist()
                lit_rw += np.array(
                    np.array_split(rw, int(np.ceil(len(rw) / 12.)))).tolist()

                #lit_rw += int(np.ceil(len(iw) / 12.)) * rw.reshape(
                #    (1, -1)).tolist()

                #for w in rho1widths[::-1]:
                #    print('%8.4f' % w, end='')
                #print('   ---- ', end='')
                #for w in rho2widths[::-1]:
                #    print('%8.4f' % w, end='')
                #print('')
                #print(lit_iw)
                #print(lit_rw)

                sfrags2 += int(np.ceil(len(iw) / 12.)) * [sfrags[frg]]
                lfrags2 += int(np.ceil(len(iw) / 12.)) * [lfrags[frg]]

            with open(litpath3He + 'frags_LIT_J%s.dat' % Jstreustring,
                      'wb') as f:
                np.savetxt(
                    f,
                    np.column_stack([sfrags2, lfrags2]),
                    fmt='%s',
                    delimiter=' ',
                    newline=os.linesep)
            f.close()

            with open(litpath3He + 'intw3heLIT_J%s.dat' % Jstreustring,
                      'wb') as f:
                for ws in lit_iw:
                    np.savetxt(f, [ws], fmt='%12.6f', delimiter=' ; ')
            f.close()
            with open(litpath3He + 'relw3heLIT_J%s.dat' % Jstreustring,
                      'wb') as f:
                for ws in lit_rw:
                    np.savetxt(f, [ws], fmt='%12.6f', delimiter=' ; ')
            f.close()

            bas = []
            dists = []
            fragstruct = {}
            print(len(lit_iw), len(lit_rw))

            iw_dist = []
            rw_dist = []

            for nz in range(len(lit_iw)):

                iw_dist.append([])
                rw_dist.append([])

                allindices = [[w1, w2] for w1 in range(len(lit_iw[nz]))
                              for w2 in range(len(lit_rw[nz]))]
                diagindices = [[w1, w1] for w1 in range(len(lit_iw[nz]))]
                allwidths = [[lit_iw[nz][ind[0]], lit_rw[nz][ind[1]]]
                             for ind in allindices]

                n0 = np.random.choice(range(len(allindices)))
                set0 = [allindices[n0]]
                allindices.remove(allindices[n0])

                while len(allindices) > 0:
                    reject = False
                    n0 = np.random.choice(range(len(allindices)))
                    e0 = allindices[n0]

                    w1 = [lit_iw[nz][e0[0]], lit_rw[nz][e0[1]]]

                    reject = True if (np.linalg.norm(w1) < eps_D) else False

                    for ntest in range(len(set0)):

                        w2 = [
                            lit_iw[nz][set0[ntest][0]],
                            lit_rw[nz][set0[ntest][1]]
                        ]
                        normdist = np.linalg.norm(np.array(w1) - np.array(w2))
                        if ((normdist < min_eucl_pair_dist) |
                            (len(w1 + w2) != len(np.flatnonzero(w1 + w2)))):

                            reject = True
                            #        exit()
                            break

                    if reject:
                        allindices.remove(allindices[n0])
                        continue
                    else:
                        set0.append(allindices[n0])
                        iw_dist[-1].append(w1)
                        rw_dist[-1].append(w2)

                    allindices.remove(e0)

                offset_int = sum([len(iws) for iws in lit_iw[:nz]])
                bas += [[el[0] + 1 + offset_int, el[1] + 1] for el in set0]

                for w1 in range(len(lit_iw[nz])):
                    fragstruct[str(w1 + 1 + offset_int)] = len(lit_rw[nz])

            print(iw_dist)
            exit()
            basvs = {}
            for bv in bas:
                try:
                    basvs[bv[0]].append(bv[1])
                except:
                    basvs[bv[0]] = [bv[1]]

            sbas = sorted(basvs.items(), key=operator.itemgetter(0))

            sbas = [[basv[0], sorted(basv[1])] for basv in sbas]

            sbasdiagonal = []

            #            bins = [0, dimC[0]]
            #
            #            for bvn in range(len(sbas)):
            #
            #                rels = [
            #                    relw for nn in range(1, len(bins)) for relw in sbas[bvn][1]
            #                    if (((relw > bins[nn - 1]) & (sbas[bvn][0] > bins[nn - 1]))
            #                        & ((relw <= bins[nn]) & (sbas[bvn][0] < bins[nn])))
            #                ]
            #
            #                relc = sbas[bvn][1][-1 + sbas[bvn][0] % len(sbas[bvn][1])]
            #                sbasdiagonal.append([
            #                    sbas[bvn][0], [np.min([relc, relc + 1])]
            #                    #list(
            #                    #    np.unique([
            #                    #        np.max([1, relc - 1]), relc,
            #                    #        np.min([len(sbas[bvn][1]), relc + 1])
            #                    #    ]))
            #                ])
            #                #if (sbas[bvn][0] < (0.5 * dimC1)):
            #                #    sbasdiagonal.append(
            #                #        [sbas[bvn][0], [int(sbas[bvn][0] + 0.5 * dimC1)]])
            #                #if sbas[bvn][0] > 0.5 * dimC1:
            #                #    sbasdiagonal.append(
            #                #        [sbas[bvn][0], [sbas[bvn][0] - int(0.5 * dimC1)]])
            #
            #            sbas = sbas  #sbasdiagonal  #

            if sbas == []:
                print('no width combination permissible.')
                exit()

            with open('LITbas_full_J%s.dat' % Jstreustring, 'wb') as f:
                np.savetxt(
                    f, [[jj[0], kk] for jj in sbas for kk in jj[1]], fmt='%d')
            f.close()

            os.system('cp LITbas_full_J%s.dat LITbas_red_J%s.dat' %
                      (Jstreustring, Jstreustring))

        visbas(
            basispath=litpath3He + 'lit_bas_rhs/',
            widthpath=litpath3He,
            exepath=BINBDGpath,
            Jstrstr=Jstreustring)

        exit()

        n3_inlu(21, fn='INLU', fr=lfrags2)
        os.system(BINBDGpath + 'DRLUD.exe')
        n3_inlu(21, fn='INLUCN', fr=lfrags2)
        os.system(BINBDGpath + 'LUDW_CN.exe')
        n3_inob(sfrags2, 20, fn='INOB')
        os.system(BINBDGpath + 'KOBER.exe')
        os.system(BINBDGpath + 'DROBER.exe')

        he3inqua(
            intwi=lit_iw,
            relwi=lit_rw,
            potf='/home/kirscher/kette_repo/sim_par/potentials/NN_pheno/AV18')

        os.system('time ' + BINBDGpath + 'QUAFL_N.exe')
        repl_line(
            'INQUA_N', 1,
            '/home/kirscher/kette_repo/sim_par/potentials/NNN_pheno/urbana9_AK_neu\n'
        )
        os.system('time ' + BINBDGpath + 'DRQUA_AK_N.exe')

        os.system('cp QUAOUT QUAOUT_tmp')
        os.system('cp DRQUAOUT DRQUAOUT_tmp')

        n3_inen_bdg(sbas, Jstreu, costr, fn='INEN', pari=0, nzop=31, tni=11)

        os.system(
            BINBDGpath + 'DR2END_AK.exe && grep -A 3 \'EIGENWER\' OUTPUT')

        cmd = "grep -c 'NaN' OUTPUT"
        nans = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        anznan = int(
            str(nans.communicate()[0]).strip().split("'")[1].split("\\")[0])
        print('%d NaN occurences in OUTPUT' % anznan)

        EBDG = get_h_ev()[0] if anznan == 0 else 666
        print('%4d: %4.4f' % (ite, EBDG))

        if use_old_widths:
            break

    curdir = os.getcwd()
    if EBDG != 666:
        plotHspec(
            inf=curdir + '/OUTPUT',
            outf=litpath3He + 'results/h_spec_%s.pdf' % suffix)

    os.system('cp QUAOUT QUAOUT_J%s' % Jstreustring)
    os.system('cp DRQUAOUT DRQUAOUT_J%s' % Jstreustring)
    print('QUAOUT_J%s saved to %s' % (Jstreustring, curdir))