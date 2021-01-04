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

        # generate widths for the two coordinates for every LIT-basis spin-/angular-momentum comfiguration
        lit_w = {}
        basisDim = 0
        # ONE : sample width set for all fragments
        for frg in range(len(lfrags)):

            contiwidths1 = []
            contiwidths2 = []
            bachanwidths1 = []
            bachanwidths2 = []

            nbrwidthsets = 1

            X, y = make_blobs(
                n_samples=150,
                n_features=2,
                centers=50,
                center_box=(0, np.random.uniform(0.000001, 0.000005)),
                cluster_std=np.random.uniform(0, 0.01),
                shuffle=True,
                random_state=None)

            km = KMeans(
                n_clusters=cluster_centers_per_zerl,
                init='random',
                n_init=10,
                max_iter=300,
                tol=1e-04,
                random_state=0)

            X = np.abs(X)

            y_km = km.fit_predict(X)

            plotti = 0
            if plotti:
                # plot the 3 clusters
                [
                    plt.scatter(
                        X[y_km == nc, 0],
                        X[y_km == nc, 1],
                        s=50,
                        marker='s',
                        edgecolor='black',
                        label='cluster %d' % nc)
                    for nc in range(cluster_centers_per_zerl)
                ]

                # plot the centroids
                plt.scatter(
                    km.cluster_centers_[:, 0],
                    km.cluster_centers_[:, 1],
                    s=250,
                    marker='*',
                    c='red',
                    edgecolor='black',
                    label='centroids')
                plt.legend(scatterpoints=1)
                plt.grid()
                plt.show()

            lit_w[frg] = km.cluster_centers_
            basisDim += len(lit_w[frg])

        # TWO: global distancing over all bvs
        lit_w_sparse = {}

        for frg in range(len(lit_w)):
            bv_init = np.random.randint(len(lit_w[frg]))
            lit_w_sparse[frg] = lit_w[frg][bv_init]
            lit_w[frg] = np.delete(lit_w[frg], bv_init, 0)

        print('LIT basis dim = %d' % basisDim)

        for ntest in range(basisDim):
            # select random bv
            frg_test = np.random.choice(list(lit_w.keys()))
            bv_test = np.random.randint(len(lit_w[frg_test]))

            reject = False
            # check its distance to each vector in the clean set
            for frg in range(len(lit_w_sparse)):
                for cleanbv in range(len(lit_w_sparse[frg])):
                    reject = True if np.linalg.norm(
                        lit_w_sparse[frg][cleanbv] - lit_w[frg_test][bv_test]
                    ) < min_eucl_pair_dist else False
                    if reject:
                        break
                if reject:
                    break

            # remove it from the old basis
            if reject == False:
                lit_w_sparse[frg_test] = np.append(
                    lit_w_sparse[frg_test], np.array(lit_w[frg_test][bv_test]))
            lit_w[frg_test] = np.delete(lit_w[frg_test], bv_test, 0)
            if len(lit_w[frg_test]) == 0:
                del lit_w[frg_test]
            if list(lit_w.keys()) == []:
                break

        lit_iw = []
        lit_rw = []

        bvBAADz = 12

        for n in range(len(lit_w_sparse)):
            sfrags2 += int(np.ceil(
                (len(lit_w_sparse[n]) / 2.) / bvBAADz)) * [sfrags[n]]
            lfrags2 += int(np.ceil(
                (len(lit_w_sparse[n]) / 2.) / bvBAADz)) * [lfrags[n]]

            lit_iw += np.array(
                np.array_split(
                    lit_w_sparse[n][0::2],
                    int(np.ceil(
                        (len(lit_w_sparse[n]) / 2.) / bvBAADz)))).tolist()

            lit_rw += np.array(
                np.array_split(
                    lit_w_sparse[n][1::2],
                    int(np.ceil(
                        (len(lit_w_sparse[n]) / 2.) / bvBAADz)))).tolist()

        sbas = []
        bv = 1
        for n in range(len(lfrags2)):
            for m in range(len(lit_rw[n])):
                sbas += [[bv, [m + 1]]]
                bv += 1

        if sbas == []:
            print('no width combination permissible.')
            exit()

        print(lfrags2)

        with open('LITbas_full_J%s.dat' % Jstreustring, 'wb') as f:
            np.savetxt(
                f, [[jj[0], kk] for jj in sbas for kk in jj[1]], fmt='%d')
        f.close()
        os.system('cp LITbas_full_J%s.dat LITbas_red_J%s.dat' % (Jstreustring,
                                                                 Jstreustring))

        with open(litpath3He + 'frags_LIT_J%s.dat' % Jstreustring, 'wb') as f:
            np.savetxt(
                f,
                np.column_stack([sfrags2, lfrags2]),
                fmt='%s',
                delimiter=' ',
                newline=os.linesep)
        f.close()

        with open(litpath3He + 'intw3heLIT_J%s.dat' % Jstreustring, 'wb') as f:
            for ws in lit_iw:
                np.savetxt(f, [ws], fmt='%12.6f', delimiter=' ; ')
        f.close()

        with open(litpath3He + 'relw3heLIT_J%s.dat' % Jstreustring, 'wb') as f:
            for ws in lit_rw:
                np.savetxt(f, [ws], fmt='%12.6f', delimiter=' ; ')
        f.close()

        visbas(
            basispath=litpath3He + 'lit_bas_rhs/',
            widthpath=litpath3He,
            exepath=BINBDGpath,
            Jstrstr=Jstreustring)

        print('LIT basis: DIM = %d   #FRGM = %d' % (len(sbas), len(lfrags2)))

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