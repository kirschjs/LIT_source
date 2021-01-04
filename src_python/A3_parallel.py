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

# define orbital [(iso)spin/flavor] structure

angu = [['000', ['he_no1', 'he_no6']], ['202', ['he_no2', 'he_no2']],
        ['022', ['he_no2']], ['111', ['he_no3']], ['112', ['he_no5']],
        ['111', ['he_no5']]]

# generate radial basis for each of the above components
basisDim = 0
lit_w = {}
lit_iw = []
lit_rw = []
lfrags = []
sfrags = []
lfrags2 = []
sfrags2 = []

for lcfg in range(len(angu)):
    sfrags = sfrags + angu[lcfg][1]
    for scfg in angu[lcfg][1]:
        lfrags = lfrags + [angu[lcfg][0]]

cluster_centers_per_zerl = 24
for frg in range(len(angu)):

    contiwidths1 = []
    contiwidths2 = []
    bachanwidths1 = []
    bachanwidths2 = []
    nbrwidthsets = 1

    X, y = make_blobs(
        n_samples=150,
        n_features=2,
        centers=50,
        center_box=(0, np.random.uniform(0.05, 2.5)),
        cluster_std=np.random.uniform(0, 0.7),
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
        lit_w_sparse[frg_test] = np.append(lit_w_sparse[frg_test],
                                           np.array(lit_w[frg_test][bv_test]))
    lit_w[frg_test] = np.delete(lit_w[frg_test], bv_test, 0)
    if len(lit_w[frg_test]) == 0:
        del lit_w[frg_test]
    if list(lit_w.keys()) == []:
        break

lit_iw = []
lit_rw = []
bvBAADz = 8
for n in range(len(lit_w_sparse)):
    sfrags2 += int(np.ceil(
        (len(lit_w_sparse[n]) / 2.) / bvBAADz)) * [sfrags[n]]
    lfrags2 += int(np.ceil(
        (len(lit_w_sparse[n]) / 2.) / bvBAADz)) * [lfrags[n]]

    lit_iw += np.sort(
        np.array(
            np.array_split(
                lit_w_sparse[n][0::2],
                int(np.ceil(
                    (len(lit_w_sparse[n]) / 2.) / bvBAADz))))).tolist()

    lit_iw[-1] = lit_iw[-1][::-1]

    lit_rw += np.sort(
        np.array(
            np.array_split(
                lit_w_sparse[n][1::2],
                int(np.ceil(
                    (len(lit_w_sparse[n]) / 2.) / bvBAADz))))).tolist()
    lit_rw[-1] = lit_rw[-1][::-1]

anzBV = sum([len(zer) for zer in lit_iw])

print('# Basisvektoren = %d' % anzBV)
sbas = []
bvn = 1
for zerl in range(len(lit_iw)):
    for bv in range(len(lit_iw[zerl])):
        sbas += [[
            bv + bvn,
            np.array([(bv + n) % 2 for n in range(1 + len(lit_rw[zerl]))]) *
            np.arange(1 + len(lit_rw[zerl]))
        ]]
    bvn += bv + 1

if sbas == []:
    print('no width combination permissible.')
    exit()

print(lfrags2)
print(sfrags2)

os.chdir('/home/kirscher/kette_repo/ComptonLIT/he3/eob/')
n3_inob(
    ['he_no1', 'he_no2', 'he_no6', 'he_no3', 'he_no5'], 8, fn='INOB', indep=+1)
os.system(BINBDGpath + 'KOBER.exe')
os.chdir('/home/kirscher/kette_repo/ComptonLIT/he3/eob-tni/')
n3_inob(
    ['he_no1', 'he_no2', 'he_no6', 'he_no3', 'he_no5'],
    15,
    fn='INOB',
    indep=+1)
os.system(BINBDGpath + 'DROBER.exe')

os.chdir('/home/kirscher/kette_repo/ComptonLIT/he3/elu/')
n3_inlu(
    8,
    fn='INLUCN',
    fr=['000', '202', '022', '110', '101', '011', '111', '112'],
    indep=+1)
os.system(BINBDGpath + 'LUDW_CN.exe')
os.chdir('/home/kirscher/kette_repo/ComptonLIT/he3/elu-tni/')
n3_inlu(
    8,
    fn='INLU',
    fr=['000', '202', '022', '110', '101', '011', '111', '112'],
    indep=+1)
os.system(BINBDGpath + 'DRLUD.exe')

os.chdir('/home/kirscher/kette_repo/ComptonLIT/he3/')

n3_inlu(8, fn='INLU', fr=lfrags2, indep=-1)
os.system(BINBDGpath + 'DRLUD.exe')
n3_inlu(8, fn='INLUCN', fr=lfrags2, indep=-1)
os.system(BINBDGpath + 'LUDW_CN.exe')

n3_inob(sfrags2, 8, fn='INOB', indep=-1)
os.system(BINBDGpath + 'KOBER.exe')
n3_inob(sfrags2, 15, fn='INOB', indep=-1)
os.system(BINBDGpath + 'DROBER.exe')

insam(len(lfrags2))
n3_inen_bdg(sbas, 0.5, costr, fn='INEN', pari=0, nzop=31, tni=11)

he3inqua(
    intwi=lit_iw,
    relwi=lit_rw,
    potf='/home/kirscher/kette_repo/sim_par/potentials/NN_pheno/AV18')

parallel_mod_of_3inqua(lfrags2, sfrags2, infile='INQUA_N', outfile='INQUA_N')

subprocess.run([
    'mpiexec', '-np', '6',
    '/home/kirscher/kette_repo/ComptonLIT/src_nucl/V18_PAR/mpi_quaf_v6'
])

subprocess.run(
    ['/home/kirscher/kette_repo/ComptonLIT/src_nucl/V18_PAR/sammel'])

he3inqua(
    intwi=lit_iw,
    relwi=lit_rw,
    potf='/home/kirscher/kette_repo/sim_par/potentials/NNN_pheno/urbana9_AK_neu'
)

parallel_mod_of_3inqua(
    lfrags2, sfrags2, infile='INQUA_N', outfile='INQUA_N', tni=1)

subprocess.run([
    'mpiexec', '-np', '6',
    '/home/kirscher/kette_repo/ComptonLIT/src_nucl/UIX_PAR/mpi_drqua_uix'
])

subprocess.run(
    ['/home/kirscher/kette_repo/ComptonLIT/src_nucl/UIX_PAR/SAMMEL-uix'])

subprocess.run(
    ['/home/kirscher/kette_repo/ComptonLIT/src_nucl/TDR2END_I_2.exe'])

os.system('grep -A 3 \'EIGENWER\' OUTPUT')
exit()