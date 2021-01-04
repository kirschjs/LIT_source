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
from scipy.spatial.distance import *

bastype = boundstatekanal
bastype = streukas[0]

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

fragfile = [
    ln for ln in open(litpath3He + 'basis_struct/frags_LIT_J%s_%s.dat' %
                      (Jstreustring, bastype))
]

lfrags = [fr.split(' ')[1].strip() for fr in fragfile]
sfrags = [fr.split(' ')[0] for fr in fragfile]

lit_iw = [
    np.array(ln.split(';')).astype(float)
    for ln in open(litpath3He + 'basis_struct/intw3heLIT_J%s_%s.dat' %
                   (Jstreustring, bastype))
]

anzLITbv = sum([len(frgm) for frgm in lit_iw])
if (len([len(frgm) for frgm in lit_iw]) != (len(fragfile))):
    print('LIT-basis fragments inconcistent!',
          len([len(frgm) for frgm in lit_iw]), (len(fragfile)))
    exit()

lit_rw = [
    np.array(ln.split(';')).astype(float)
    for ln in open(litpath3He + 'basis_struct/relw3heLIT_J%s_%s.dat' %
                   (Jstreustring, bastype))
]

anzBV = sum([len(zer) for zer in lit_iw])
print('# Basisvektoren = %d' % anzBV)
sbas = []
bvn = 1
bvv = 1
wgri = []
wgridx = []
for zerl in range(len(lit_iw)):
    for bv in range(len(lit_iw[zerl])):
        for rw in range(len(lit_rw[zerl])):
            wgri += [[lit_iw[zerl][bv], lit_rw[zerl][rw]]]
            wgridx += [[bvn, bvv, rw]]
            bvn += 1
        bvv += 1

IDXs = []


def push_back(idx, dist):
    IDXs.append(idx)


#metrics = ['braycurtis’, 'canberra’, 'chebyshev’, 'cityblock’, 'correlation’, 'cosine’, 'dice’, 'euclidean’, 'hamming’, 'jaccard’, 'jensenshannon’, 'kulsinski’, 'mahalanobis’, 'matching’, 'minkowski’, 'rogerstanimoto’, 'russellrao’, 'seuclidean’, 'sokalmichener’, 'sokalsneath’, 'sqeuclidean’, 'wminkowski’, 'yule’]
#def mat_dist(MA, MB):
#    return cdist(MA, MB, metric=metrics[1])  #'seuclidean')  #

ns = int(len(wgri) / 14)
wgri_strat = subset.select_greedy_maxisum(wgri, ns, callback=push_back)
#wgri_strat = subset.select_greedy_maximin(wgri, ns, dist_matrix_function=mat_dist, callback=push_back)
#wgri_strat = subset.select_greedy_energy(
#    wgri, ns, dist_matrix_function=mat_dist, callback=push_back)
wgri_strat2 = np.take(wgridx, IDXs[-1], axis=0)

fig = plt.figure(figsize=[12, 6])
#fig.subplots_adjust(hspace=0.4, wspace=0.4)
ax = plt.subplot(111)

ax.plot(
    np.array(wgri_strat)[:, 0],
    np.array(wgri_strat)[:, 1],
    'o',
    alpha=1.,
    color='black',
    label=r'reduced basis set',
    marker='.',
    markersize=14)

ax.plot(
    np.array(wgri)[:, 0],
    np.array(wgri)[:, 1],
    'o',
    alpha=1.,
    color='orange',
    label=r'full basis set',
    marker='.',
    markersize=6)

plt.xlabel(r'$\gamma(\rho_1)\;\;\; [fm^{-2}]$')
plt.ylabel(r'$\gamma(\rho_2)\;\;\; [fm^{-2}]$')

plt.legend(loc='best')

fig.savefig('WidthXY.pdf')

sbas = np.array(wgri_strat2)
basi = {}
for bv in np.unique(sbas[:, 1]):
    basi[bv] = [bv, []]
    for bvv in sbas[:, 1:]:
        if bvv[0] == bv:
            basi[bv][1].append(bvv[1])

sbas = list(basi.values())

if sbas == []:
    print('no width combination permissible.')
    exit()
os.chdir(v18uixpath)
n3_inen_bdg(sbas, Jstreu, costr, fn='INEN', pari=0, nzop=31, tni=11)

subprocess.run([home + '/kette_repo/ComptonLIT/src_nucl/TDR2END_AK.exe'])

ew_threshold = 6**(-5)

matout = [line for line in open('MATOUT')]

goodEVs = smart_ev(matout, ew_threshold)

print('the good guys: ', np.real(goodEVs[:4]), ' ... ', np.real(goodEVs[-4:]))