from bridge import *
from three_particle_functions import *
from triton_width_gen import *
from parameters_and_constants import *

from matplotlib import cm

from sklearn.datasets import make_blobs
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import brier_score_loss
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import train_test_split

bastype = streukas[0]  #boundstatekanal  #

angu = channels[bastype]

Jay = 0.5
Jstreustring = '%s' % str(Jay)

if os.path.isdir(litpath3He) == False:
    os.mkdir(litpath3He)
os.chdir(litpath3He)

if os.path.isdir(litpath3He + 'basis_struct/') == False:
    os.mkdir(litpath3He + 'basis_struct/')
os.chdir(litpath3He + 'basis_struct/')

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

# total number of (internal,relative) pairs"The solid earth
int_samples = 500
rel_samples = 500

int_centers = [np.random.uniform(0.01, 1.3) for n in range(len(lfrags))]
rel_centers = [np.random.uniform(0, 0.003) for n in range(len(lfrags))]

I, label_i = make_blobs(
    n_samples=int_samples,
    #    centers=int_centers,
    shuffle=False,
    random_state=42,
    n_features=1)

R, label_r = make_blobs(
    n_samples=rel_samples,
    #    centers=rel_centers,
    shuffle=False,
    random_state=42,
    n_features=1)

sample_weight_i = np.random.RandomState(42).rand(label_i.shape[0])
sample_weight_r = np.random.RandomState(42).rand(label_r.shape[0])

# split train, test for calibration
int_train, int_test, lab_int_train, lab_int_test, sw_i_train, sw_i_test = \
    train_test_split(I, label_i, sample_weight_i, test_size=0.9, random_state=42)
rel_train, rel_test, lab_rel_train, lab_rel_test, sw_r_train, sw_r_test = \
    train_test_split(R, label_r, sample_weight_r, test_size=0.9, random_state=42)

# Plot the data
plt.figure()
y_unique = np.unique(label_i)

colors = cm.rainbow(np.linspace(0.0, 1.0, y_unique.size))

for this_y, color in zip(y_unique, colors):
    this_X = int_train[lab_int_train == this_y]
    this_sw_i = sw_i_train[lab_int_train == this_y]
    this_Y = rel_train[lab_rel_train == this_y]
    this_sw_r = sw_r_train[lab_rel_train == this_y]
    plt.scatter(
        np.abs(this_X),
        np.abs(this_Y),
        #        s=this_sw * 50,
        c=color[np.newaxis, :],
        alpha=0.5,
        edgecolor='k',
        label="Fragment %s" % this_y)
plt.legend(loc="best")
plt.title("Data")

plt.show()