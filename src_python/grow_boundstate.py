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

os.chdir(v18uixpath)

inen = [line for line in open('INEN')]

bastype = streukas[0]
bastype = boundstatekanal

Jstreu = float(bastype.split('^')[0][-3:])
Jstreustring = '%s' % str(Jstreu)[:3]

ite = 0
Eth = -7.93
#while ((abs(EBDG + 2.22) > 0.05) & (ite < 10)):
while (ite < 1):
    ite += 1
    # read width sets of the ground-state (3-helium) basis
    he_iw, he_rw, frgs = retrieve_he3_widths('INQUA_N')
    if len(he_iw) != len(he_rw):
        print('anz rw sets > anz iw sets.')
        exit()
    bvnr = 1
    bvs = []
    for n in range(len(he_iw)):
        for bv in he_iw[n]:
            bvs.append([bvnr, len(he_rw[n])])
            bvnr += 1
    #np.random.shuffle(bvs)
    sout = ''
    for line in inen[:8]:
        sout += line
    bvnr = 1
    for bv in bvs:
        rwunglz = False
        soutmp = sout + '   1%4d\n' % bv[0]
        rwnr = 1
        for relw in range(bv[1]):
            soutmp += '%3d' % np.random.randint(0, 2)
            with open('INEN', 'w') as outfile:
                outfile.write(soutmp)
            repl_line('INEN', 7,
                      '%4s%4d%4s%4s%4s\n' % (inen[7][:4], bvnr, inen[7][8:12],
                                             inen[7][12:16], inen[7][16:20]))
            subprocess.run(
                [home + '/kette_repo/ComptonLIT/src_nucl/TDR2END_AK.exe'])
            if bvnr != 1:
                diag = open('OUTPUT',
                            'r').read().find('DIAGONALISIERUNG FEHLERHAFT')
                if ((float(get_h_ev(4, 'OUTPUT')[0]) < Eth) | (diag > 0)):
                    soutmp = soutmp[:-3] + '  0'
                else:
                    print(get_h_ev(4, 'OUTPUT')[0])
                    rwunglz = True
            else:
                rwunglz = True
        if rwunglz:
            soutmp += '\n'
            sout = soutmp
            bvnr += 1
    exit()