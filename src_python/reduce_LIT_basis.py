from bridge import *
from three_particle_functions import *
from BasisVisualization import visbas
import operator

if os.path.isdir(litpath3He + 'lit_bas_rhs/') == False:
    print('no basis set, yet. run <v18uix_LITbasis.py> ')
    exit()
os.chdir(litpath3He + 'lit_bas_rhs/')

for streukanal in streukas:
    Jstreu = float(streukanal.split('^')[0])
    Jstreustring = '%s' % str(Jstreu)[:3]

    os.system('cp QUAOUT_J%s QUAOUT' % Jstreustring)
    os.system('cp DRQUAOUT_J%s DRQUAOUT' % Jstreustring)

    # -- load available basis states and diagonalize H in this `full` space
    bas = np.loadtxt(
        litpath3He + 'lit_bas_rhs/LITbas_full_J%s.dat' % Jstreustring).astype(
            int)

    basvs = {}
    for bv in bas:
        try:
            basvs[bv[0]].append(bv[1])
        except:
            basvs[bv[0]] = [bv[1]]
    sbas = sorted(basvs.items(), key=operator.itemgetter(0))

    n3_inen_bdg(sbas, Jstreu, costr, fn='INEN', pari=0, nzop=31, tni=11)

    os.system(BINBDGpath + 'DR2END_AK.exe && grep -A 3 \'EIGENWER\' OUTPUT')
    os.system('cp OUTPUT OUTPUT_full_J%s.dat' % Jstreustring)
    # --

    intwLIT = [
        np.array(ln.split(';')).astype(float).tolist()
        for ln in open(litpath3He + 'intw3heLIT_J%s.dat' % Jstreustring)
    ]

    relwLIT = [
        np.array(ln.split(';')).astype(float).tolist()
        for ln in open(litpath3He + 'relw3heLIT_J%s.dat' % Jstreustring)
    ]

    iws = []
    rws = []

    bass = []

    #print(' BV REL          wi          wr')
    th1_low = 0.0000003
    th2_low = 0.0000003
    th1_up = 0.09
    th2_up = 0.09
    w_max = max(th1_up, th2_up)
    # -- loop through the basis
    for bv in bas:
        for fr in range(len(intwLIT)):
            if bv[0] <= sum([len(ws) for ws in intwLIT[:fr + 1]]):
                iw = bv[0] - sum([len(ws) for ws in intwLIT[:fr]])
                iws.append(float(intwLIT[fr][iw - 1]))
                rws.append(float(relwLIT[fr][bv[1] - 1]))
                #print('%3d%3d%12.4f%12.4f' % (bv[0], bv[1], iws[-1], rws[-1]))
                # include bv in reduced basis if ``a criterium'' is satisfied

                if ((th1_low < iws[-1] < th1_up) &
                    (th2_low < rws[-1] < th2_up)):
                    bass.append([bv[0], bv[1]])

                #if (((th1_low < iws[-1] < th1_up) & (rws[-1] < th2_up)) |
                #    ((th2_low < rws[-1] < th2_up) & (iws[-1] < th1_up))):
                #    bass.append([bv[0], bv[1]])

                #if (((iws[-1] < th1_up) | (rws[-1] < th2_up)) &
                #    (w_max > iws[-1] > th1_low) & (w_max > rws[-1] > th2_low)):
                #    bass.append([bv[0], bv[1]])
                break
# --

    basvs = {}

    for bv in bass:
        try:
            basvs[bv[0]].append(bv[1])
        except:
            basvs[bv[0]] = [bv[1]]

    sbas = sorted(basvs.items(), key=operator.itemgetter(0))

    with open(litpath3He + 'lit_bas_rhs/LITbas_red_J%s.dat' % Jstreustring,
              'w') as f:
        np.savetxt(f, [[jj[0], kk] for jj in sbas for kk in jj[1]], fmt='%d')
    f.close()

    n3_inen_bdg(sbas, Jstreu, costr, fn='INEN', pari=0, nzop=31, tni=11)

    os.system(BINBDGpath + 'DR2END_AK.exe && grep -A 3 \'EIGENWER\' OUTPUT')
    os.system('cp OUTPUT OUTPUT_red_J%s.dat' % Jstreustring)

    EBDG = get_h_ev()[0]
    np.savetxt(
        litpath3He + '/Eth_J%s.dat' % Jstreustring,
        np.array([EBDG]),
        fmt='%12.4f')

    visbas(
        basispath=litpath3He + 'lit_bas_rhs/',
        widthpath=litpath3He,
        exepath=BINBDGpath,
        Jstrstr=Jstreustring)