import subprocess
import numpy as np

anzStrBas = 3
StrBas = np.arange(1, 1 + anzStrBas)

for StreubasNR in StrBas:

    subprocess.call('python3 NextToNewestGeneration.py %d %d' %
                    (StreubasNR, StreubasNR),
                    shell=True)
    subprocess.call('python3 A3_lit_M.py %d %d' % (StreubasNR, StreubasNR),
                    shell=True)