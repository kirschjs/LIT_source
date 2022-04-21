import subprocess
import numpy as np
from bridgeA3 import *

anzStrBas = 21
StrBas = np.arange(1, 1 + anzStrBas)

if os.path.isdir(bkpdir) == False:
    subprocess.check_call(['mkdir', '-p', bkpdir])

for StreubasNR in StrBas:

    subprocess.call('./prun.sh %d' % StreubasNR, shell=True)
#    subprocess.call('python3 NextToNewestGeneration.py %d %d' %
#                    (StreubasNR, StreubasNR),
#                    shell=True)
#    subprocess.call('python3 A3_lit_M.py %d %d' % (StreubasNR, StreubasNR),
#                    shell=True)