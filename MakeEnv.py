import os
import subprocess

home = os.getenv("HOME")
source_dirs = [
    home + '/kette_repo/ComptonLIT/source/src_elma_new',
    home + '/kette_repo/ComptonLIT/source/src_nucl',
    home + '/kette_repo/ComptonLIT/source/src_nucl/V18_PAR',
    home + '/kette_repo/ComptonLIT/source/src_nucl/UIX_PAR'
]

for srcd in source_dirs:
    os.chdir(srcd)
    subprocess.call('make all', shell=True)