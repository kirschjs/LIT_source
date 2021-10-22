import numpy as np


def polyWidths(wmin=10**-2, wmax=10, nbrw=10, npoly=4):

    wds = np.flip([(2. / x)**npoly * (wmin - wmax) / (1 - 2**npoly) + wmax +
                   (wmax - wmin) / (-1. + 1. / 2**npoly)
                   for x in np.linspace(1, 2, nbrw + 2)])
    return wds[1:-1]


wmi = 0.001
wma = 20.0
nbr = 6
npol = 12

wgeo = np.geomspace(
    #lit_w_tmp = np.abs(np.linspace(
    start=wmi,
    stop=wma,
    num=nbr,
    endpoint=True,
    dtype=None)

wset = polyWidths(wmin=wmi, wmax=wma, nbrw=nbr, npoly=npol)
print(wset)
print(wgeo)