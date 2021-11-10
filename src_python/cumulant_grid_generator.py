import matplotlib.pyplot as plt

import numpy as np
from scipy.optimize import minimize, fmin
from scipy import special


def widthEq(upperIntBound, centerArray, widthArray, yTic):

    LHS = 0.
    for n in range(len(centerArray)):
        LHS += special.erf(
            np.sqrt(widthArray[n]) * (upperIntBound - centerArray[n]))
    return (LHS - (2. * yTic - 1.) * len(centerArray))**2


def cumWidths(anza=10, centers=[1., 6.], widths=[1., 1.]):

    v0 = 2.1

    ww = []
    for yT in np.linspace(0.0, 1.0, anza + 1):
        newW = fmin(widthEq, v0, args=(centers, widths, yT), disp=False)[0]
        if newW >= 0.:
            ww.append(newW)
            v0 = ww[-1]

    return ww


anzW = 10
cw = cumWidths(anza=anzW)

binns = np.linspace(0., np.max(cw), anzW)

plt.hist(cw, bins=binns)
plt.show()