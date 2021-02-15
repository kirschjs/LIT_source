import numpy as np
from numpy import linalg as LA


def smart_ev(matout, threshold=10**-7):

    dim = int(np.sqrt(len(matout) * 0.5))

    # read Norm and Hamilton matrices
    normat = np.reshape(np.array(matout[:dim**2]).astype(float), (dim, dim))
    hammat = np.reshape(np.array(matout[dim**2:]).astype(float), (dim, dim))

    # normalize the matrices with the Norm's diagonal
    normdiag = [normat[n, n] for n in range(dim)]
    umnorm = np.diag(1. / np.sqrt(normdiag))
    nm = np.dot(np.dot(np.transpose(umnorm), normat), umnorm)
    hm = np.dot(np.dot(np.transpose(umnorm), hammat), umnorm)

    # diagonalize normalized norm
    ew, ev = LA.eig(nm)
    idx = ew.argsort()[::-1]
    ew = [eww for eww in ew[idx]]

    # project onto subspace with ev > threshold
    ew = [eww for eww in ew if np.real(eww) > threshold]
    dimRed = len(ew)
    ev = ev[:, idx][:, :dimRed]

    # transormation matric for (H-E*N)PSI=0 such that N->id
    Omat = np.dot(ev, np.diag(1. / np.sqrt(ew)))

    # diagonalize the projected Hamiltonian
    Hgood = np.dot(np.dot(np.transpose(Omat), hm), Omat)
    ewGood, evGood = LA.eig(Hgood)

    idx = ewGood.argsort()[::-1]
    ewGood = [eww for eww in ewGood[idx]]
    evGood = evGood[:, idx]

    print('dim0 = %d\ndimR = %d' % (dim, dimRed))
    # return the ordered eigenvalues
    return ewGood