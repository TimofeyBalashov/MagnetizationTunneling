import numpy as np
from mpmath import mp


def T(J):
    JJ = int(2*J+1)
    if J == int(J):
        M = mp.zeros(JJ)
        if J % 2 == 0:
            f = 1
        else:
            f = -1
    else:
        M = mp.zeros(JJ, dtype=complex)
        if (2*J-1)/2 % 2 == 0:
            f = 1j
        else:
            f = -1j
    for i in range(JJ):
        M[JJ-i-1, i] = f
        f *= -1

    return M


def timeInverse(Xs):
    # Xs is a matrix representing some set of eigenstates.
    # Every one will be time-reversed
    JJ = Xs.shape()[0]
    J = (JJ-1)/2.
    Tm = T(J)
    return np.conjugate(np.dot(Tm, Xs))
