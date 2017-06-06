import numpy as np
from mpmath import mp
import math
from .Setup import SetupClass, resultmethod

kB = 1/11.604


#def sortedEig(A):
    #"""Calculates eigenvectors and eigenvalues of A, and
       #returns a tuple (energies, eigenvectors), with states
       #sorted by energy."""
    #E, X = mp.eigh(A)

    #si = np.argsort(E)
    #E = E[si]
    #X = X[:, si]

    #return E, X


#def normalizedEig(A):
    #"""Calculates eigenvectors and eigenvalues of A, and
       #returns a tuple (norm_energies, eigenvectors), with states
       #sorted by energy, and energies normalized to cover [0, 1] range."""
    #E, X = mp.eigh(A)

    #si = np.argsort(E)
    #emin = np.min(E)
    #emax = np.max(E)

    #E = (E[si]-emin)/(emax-emin)
    #X = X[:, si]

    #return E, X

def get_diag(M):
    imax = min(M.rows, M.cols)
    ans = mp.zeros(imax, 1)
    for i in range(imax):
        ans[i] = M[i,i]
    return ans

class QuantumSystem(SetupClass):
    """A QuantumSystem is a general representation of a Schroedinger
       equation with Hamiltonian QS.H, energies QS.Es and an eigenstate
       matrix QS.Xs. The number of states (and size of all the matrixes)
       can be accessed at QS.nstates.

       The expectation values for any operator can be obtained by calling
       QS.spectrum(). For transitions between various states
       see QS.transitions().

       This is an abstract class. At least the Hamiltonian construction
       self._buildH(), that sets the inner variable self._H
       has to be implemented in a subclass.
       """

    def __init__(self, parent=None):
        SetupClass.__init__(self, parent)
        self.nstates = 0

    @property
    def H(self):
        """Calculates and returns the Hamiltonian of the system"""
        if not self.ready:
            self._buildH()
        return self._H

    @property
    @resultmethod
    def Es(self):
        return np.real(np.array(self._Es.tolist(), dtype=complex).flatten())

    @property
    @resultmethod
    def Xs(self):
        return self._Xs

    def _build(self):
        """Calculates the energies and eigenstates of the system."""
        self._Es, self._Xs = mp.eigh(self.H)

    @resultmethod
    def spectrum(self, *ops, N=None):
        """QS.spectrum(op1, op2, op3 ..., N=None)
        Calculates expectation values in every eigenstate
        for every operator passed to the function .

        Returns a tuple of NumPy arrays with the expectation values.

        By default calculates the energies of the eigenstates.

        Pass an N to the function to obtain only the values
        for the lowest N states."""

        if N is None or N > self.nstates:
            N = self.nstates

        if ops is None or len(ops) == 0:
            ops = (self.H, )
        return tuple(get_diag(M)[:N] for M in self.transitions(*ops))

    def print_spectrum(self, *ops, N=None, prefix="", format='10.3g'):
        """QS.print_spectrum(op1, op2, op3 ..., N=None, prefix='', format='10.3g')
        Calculates expectation values for every operator
        passed to the function in every eigenstate.
        Prints a table of the expectation values for the lowest N states.

        Every line will be prefixed by prefix."""
        S = list(self.spectrum(*ops, N=N))
        N = len(S[0])

        format_line = "{{:{}}}".format(format)

        for i in range(N):
            print(prefix + " ".join(format_line.format(float(M[i])) for M in S))

    @resultmethod
    def transitions(self, *ops):
        """QS.transitions(op1, op2, op3 ...)
        Calculates transition matrix elements for every operator
        passed to the function between every pair of eigenstates.
        Returns a tuple of NumPy matrices with the expectation values.

        The initial states span the columns of the matrix, the final states
        span the rows."""
        return tuple(self.Xs.H*op*self.Xs for op in ops)

    @resultmethod
    def expectBolzmann(self, T, *ops):
        """Calculates the expectation values of arbitrary operators
        at a finite temperature *T*. Assumes Bolzmann distribution.
        """
        Es_corr = self.Es - self.Es[0]  # States ordered by energy
        Z = sum(math.exp(-e/kB/T) for e in Es_corr)
        return tuple(
            sum(v*math.exp(-e/kB/T) for v, e in zip(Vs, Es_corr)) / Z
            for Vs in self.spectrum(*ops))
