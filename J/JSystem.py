import numpy as np
from mpmath import mp
from ..core.AngularMomentum import J2, Jz, Jplus, Jminus  # , Jrange #, J2range
from ..core.Setup import resultmethod
from ..core.QuantumSystem import QuantumSystem, get_diag


class JSystem(QuantumSystem):

    def __init__(self, J, parent=None):

        QuantumSystem.__init__(self, parent)

        self.J = J
        self.J2 = J2(self.J)
        self.Jz = Jz(self.J)
        self.Jp = Jplus(self.J)
        self.Jm = Jminus(self.J)

        self.nstates = int(2*J+1)

    def _build(self):
        QuantumSystem._build(self)
        self._Js = get_diag(self._Xs.H*self.Jz*self._Xs)

    @property
    @resultmethod
    def Js(self):
        return np.real(np.array(self._Js.tolist(), dtype=complex).flatten())

    @resultmethod
    def J_transitions(self, power=1, N=None):
        def multiply(m1,m2):
            return mp.matrix(np.asarray(m1.tolist())*np.asarray(m2.tolist()))

        if power == 0:
            T = self.transitions(mp.matrix(np.eye(self.nstates)))[0]
            return multiply(T,T.conjugate())[:N, :N]
        
        JZ, JP, JM = self.transitions(self.Jz, self.Jp, self.Jm)

        if N is None or N > self.nstates:
            N = self.nstates

        return (((2*multiply(JZ,JZ.conjugate()) +
                multiply(JP,JP.conjugate()) +
                multiply(JM,JM.conjugate()))/2/self.J/(self.J+1))**power)[:N, :N]
        
