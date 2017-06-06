import numpy as np
from mpmath import mp
from ..core.Setup import setupmethod
from ..core.MagneticField import ZeemanTerm as ZeemanTermCore


class ZeemanTerm(ZeemanTermCore):
    """ ZeemanTerm is a term of the form $\vec{B}\vec{J}$.

    The basis is assumed to be $\left|J_z\right>$, so the matrix
    ZT.B is a square matrix of size 2*ZT.J+1.

    The x, y and z components of the magnetic field can be accessed at any time
    with ZT.Bx, ZT.By and ZT.Bz. The magnitude and the direction of the field
    is available as ZT.Br, ZT.Bth and ZT.Bph

    The field components can be set either in Carthesian coordinates
    with setBx(), setBy(), setBz() and setBxyz() or in spherical with
    setBr(), setBtheta(), setBphi() and setBrtp().

    By default, the magnetic field units are meV. This can be changed
    by using setBFactor() to set an arbitrary conversion factor or
    setg() to set the Lande g-factor and convert magnetic field to Tesla.
    """
    @setupmethod
    def setg(self, *args):
        """The Zeeman term is $g\mju_B\vec{B}\vec{J}$,
           where g is the Lande g-factor.

           You can either set the g-factor directly with setg(g)
           or calculate g from L ans S with setg(L, S).
        """
        if len(args) == 1:  # Just g
            ZeemanTermCore.setg(self, args[0])
        elif len(args) == 2:  # L & S
            L, S = args
            J = L + S
            self.setg((3.002319*J*(J+1) +
                       1.002319*S*(S+1) -
                       1.002319*+L*(L+1))
                      / 2 / J / (J+1))
        else:
            raise ValueError("Wrong number of arguments to ZT.setg")

    def _build(self):
        '''Construct the operator matrix. If the y component of the field is
           zero, the matrix will be real, otherwise, it will be complex.'''

        if self.By == 0:
            self.B = self.parent.Jz*self.Bz \
                + (self.parent.Jp + self.parent.Jm)*self.Bx/2
        else:
            Bplus = mp.mpc(complex(self.Bx, self.By))
            Bminus = mp.mpc(complex(self.Bx, -self.By))
            self.B = self.parent.Jz*self.Bz \
                + (Bplus*self.parent.Jp + Bminus*self.parent.Jm)/2
        self.B *= self.conv
