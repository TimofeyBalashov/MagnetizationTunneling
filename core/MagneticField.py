from .Setup import SetupClass, setupmethod
from math import sin, cos, atan2, sqrt


class ZeemanTerm(SetupClass):
    """ ZeemanTerm is a term of the form $\vec{B}\vec{J}$.

    The basis is assumed to be $\left|J_z\right>$, so the matrix
    ZT.B is a square matrix of size 2*ZT.J+1.

    The x, y and z components of the magnetic field can be accessed
    at any time with ZT.Bx, ZT.By and ZT.Bz. The magnitude and the
    direction of the field is available as ZT.Br, ZT.Bth and ZT.Bph

    The field components can be set either in Carthesian coordinates
    with setBx(), setBy(), setBz() and setBxyz() or in spherical with
    setBr(), setBtheta(), setBphi() and setBrtp().

    By default, the magnetic field units are meV. This can be changed
    by using setBFactor() to set an arbitrary conversion factor or
    setg() to set the Lande g-factor and convert magnetic field to Tesla.
    """
    def __init__(self, parent=None):

        SetupClass.__init__(self, parent)

        self.Bx = 0
        self.By = 0
        self.Bz = 0
        self.Br = 0
        self.Bth = 0
        self.Bph = 0

        self.conv = 1

    @setupmethod
    def setBFactor(self, f):
        """The Zeeman term is $f\vec{B}\vec{J}$"""
        self.conv = f

    @setupmethod
    def setg(self, g):
        """The Zeeman term is $g\mju_B\vec{B}\vec{J}$,
           where g is the Lande g-factor.
        """
        self.setBFactor(0.057883818066*g)

    @setupmethod
    def setBx(self, Bx):
        """Set the x component of B in carthesian coordinates to Bx.
                   y- and z- components are not affected."""
        if (Bx == self.Bx):
            return
        self.Bx = Bx
        self._setRTP()

    @setupmethod
    def setBy(self, By):
        """Set the y component of B in carthesian coordinates to By.
                   x- and z- components are not affected."""
        if (By == self.By):
            return
        self.By = By
        self._setRTP()

    @setupmethod
    def setBz(self, Bz):
        """Set the z component of B in carthesian coordinates to Bz.
                   x- and y- components are not affected."""
        if (Bz == self.Bz):
            return
        self.Bz = Bz
        self._setRTP()

    @setupmethod
    def setBxyz(self, Bx, By, Bz):
        """Set all the components of B in carthesian coordinates."""
        self.Bx = Bx
        self.By = By
        self.Bz = Bz
        self._setRTP()

    @setupmethod
    def setBr(self, Br):
        """Set the magnitude of B to Br.
                   The direction of the field is not affected."""
        if (Br == self.Br):
            return
        self.Br = Br
        self._setXYZ()

    @setupmethod
    def setBtheta(self, Bth):
        """Set the azimuthal angle (between B and the z axis) to Bth.
                   The magnitude and the polar angle are not affected."""
        if (Bth == self.Bth):
            return
        self.Bth = Bth
        self._setXYZ()

    @setupmethod
    def setBphi(self, Bph):
        """Set the polar angle (between B and the x axis) to Bph.
                   The magnitude and the azimuthal angle are not affected."""
        if (Bph == self.Bph):
            return
        self.Bph = Bph
        self._setXYZ()

    @setupmethod
    def setBrtp(self, Br, Btheta, Bphi):
        """Set all the components of B in spherical coordinates."""
        self.Br = Br
        self.Bth = Btheta
        self.Bph = Bphi
        self._setXYZ()

    def _setXYZ(self):
        '''Recalculate Bx, By and Bz using Br, Bth and Bph.'''
        self.Bx = self.Br * sin(self.Bth) * cos(self.Bph)
        self.By = self.Br * sin(self.Bth) * sin(self.Bph)
        self.Bz = self.Br * cos(self.Bth)

    def _setRTP(self):
        '''Recalculate Br, Bth and Bph using Bx, By and Bz.'''
        self.Br = sqrt(self.Bx**2 + self.By**2 + self.Bz**2)
        if self.Br == 0:
            self.Bth = 0
            self.Bph = 0
        else:
            self.Bth = atan2(sqrt(self.Bx**2 + self.By**2), self.Bz)
            if (self.By == 0):
                self.Bph = 0
            else:
                self.Bph = atan2(self.Bx, self.By)
