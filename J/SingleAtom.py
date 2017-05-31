# Numpy module
import numpy as np

from ..core.CrystalField import CrystalField
from ..core.Setup import buildmethod

from .MagneticField import ZeemanTerm
from .JSystem import JSystem


class SingleAtom(JSystem):
    def __init__(self, J, orbital, parent=None):

        JSystem.__init__(self, J, parent)

        self.CF = CrystalField(J, orb=orbital, parent=self)
        self.ZT = ZeemanTerm(self)

    @buildmethod
    def _buildH(self):
        self.CF.makeReady()
        self.ZT.makeReady()
        self._H = self.CF.CF + self.ZT.B

if __name__ == "__main__":

    def print_np_matrix(m):
        for xs in np.array(m):
            print(" ".join("{:8.2g}".format(xxs) for xxs in xs))

    import sys
    import re

    from ..core.AngularMomentum import Jrange

    j = 8
    if len(sys.argv) > 1:
        m = re.match('([0-9]+)(/2)?', sys.argv[1])
        if m is not None:
            j2 = int(m.group(1))
            if m.group(2) is None:
                j = j2
            else:
                j = j2/2.0
    sa = SingleAtom(j, 3)
    sa.CF.setSymmetry('C3v')
    try:
        sa.CF.setCoefficient(2, 0, -1)
        sa.CF.setCoefficient(4, 3, -0.001)
#        sa.CF.setCoefficient(6, 3, -0.01)
        sa.CF.setCoefficient(6, 6, -0.00001)
    except Exception:
        pass
    sa.ZT.setBz(0.001)
    sa.spectrum()
    print()
    print_np_matrix(sa.CF.CF)
    print()
    print_np_matrix(sa.Xs)

    print(sa.Xs.H[1, :]*np.diag(Jrange(sa.J))*sa.Xs[:, 0])
