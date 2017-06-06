from . import StevensOperators as StOp
from mpmath import mp
from .Setup import SetupClass, setupmethod


class CrystalField(SetupClass):

    def __init__(self, J=0, orb=0, ssym='', coeff=[],
                 no_constant_term=False, parent=None):

        SetupClass.__init__(self, parent)

        self.no_constant_term = no_constant_term

        self.J = 0
        self._sz = int(2*self.J+1)
        self.orbital = 0
        self.symmetry_string = ''
        self.symmetry = (0, True)
        self.coeff = []
        self.orders = []
        self.xorders = []  # additional orders, not user-controlled

        self.setJ(J)
        self.setOrbital(orb)
        self.setSymmetry(ssym)
        self.setCoefficients(coeff)

    @setupmethod
    def setSymmetry(self, ssym):
        def _str2sym(ssym):
            def cn(n):
                return (n, False, False)

            def cnv(n):
                return (n, True)

            def dn(n):
                return (n, False, True)

            if ssym == "Ci":
                return cn(1)
            elif ssym == "Cs":
                return cnv(1)

            # First symbol - C, D or S (or T or O)
            # (can also be I, but we're not doing quasicrystals
            stype = ssym[0]
            # Second symbol - 1-2-3-4-6 or i or s (or d or h)
            n = int(ssym[1])
            # Third symbol - None, d, v or h
            if len(ssym) > 2:
                subtype = ssym[2]
            else:
                subtype = None

            # There are only three different symmetry types: Cn, Cnv and Dn
            # They are coded as (n, False, False), (n, True)
            # and (n, False, True)
            # (n, can exclude '-', can exclude some '-')

            if stype == 'D':
                if n % 2 == 0:
                    if subtype is None:
                        return cnv(n)
                    elif subtype == 'h':
                        return cnv(n)
                    elif subtype == 'd':
                        return cnv(2*n)
                else:
                    if subtype is None:
                        return dn(n)
                    elif subtype == 'h':
                        return cnv(2*n)
                    elif subtype == 'd':
                        return cn(n)

            elif stype == 'S':
                if n % 4 == 0:
                    return cn(n)
                elif n % 2 == 0:
                    return cn(n/2)

            elif stype == 'C':
                if subtype is None:
                    return cn(n)
                elif subtype == 'v':
                    return cnv(n)
                elif subtype == 'h':
                    if n % 2 == 0:
                        return cn(n)
                    else:
                        return cn(2*n)

            raise ValueError('Unsupported symmetry string "{}"'.format(ssym))

        if ssym == self.symmetry_string:
            return

        if len(ssym) < 2:
            raise ValueError('Symmetry string "{}" too short'.format(ssym))

        self.symmetry_string = ssym
        self.symmetry = _str2sym(ssym)

        self._rebuildOrders()

    @setupmethod
    def _rebuildOrders(self):

        if not len(self.symmetry_string):
            self.orders = []
            self.ops = []
            return

        n = self.symmetry[0]
        nmax = 1+min(int(self.J), self.orbital)

        oldorders = self.orders

        # compute CF orders

        self.orders = []
        self.xorders = []

        for nn in range(1, nmax, 1):
            self.orders += [(nn*2, 0)]
            if self.symmetry[1]:  # Cnv
                for qq in range(1, int(2.0*nn/n + 1)):
                    self.orders += [(nn*2, qq*n)]
                    self.xorders += [(nn*2, qq*n, [(1, (nn*2, -qq*n))]),
                                     (nn*2, -qq*n, 0, nn*2, qq*n)]
            elif self.symmetry[2]:  # Dn
                for qq in range(1, 1+int(2.0*nn/n), 1):
                    self.orders += [(nn*2, qq*n*(1-2*qq % 2))]
                    self.xorders += [(nn*2,
                                     -qq*n*(1 - 2 * qq % 2),
                                      0,
                                      nn*2,
                                      qq*n*(1 - 2 * qq % 2))]
            else:  # Cn
                for qq in range(1, 1+int(nn*2.0/n), 1):
                    self.orders += [(nn*2, qq*n), (nn*2, -qq*n)]

        # copy coefficients

        oldcoeffs = self.coeff

        self.coeff = [0]*len(self.orders)

        for (i, o) in enumerate(oldorders):
            if o in self.orders:
                self.coeff[self.orders.index(o)] = oldcoeffs[i]

        # make matrices

        self.ops = [
            StOp.O(self.J, nn, q, self.no_constant_term)
            for (nn, q) in self.orders]

    def rotate(self, angle):
        pass

    @setupmethod
    def setJ(self, J):
        if J == self.J:
            return

        self.J = J
        self._sz = int(2*J+1)

        self._rebuildOrders()

    @setupmethod
    def setOrbital(self, orb):

        if isinstance(orb, str):
            try:
                orb = 'spdf'.index(orb)
            except ValueError:
                raise ValueError('Unsupported orbital')
        else:
            try:
                orb = int(orb)
            except TypeError:
                raise ValueError('Unsupported orbital')

        if orb == self.orbital:
            return

        self.orbital = orb

        self._rebuildOrders()

    @setupmethod
    def setCoefficients(self, coeff):
        if coeff == self.coeff:
            return

        if len(coeff) != len(self.orders):
            raise ValueError("Wrong number of parameters for CF")

        self.coeff = coeff

    @setupmethod
    def setCoefficient(self, n, q, coeff):
        if (n, q) not in self.orders:
            raise ValueError('Stevens operator not corresponding to symmetry')
        self.coeff[self.orders.index((n, q))] = coeff

    def _build(self):
        # traceback.print_stack()
        if any(q < 0 for (n, q) in self.orders):
            self.CF = mp.zeros(self._sz)
        else:
            self.CF = mp.zeros(self._sz)

        for (i, op) in enumerate(self.ops):
            self.CF += op*self.coeff[i]

if __name__ == "__main__":
    cf = CrystalField()
    cf.setJ(1)
    cf.setOrbital('p')
    cf.setSymmetry('C3v')
