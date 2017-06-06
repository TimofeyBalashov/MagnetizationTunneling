import numpy as np
from mpmath import mp

# EVERYTHING WORKS ONLY FOR SQUARE MATRICES !!!


def shift_slice(s, n):
    if s.start is None:
        return slice(n, s.stop + n, s.step)
    else:
        return slice(s.start+n, s.stop + n, s.step)


def bloat_slice(s, factor):
    if s.step is None:
        if s.start is None:
            return slice(0, s.stop*factor, factor)
        else:
            return slice(s.start,
                         s.start + (s.stop - s.start) * factor,
                         factor)
    else:
        if s.start is None:
            return slice(0, s.stop*factor, s.step*factor)
        else:
            return slice(s.start,
                         s.start + (s.stop - s.start) * factor,
                         s.step*factor)


def _expand(op, n, shifts, minor):
    """
    op - matrix
    n x n - size of the result
    imap - list of lists for mapping
    """
    expanded = mp.zeros(n, n, dtype=op.dtype)

    for s in shifts:
        for j in range(minor):
            expanded[
                shift_slice(bloat_slice(slice(op.shape[0]), minor), s+j),
                shift_slice(bloat_slice(slice(op.shape[1]), minor), s+j)
                ] = op

    return expanded


def _collect(lists):
    for l in lists:
        for q in l:
            yield q


def _srange(shift, num):
    for n in range(num):
        yield shift + n


class SpaceMixer:
    def __init__(self):
        self._spaces = []
        self._maps = []
        self._tags = {}

    def addSubspace(self, size, tag=None):
        if tag is not None:
            if tag in self._tags:
                raise ValueError(
                    "Tag {} already in use for subspace #{}".
                    format(tag, self._tags[tag]))
            self._tags[tag] = len(self._spaces)
        self._spaces.append(size)
        self._maps = [self._getMap(ssi) for ssi in range(len(self._spaces))]

    def size(self):
        return int(np.prod(self._spaces))

    def _getMap(self, ssi):
        elems = int(np.prod(self._spaces[0:ssi]))
        jump = int(np.prod(self._spaces[ssi:]))
        return [i*jump for i in range(elems)]

    def convertFromSubspace(self, ssi, op):
        """ssi - sub_space_index or tag"""
        if ssi in self._tags:
            ssi = self._tags[ssi]
        return _expand(op,
                       self.size(),
                       self._maps[ssi],
                       int(np.prod(self._spaces[ssi+1:])))


if __name__ == '__main__':
    import unittest

    class testMixer(unittest.TestCase):
        def setUp(self):
            self.m1 = mp.matrix('1 2;3 4')
            self.x = SpaceMixer()
            self.x.addSubspace(2)
            self.x.addSubspace(2)
            self.y = SpaceMixer()
            self.y.addSubspace(2)
            self.y.addSubspace(3)
            self.y.addSubspace(4)

        def assertEqualMatrices(self, m1, m2, message=None):
            self.assertTrue((m1 == m2).all(), message)

        def test_expand(self):
            """ make sure the trivial expansion works"""
            self.assertEqualMatrices(self.m1, _expand(self.m1, 2, [0], 1))

        def test_collect(self):
            self.assertEqual([i for i in _collect([[1],
                                                   [2, 3, 4],
                                                   [5, 6],
                                                   [],
                                                   [7, 8]])],
                             [1, 2, 3, 4, 5, 6, 7, 8])

        def test_getmap(self):
            self.assertEqual(self.x._getMap(0), [0])
            self.assertEqual(self.x._getMap(1), [0, 2])

            self.assertEqual(self.y._getMap(1), [0, 12])

        def test_convert(self):
            self.assertEqualMatrices(
                self.x.convertFromSubspace(0, self.m1),
                mp.matrix('1 0 2 0; 0 1 0 2; 3 0 4 0; 0 3 0 4'),
                self.x.convertFromSubspace(0, self.m1))
            self.assertEqualMatrices(
                self.x.convertFromSubspace(1, self.m1),
                mp.matrix('1 2 0 0; 3 4 0 0; 0 0 1 2; 0 0 3 4'),
                "x 1")

    unittest.main()
