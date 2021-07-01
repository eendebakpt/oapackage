""" Orthogonal Array package test functions
"""

import unittest

import numpy as np

import oapackage
import oapackage.graphtools
import oapackage.scanf


class TestGraphtools(unittest.TestCase):
    """ Test functionality related to orthogonal array files """

    def test_selectIsomorphismClasses(self):
        ll = [oapackage.exampleArray(0), oapackage.exampleArray(0)]
        indices, mm = oapackage.graphtools.selectIsomorphismClasses(ll, verbose=0)
        assert(indices[0] == indices[1])
        self.assertEqual(len(mm), 2)
        assert(np.all(mm[0] == mm[1]))


if __name__ == '__main__':
    unittest.main()
