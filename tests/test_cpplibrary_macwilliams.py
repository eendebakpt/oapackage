import unittest

import numpy as np

import oapackage


class TestMacWilliams(unittest.TestCase):

    def test_macwilliams_transform(self):

        array = oapackage.exampleArray(1, 0)
        array = array.selectFirstColumns(3)
        arrayclass = oapackage.arraylink2arraydata(array)

        Dm = oapackage.distance_distribution_mixed(array, 0)
        D = np.array(Dm)

        N = array.n_rows
        s = 2
        factor_levels_for_groups = arrayclass.factor_levels_column_groups()

        Bp = oapackage.macwilliams_transform(D, N, s)
        Bpm = oapackage.macwilliams_transform_mixed(Dm, N, factor_levels_for_groups, verbose=0)

        np.testing.assert_array_equal(Bp, [1, 0, 0, .25])
        np.testing.assert_array_equal(Bp, Bpm)

        oapackage.GWLP(array)
        oapackage.GWLPmixed(array, verbose=2)


if __name__ == '__main__':
    """ Test code """
    unittest.main()
