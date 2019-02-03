# -*- coding: utf-8 -*-
"""
Pieter Eendebak <pieter.eendebak@gmail.com>

"""

#%%
import numpy as np
import unittest

import oapackage
import oapackage.conference

#%%


def _statistics_equal(statistics, expected):
    np.testing.assert_array_almost_equal(np.array(statistics), np.array(expected))

class TestJcharacteristics(unittest.TestCase):
    
    def test_Jcharacteristics(self):
        array = oapackage.exampleArray(46, 0)
        
        JJ = oapackage.Jcharacteristics_conference(array, number_of_columns=3)
        self.assertEqual( len(JJ), oapackage.choose(array.n_columns, 3))

        JJ = oapackage.Jcharacteristics_conference(array, 4)        
        self.assertEqual( len(JJ), oapackage.choose(array.n_columns, 4))
        hist=np.histogram(np.abs(JJ), bins=np.arange(0,21,2))[0]
        self.assertEqual(list(hist), [ 9,  0, 52,  0,  8,  0,  1,  0,  0,  0])

class TestConferenceStatistics(unittest.TestCase):

    def test_momentMatrix(self):
        M1 = np.array([[1., 0., 1. / 3],
                       [0., 1. / 3, 0.],
                       [1. / 3, 0., 0.2]])
        M = oapackage.conference.momentMatrix(1)
        np.testing.assert_array_equal(M, M1)
        M = oapackage.conference.momentMatrix(4)
        self.assertEqual(M.shape, (15, 15))

    def test_conference_model_statistics(self):
        arrays = [oapackage.exampleArray(idx, 0).selectFirstColumns(4) for idx in [46, 47, 48]]
        dsds = [oapackage.conference2DSD(array) for array in arrays]
        statistics = [oapackage.conference.modelStatistics(al, verbose=2) for al in dsds]
        expected = [(0, 0, 0), (0, 0, 0), (1, 17.081221494279234, 1.7610790263409104)]
        _statistics_equal(statistics, expected)

        arrays = [oapackage.exampleArray(idx, 0).selectLastColumns(4) for idx in [46, 47, 48]]
        dsds = [oapackage.conference2DSD(array) for array in arrays]
        statistics = [oapackage.conference.modelStatistics(al, verbose=0) for al in dsds]
        expected = [(1, 17.312358463836613, 1.8538666749267378), (1, 17.312358463836613,
                                                                  1.8538666749267378), (1, 17.08122149427922, 1.7610790263409093)]
        _statistics_equal(statistics, expected)

    def test_conference_projection_statistics(self):
        arrays = [oapackage.exampleArray(idx, 0) for idx in [45, 46, 47, 48]]
        statistics = [oapackage.conference.conferenceProjectionStatistics(array, verbose=1) for array in arrays]
        expected = [(0.9857142857142858, 16.831420420862223, 1.7523727421346018), (0.9857142857142858, 16.782711360159983, 1.7515927378529395),
                    (0.9714285714285714, 16.594006884359878, 1.7298655460036123), (1.0, 16.80031370680457, 1.7750342960753236)]
        _statistics_equal(statistics, expected)

    def test_leftDivide(self):
        A = np.array([[1, 2], [3, 4]])
        B = np.array([[51, -2], [3, 4]])
        C = oapackage.conference._leftDivide(A, B)
        np.testing.assert_array_almost_equal(C, np.array([[-99.,   8.], [75.,  -5.]]))


if __name__ == '__main__':
    unittest.main()
