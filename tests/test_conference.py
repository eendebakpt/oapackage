# -*- coding: utf-8 -*-
"""
Pieter Eendebak <pieter.eendebak@gmail.com>

"""

# %%
import sys
import unittest
import numpy as np

import oapackage
import oapackage.conference

if sys.version_info >= (3, 4):
    import io
    import unittest.mock as mock
    from unittest.mock import patch
else:
    io = None
    mock = None
    patch = None


def only_python3(function):
    """ Decorator to only execute a test function in Python 3 """
    python3 = sys.version_info >= (3, 4)
    if python3:
        def only_python3_function(*args, **kwargs):
            return function(*args, **kwargs)
    else:
        def only_python3_function(*args, **kwargs):
            return None
    return only_python3_function

# %%


def _statistics_equal(statistics, expected):
    np.testing.assert_array_almost_equal(np.array(statistics), np.array(expected))


class TestFvalues(unittest.TestCase):

    def test_FvaluesConference(self):
        array = oapackage.exampleArray(47, 0)
        js = oapackage.jstructconference_t(array)
        self.assertEqual(js.Jvalues(), (16, 12, 8, 4, 0))
        self.assertEqual(js.calculateF(), (0, 2, 4, 51, 13))

    def test_FvaluesConference_odd_runs(self):
        array = oapackage.exampleArray(47, 0)
        definitive_screening_design = oapackage.conference2DSD(array)
        with self.assertRaises(RuntimeError):
            _ = oapackage.jstructconference_t(definitive_screening_design)


class TestJcharacteristics(unittest.TestCase):

    def test_Jcharacteristics(self):
        array = oapackage.exampleArray(46, 0)

        JJ = oapackage.Jcharacteristics_conference(array, number_of_columns=3)
        self.assertEqual(len(JJ), oapackage.choose(array.n_columns, 3))

        JJ = oapackage.Jcharacteristics_conference(array, 4)
        self.assertEqual(len(JJ), oapackage.choose(array.n_columns, 4))
        hist = np.histogram(np.abs(JJ), bins=np.arange(0, 21, 2))[0]
        self.assertEqual(list(hist), [9, 0, 52, 0, 8, 0, 1, 0, 0, 0])


class TestConferenceStatistics(unittest.TestCase):

    def test_momentMatrix(self):
        M1 = np.array([[1., 0., 1. / 3],
                       [0., 1. / 3, 0.],
                       [1. / 3, 0., 0.2]])
        M = oapackage.conference.momentMatrix(1)
        np.testing.assert_array_equal(M, M1)
        M = oapackage.conference.momentMatrix(4)
        self.assertEqual(M.shape, (15, 15))

    @only_python3
    def test_conference_model_statistics(self):
        arrays = [oapackage.exampleArray(idx, 0).selectFirstColumns(4) for idx in [46, 47, 48]]
        dsds = [oapackage.conference2DSD(array) for array in arrays]
        with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            statistics = [oapackage.conference.modelStatistics(al, verbose=2) for al in dsds]
            stdout = mock_stdout.getvalue()
            self.assertIn('modelStatistics: condition number: ', stdout)

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
        statistics = [oapackage.conference.conferenceProjectionStatistics(array, verbose=0) for array in arrays]
        expected = [(0.9857142857142858, 16.831420420862223, 1.7523727421346018), (0.9857142857142858, 16.782711360159983, 1.7515927378529395),
                    (0.9714285714285714, 16.594006884359878, 1.7298655460036123), (1.0, 16.80031370680457, 1.7750342960753236)]
        _statistics_equal(statistics, expected)

    def test_leftDivide(self):
        A = np.array([[1, 2], [3, 4]])
        B = np.array([[51, -2], [3, 4]])
        C = oapackage.conference._leftDivide(A, B)
        np.testing.assert_array_almost_equal(C, np.array([[-99., 8.], [75., -5.]]))


if __name__ == '__main__':
    unittest.main()
