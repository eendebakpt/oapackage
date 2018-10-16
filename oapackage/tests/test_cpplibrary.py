""" Orthogonal Array package test functions
"""

import sys
import os
import numpy as np
import numpy
import unittest
if sys.version_info >= (3, 4):
    import unittest.mock as mock
    import io
    from unittest.mock import patch
    python3 = True
else:
    import mock
    python3 = False


import oapackage


class TestCppLibrary(unittest.TestCase):

    def test_splits(self):
        self.assertEqual(oapackage.splitTag([10, 12]), '10.12')
        self.assertEqual(oapackage.splitFile([]), '')
        self.assertEqual(oapackage.splitDir([1, 2]), 'sp0-split-1' + os.path.sep + 'sp1-split-2' + os.path.sep)

    def test_array_transformation_t(self):
        at = oapackage.array_transformation_t()
        if python3:
            with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
                at.show()
                std_output = mock_stdout.getvalue()
                self.assertEqual(std_output, 'array transformation: no class defined\n')

        al = oapackage.exampleArray(0)
        arrayclass = oapackage.arraylink2arraydata(al)
        at = oapackage.array_transformation_t(arrayclass)
        at.setcolperm([1, 0])
        _ = at.show()
        self.assertEqual(at.colperm(), (1, 0))

    def test_arraylink2arraydata(self):
        #ll = [oapackage.exampleArray(0), oapackage.exampleArray(0)]
        al = oapackage.exampleArray(0)
        adata = oapackage.arraylink2arraydata(al)
        self.assertEqual(str(adata), r'arrayclass: N 8, k 2, strength 2, s {2,2}, order 0')
        al = oapackage.array_link(4, 4, 0)
        al.setconstant(0)
        adata = oapackage.arraylink2arraydata(al)
        al.setconstant(-1)
        with self.assertRaises(RuntimeError):
            _ = oapackage.arraylink2arraydata(al)

    def test_exception_handling(self):
        with self.assertRaises(RuntimeError):
            oapackage.mycheck_handler("file", "function", 10, 0, "hi")
        oapackage.mycheck_handler("file", "function", 10, 1, "hi")

        with self.assertRaises(RuntimeError):
            oapackage.throw_runtime_exception("dsfs")

        al = oapackage.oalib.exampleArray(18, 1)
        with self.assertRaisesRegex(RuntimeError, "array cannot have negative elements"):
            _ = oapackage.array2eigenModelMatrixMixed(al, 2)

    def test_selectFirstColumns(self):
        al = oapackage.exampleArray(41, 1)
        al = al.selectFirstColumns(3)
        assert(al.n_columns == 3)

        al = oapackage.exampleArray(1000, 1)

        with self.assertRaises(RuntimeError):
            al = al.selectFirstColumns(1)

    def test_mycheck_handler(self):
        oapackage.mycheck_handler('a', 'b', 1, 1, 'bla')
        with self.assertRaises(RuntimeError):
            oapackage.mycheck_handler('a', 'b', 1, 0, 'bla')

    def distance_distribution(self):
        al=oapackage.array_link(2,2, 0)
        distance_distrib = oapackage.distance_distribution(al)
        self.assertEqual(distance_distrib, (2.0, 0.0, 0.0))

        al[1,0]=1
        distance_distrib = oapackage.distance_distribution(al)
        self.assertEqual(distance_distrib, (1.0, 1.0, 0.0))

        al = oapackage.exampleArray(2, 1)

        distance_distrib = oapackage.distance_distribution(al)
        self.assertEqual(distance_distrib, (1.25, 0.75, 1.5, 6.5, 5.25, 0.75, 0.0))
        
    def test_projection_efficiencies(self):
        al = oapackage.exampleArray(11, 1)
        d = oapackage.projDeff(al, 3, 1)
        D = al.selectFirstColumns(3).Defficiency()
        assert(D == d[0])
        numpy.testing.assert_almost_equal(numpy.mean(d), 0.99064112542249538329031111061340197921)

        pec_seq = oapackage.PECsequence(al)
        numpy.testing.assert_equal(pec_seq, (1.0,) * len(pec_seq))
        pic_seq = oapackage.PICsequence(al)
        numpy.testing.assert_equal(pic_seq, (0.9985780064264659, 0.9965697009006985, 0.9906411254224957, 0.9797170906488152, 0.9635206782887167, 0.9421350381959234, 0.9162739059686846, 0.8879176205539139))

    def test_arraylink(self):
        al=oapackage.exampleArray(0)
        self.assertEqual(al.at(0,1), 0)
        self.assertEqual(al.at(0), 0)
        with self.assertRaises(IndexError):
            al.at(-1,0)
        with self.assertRaises(IndexError):
            al.at(-1)
        with self.assertRaises(IndexError):
            al.at(0,al.n_columns)
            
    def test_arraylink_slicing(self):
        numpy_array = np.arange(0, 6 * 10).reshape((6, 10))

        al = oapackage.makearraylink(numpy_array)
        assert(al[0] == numpy_array.flatten()[0])
        assert(al[0, 1] == numpy_array[0, 1])
        assert(al[4, 2] == numpy_array[4, 2])
        np.testing.assert_equal(al[0:4, 1:5], np.array(al)[0:4, 1:5])
        np.testing.assert_equal(al[0:1, 0:10:2], np.array(al)[0:1, 0:10:2])
        np.testing.assert_equal(al[3, 3::], np.array(al)[3:4, 3::])
        np.testing.assert_equal(al[2, :8:2], np.array(al)[2:3, :8:2])
        np.testing.assert_equal(al[2:3, :8:2], np.array(al)[2:3, :8:2])

        with self.assertRaises(IndexError):
            al[-1, 1]

    def test_conference_generation(self):
        
        al=oapackage.exampleArray(42)
        lst=[al]
        conference_type=oapackage.conference_t(al.n_rows, al.n_rows, 0)
        
        extensions = oapackage.extend_conference_restricted (lst, conference_type, verbose=1)
        self.assertEqual(len(extensions), 10)

        self.assertEqual(extensions[0].md5(), 'f759e75d3ce6adda5489fed4c528a6fb')
        self.assertEqual(extensions[-1].md5(), 'ec13ed02c70dbc9975a99e70e23154b3')


if __name__ == '__main__':
    """ Test code """
    unittest.main()
    # t=TestDoptimize()
