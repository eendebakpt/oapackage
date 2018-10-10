""" Orthogonal Array package test functions
"""

import sys
import os
import numpy as np
import numpy
import unittest
if sys.version_info >= (3, 4):
    import unittest.mock as mock
else:
    import mock

import oapackage



class TestCppLibrary(unittest.TestCase):

    def test_splits(self):
        self.assertEqual(oapackage.splitTag([10,12]), '10.12')
        self.assertEqual(oapackage.splitFile([]), '')
        self.assertEqual(oapackage.splitDir([1,2]), 'sp0-split-1' + os.path.sep + 'sp1-split-2' + os.path.sep)
        
    def test_exception_handling(self):
        with self.assertRaises(RuntimeError):
            oapackage.mycheck_handler("file","function", 10, 0, "hi")
        oapackage.mycheck_handler("file","function", 10, 1, "hi")       
    
        with self.assertRaises(RuntimeError):
            oapackage.throw_runtime_exception("dsfs")
    
        al=oapackage.oalib.exampleArray(18,1)
        with self.assertRaisesRegexp(RuntimeError, "array cannot have negative elements"):
            _ = oapackage.array2eigenModelMatrixMixed (al, 2);

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

    def test_projection_efficiencies(self):
        al = oapackage.exampleArray(11, 1)
        d = oapackage.projDeff(al, 3, 1)
        D = al.selectFirstColumns(3).Defficiency()
        assert(D == d[0])
        numpy.testing.assert_almost_equal(numpy.mean(d), 0.99064112542249538329031111061340197921)

        seq = oapackage.PECsequence(al)
        numpy.testing.assert_equal(seq, (1.0,) * len(seq))

    def test_arraylink_slicing(self):
        numpy_array=np.arange(0, 6*10).reshape( (6,10))
        
        al=oapackage.makearraylink(numpy_array)
        assert(al[0]==numpy_array.flatten()[0])
        assert(al[0,1]==numpy_array[0,1])
        assert(al[4,2]==numpy_array[4,2])   
        np.testing.assert_equal(al[0:4,1:5], np.array(al)[0:4, 1:5] )
        np.testing.assert_equal(al[0:1,0:10:2], np.array(al)[0:1,0:10:2] )
        np.testing.assert_equal(al[3,3::], np.array(al)[3:4,3::] )
        np.testing.assert_equal(al[2,:8:2], np.array(al)[2:3,:8:2] )
        np.testing.assert_equal(al[2:3,:8:2], np.array(al)[2:3,:8:2] )

        with self.assertRaises(IndexError):
            al[-1,1]


if __name__ == '__main__':
    """ Test code """
    unittest.main()
    # t=TestDoptimize()
