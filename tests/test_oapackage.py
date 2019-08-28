""" Orthogonal Array package test functions
"""

import sys
import os
import numpy as np
import tempfile
import logging
import importlib
import unittest
if sys.version_info >= (3, 4):
    import unittest.mock as mock
    from unittest.mock import patch
    python3 = True
else:
    try:
        import mock
        from mock import patch
    except ImportError as ex:
        logging.exception(ex)
        raise Exception('to perform tests with python2 install the mock package (see https://pypi.org/project/mock/)')
    python3 = False

import oalib
import oapackage
import oapackage.scanf
import oapackage.graphtools


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


def autodoctest():
    """ Test the module using autodoc
    Example:
      >>> import oapackage
      >>> arrayclass=oapackage.arraydata_t(2, 40, 0, 7)
      >>> print(arrayclass)
      arrayclass: N 40, k 7, strength 0, s {2,2,2,2,2,2,2}, order 0

    """
    return


class TestMisc(unittest.TestCase):

    def test_reduceGraphNauty(self):
        G = np.zeros((5, 5), dtype=int)
        G[1, 0] = G[0, 1] = 1
        v = oapackage.reduceGraphNauty(G)
        self.assertTrue(len(v) == G.shape[0])

    def test_exampleArray(self):
        # test a selection of the example arrays
        al = oapackage.exampleArray(5, 0)
        self.assertTrue(al.md5() == '3885c883d3bee0c7546511255bb5c3ae')
        al = oapackage.exampleArray(20, 0)
        self.assertTrue(np.array(al).shape == (24, 3))
        al = oapackage.exampleArray(40, 0)
        self.assertTrue(np.array(al).shape == (14, 5))

    def test_scanf(self):
        r = oapackage.scanf.sscanf('1', '%d')
        self.assertTrue(r[0] == 1)

    def test_oa2graph(self):
        al = oapackage.exampleArray(2, 0)
        adata = oapackage.arraylink2arraydata(al)
        g = oapackage.graphtools.oa2graph(al, adata)
        self.assertTrue(g[0].shape == (34, 34))


def test_numpy_interface(verbose=0):
    A = np.eye(3, 4).astype(int)
    A[0, :] = [10, 20, 30, 50]
    A[2, :] = [-3, -4, -5, -6]

    if verbose:
        print('makearraylink')
    al = oapackage.makearraylink(A)

    np.testing.assert_array_equal(np.array(al), A)

    if verbose:
        al.showarray()

    if verbose:
        print('direct')
    al = oapackage.array_link(A)
    if verbose:
        al.showarray()
    Ax = np.array(al)
    if verbose:
        print(A)
        print(Ax)

    with np.testing.assert_raises(TypeError):
        # not possible right now...
        if verbose:
            print('direct float')
        A = np.eye(3).astype(float)
        al = oapackage.array_link(A)
        if verbose:
            al.showarray()


def test_nauty(verbose=0):
    if verbose:
        print('test_nauty: test reduction to normal form')
    al = oapackage.exampleArray(0, verbose)
    alr = al.randomperm()
    tr = oapackage.reduceOAnauty(alr)
    alx = tr.apply(alr)
    assert(alx == al)


def miscunittest(verbose=1):
    """ Perform some unit testing, return True if succesfull """
    if verbose:
        print('oapackage: unittest: oalib version %s' % oalib.version())
    al = oalib.array_link()
    ii = 0
    al = oalib.exampleArray(ii, 0)

    arrayclass = oalib.arraydata_t(2, 16, 2, 6)  # fine
    alist = [arrayclass.create_root()]
    r = oalib.extend_arraylist(alist, arrayclass)
    if len(r) != 3:
        raise Exception('extension generation for %s failed' % (arrayclass, ))

    if not isinstance(al.getarray(), np.ndarray):
        print(
            'oapackage: unittest: error: array interface not working properly')
    else:
        if not al[2, 0] == al.getarray()[2, 0]:
            print(
                'oapackage: unittest: error: array interface not working properly')

    arrayclass = oalib.arraylink2arraydata(al)

    if verbose >= 2:
        print('unittest: calculate efficiencies')
    Deff = al.Defficiency()
    aa = oalib.Aefficiencies(al)
    assert(aa[0] == 1.0)
    assert(aa[1] == 1.0)
    assert(aa[2] == 1.0)

    if verbose >= 2:
        print('## oapackage test: example array %d: Deff %.3f' % (ii, Deff))

    # DOP reduction
    if verbose >= 2:
        print('unittest: test delete-one-factor GWLP reduction')
    al = oalib.exampleArray(5, verbose)
    al2 = al.randomperm()

    alr = al.reduceDOP()
    al2r = al2.reduceDOP()
    if not alr == al2r:
        print('error: DOP reduced arrays unequal!: %d' % (alr == al2r))
        print('alr')
        alr.showarraycompact()
        print('al2r')
        al2r.showarraycompact()
        return False

    transformation = oalib.reductionDOP(al)
    check = transformation.apply(al) == al.reduceDOP()
    if not check:
        print('error: DOP reduction transformation is invalid')

    # test graphtools
    if verbose >= 2:
        print('unittest: test graphtools')
    from oapackage.graphtools import oa2graph
    arrayclass = oalib.arraylink2arraydata(al)
    _ = oa2graph(al, arrayclass)

    test_numpy_interface()
    test_nauty()

    return True

# %%


class TestOAfiles(unittest.TestCase):
    """ Test functionality related to orthogonal array files """

    def test_misc_file_operations(self):
        _, array_filename = tempfile.mkstemp(suffix='.oa')
        lst = [oapackage.exampleArray(4, 0)]
        oapackage.writearrayfile(array_filename, lst)
        assert(oapackage.oahelper.oaIsBinary(array_filename) is False)
        oapackage.writearrayfile(array_filename, oapackage.exampleArray(4, 0), oapackage.ABINARY)
        assert(oapackage.oahelper.oaIsBinary(array_filename))

        oapackage.oahelper.oainfo(array_filename)

        _ = oapackage.oahelper.compressOAfile(array_filename)

    def test_selectArrays(self):
        _, array_filename = tempfile.mkstemp(suffix='.oa', dir=tempfile.tempdir)
        _, array_filename_out = tempfile.mkstemp(suffix='.oa', dir=tempfile.tempdir)
        oapackage.writearrayfile(array_filename, [oapackage.exampleArray(4, 0), oapackage.exampleArray(4, 0)])
        oapackage.oahelper.selectArraysInFile(array_filename, array_filename_out, [
            1], afmode=oalib.ABINARY, verbose=1, cache=0)
        oapackage.oahelper.selectArraysInFile(array_filename, array_filename_out,
                                              [-1], afmode=oalib.ABINARY, verbose=1, cache=0)

    def test_nArrayFile(self):
        _, array_filename = tempfile.mkstemp(suffix='.oa', dir=tempfile.tempdir)
        oapackage.writearrayfile(array_filename, [oapackage.exampleArray(4, 0)])
        self.assertEqual(oapackage.oahelper.nArrayFile(array_filename), 1)
        self.assertEqual(oapackage.oahelper.nArrayFile('notavalidfile.oa'), -1)

    def test_findfiles(self):
        _, array_filename = tempfile.mkstemp(suffix='.oa', dir=tempfile.tempdir)
        lst = [oapackage.exampleArray(4, 0)]
        oapackage.writearrayfile(array_filename, lst)
        lst = oapackage.oahelper.findfiles(tempfile.tempdir, '.*oa')
        self.assertIn(os.path.split(array_filename)[-1], lst)

    def test_findfilesR(self):
        _ = oapackage.oahelper.findfilesR(tempfile.tempdir, '.*oa')

    def test_checkArrayFile(self):
        a = tempfile.mktemp(suffix='.oa')
        self.assertFalse(oapackage.oahelper.checkArrayFile(a))
        self.assertTrue(oapackage.oahelper.checkArrayFile(a, -1))

    def test_finddirectories(self):
        _ = oapackage.oahelper.finddirectories(os.getcwd())


class TestParetoFunctionality:

    def test_selectParetoArrays(self):

        arrays = [oapackage.array_link(np.array([[ii]])) for ii in range(5)]
        pareto_object = oapackage.ParetoLongLong()

        for ii in range(len(arrays)):
            value = [ii, ii % 2]
            pareto_object.addvalue(value, ii)
        pareto_object.show(2)

        selected = oapackage.oahelper.selectParetoArrays(arrays, pareto_object)
        self.assertEqual(selected, arrays[4, 5])


class TestOAhelper(unittest.TestCase):
    """ Test functionality contained in oahelper module """

    def setUp(self):
        self.test_array = oapackage.exampleArray(1, 0)

    def test_helmert_contrasts(self):
        hc = oapackage.oahelper.helmert_contrasts(2, verbose=0)
        np.testing.assert_array_almost_equal(hc, np.array([[-1.], [1.]]))

        hc = oapackage.oahelper.helmert_contrasts(3, verbose=0)
        np.testing.assert_array_almost_equal(hc, np.array([[-1.22474487, -0.70710678],
                                                           [1.22474487, -0.70710678],
                                                           [0., 1.41421356]]))

        hc = oapackage.oahelper.helmert_contrasts(10, verbose=0)
        np.testing.assert_array_almost_equal(hc[0], np.array([-2.23606798, -1.29099445, -0.91287093, -0.70710678, -0.57735027,
                                                              -0.48795004, -0.42257713, -0.372678, -0.33333333]))

        for num_levels in [4, 10, 12]:
            hc = oapackage.oahelper.helmert_contrasts(num_levels, verbose=0)
            np.testing.assert_array_almost_equal(hc.T.dot(hc), num_levels * np.eye(num_levels - 1))

    def test_array2latex(self):

        latex_str = oapackage.oahelper.array2latex(np.array(self.test_array))
        self.assertEqual(latex_str[0:15], r'\begin{tabular}')
        latex_str = oapackage.oahelper.array2latex(np.array(self.test_array), mode='psmallmatrix')
        self.assertEqual(latex_str[0:15], r'\begin{psmallma')
        latex_str = oapackage.oahelper.array2latex(np.array(self.test_array), mode='pmatrix')
        self.assertEqual(latex_str[0:15], r'\begin{pmatrix}')

    @only_python3
    def test_gwlp2str(self):
        with self.assertWarns(UserWarning):
            self.assertEqual(oapackage.oahelper.gwlp2str([1, 2, 3]), '')
        self.assertEqual(oapackage.oahelper.gwlp2str([1, 0, .1]), '1.00,0.00,0.10')

    def test_parseProcessingTime(self):
        filename = tempfile.mktemp(suffix='.txt')
        with open(filename, 'wt') as fid:
            fid.write('# my logfile\n')
            fid.write('#time start: 2012-01-19 17:21:00\n')
            fid.write('#time end: 2012-01-19 18:21:00\n')
        dtt = oapackage.oahelper.parseProcessingTime(filename, verbose=2)
        self.assertEqual(dtt, 3600.)

    def test_argsort(self):
        idx = oapackage.oahelper.argsort([1, 2, 3])
        assert(idx == [0, 1, 2])
        idx = oapackage.oahelper.argsort([2, 2, 1])
        assert(idx == [2, 0, 1])

    @only_python3
    def test_plot2Dline(self):
        if importlib.util.find_spec('matplotlib') is not None:
            with mock.patch('matplotlib.pyplot.plot') as MockPlt:
                oapackage.oahelper.plot2Dline([1, 0, 0])
                self.assertTrue(MockPlt.called)

    @only_python3
    def test_deprecated(self):
        def func():
            return 'hi'
        deprecated_function = oapackage.oahelper.deprecated(func)
        with self.assertWarns(Warning):
            _ = deprecated_function()

    def test_formatC(self):
        c_code = oapackage.oahelper.formatC(self.test_array)
        self.assertEqual(
            c_code, '\tarray_link array ( 16,5, 0 );\n\tint array_data_tmp[] = {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,1,0,1,1,1,0,1,1,1,0,0,0,1,0,0,1,0,1,0,1,1,1,0,1,1,0,0,1,0,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0};\tarray.setarraydata(array_data_tmp, array.n_rows * array.n_columns);\n')

    def test_runExtend(self):
        N = 24
        k = 5
        t = 2
        r = oapackage.oahelper.runExtend(N, 4, t=t, l=2, verbose=3)
        rx = oapackage.oahelper.runExtend(N, k, t=t, l=2, initsols=r, verbose=3)
        self.assertTrue(len(r) == 10)
        self.assertTrue(len(rx) == 63)

    def test_joinArrayLists(self):
        l1 = [oapackage.exampleArray(2, 0)]
        l2 = [oapackage.exampleArray(2, 0), oapackage.exampleArray(2, 0)]
        l = oapackage.oahelper.joinArrayLists([l1, l2])
        assert(len(l) == len(l1) + len(l2))

    def test_checkFiles(self):
        self.assertTrue(oapackage.oahelper.checkFiles([], -1))
        self.assertTrue(oapackage.oahelper.checkFiles(['nosuchfile'], -1))
        self.assertTrue(oapackage.oahelper.checkFiles('nosuchfile', -1))

    def test_checkFilesOA(self):
        r = oapackage.oahelper.checkFilesOA([], cache=1, verbose=0)
        self.assertTrue(r)

    def test_floatformat(self):
        s = oapackage.oahelper.floatformat(3.14, mind=2, maxd=4)
        self.assertEqual(s, '3.14')
        s = oapackage.oahelper.floatformat(3.14, mind=1, maxd=2)
        self.assertEqual(s, '3.1')

    def test_safemin_safemax(self):
        r = oapackage.oahelper.safemin(np.array([1, -2, 3]), default=0)
        self.assertEqual(r, -2)
        r = oapackage.oahelper.safemin(np.array([]), default=-2)
        self.assertEqual(r, -2)
        r = oapackage.oahelper.safemax(np.array([1, -2, 3]), default=0)
        self.assertEqual(r, 3)
        r = oapackage.oahelper.safemax([], default=-23)
        self.assertEqual(r, -23)

    def test_choose(self):
        test_cases = [((3, 2), 3), ((10, 1), 10), ((5, 2), 10), ((1, 0), 1), ((-1, 0), 1)]
        for args, expected in test_cases:
            result = oapackage.oahelper.choose(*args)
            self.assertEqual(result, expected)

    def test_create_pareto_element(self):
        values = [1, 2, 3]
        p = oapackage.oahelper.create_pareto_element(values, pareto=None)
        self.assertEqual(p, values)

        self.assertRaises(Exception, oapackage.oahelper.create_pareto_element, dict())

        pareto = oalib.ParetoMultiDoubleLong()
        pareto_element = oapackage.oahelper.create_pareto_element([[1, 2], [3, 4]], pareto)
        self.assertEqual(pareto_element.size(), 2)

        pareto = oalib.ParetoMultiLongLong()
        pareto_element = oapackage.oahelper.create_pareto_element([[1, 2], [3, 4]], pareto)
        self.assertEqual(pareto_element.size(), 2)

    def test_create_pareto_element_invalid_type(self):
        with self.assertRaises(Exception):
            oapackage.oahelper.create_pareto_element([1.], 1.)

    def test_designStandardError(self):
        al = oapackage.exampleArray(14, 0)
        v = oapackage.oahelper.designStandardError(al)
        self.assertAlmostEqual(v[0], 0.3747931073686535)
        al = oapackage.exampleArray(9, 0)
        v = oapackage.oahelper.designStandardError(al)
        np.testing.assert_array_almost_equal(v[1], np.array([0.1679305, 0.17229075, 0.17286095, 0.17287786, 0.17303912,
                                                             0.17353519, 0.17548291]))

    def test_fac(self):
        self.assertEqual(oapackage.oahelper.fac(4), 24)

    def test_testHtml(self):
        import webbrowser
        with patch.object(webbrowser, "open", return_value=None):
            oapackage.oahelper.testHtml('<p>hi</p>')

    def test_bounds(self):
        b = oapackage.oahelper.DefficiencyBound(.8, 4, 6)
        self.assertAlmostEqual(b, 0.8944271909999)

    def test_array2html(self):
        X = np.array([[1, 2], [3, 4]], dtype='U200')
        h = oapackage.oahelper.array2html(X, header=1, tablestyle='border-collapse: collapse;',
                                          trclass='', tdstyle='', trstyle='', thstyle='')
        assert('table' in str(h))

    def test_sortrows(self):
        a = np.array([[1, 1], [-2, 2], [-3, 3], [-4, 4], [5, 5]])
        s = oapackage.oahelper.sortrows(a)
        self.assertTrue(np.all(s == [3, 2, 1, 0, 4]))

        x = -np.array([[0, 0, 1], [0, 1, 0], [0, 1, 1]])
        idx = oapackage.oahelper.sortrows(x)
        assert(np.all(idx == np.array([2, 1, 0], dtype=np.int64)))

        x = np.array([])
        idx = oapackage.oahelper.sortrows(x)
        self.assertTrue(len(idx) == 0)

    def test_sortcols(self):
        a = np.array([[1, 1], [-2, 2], [-3, 3], [-4, 4], [5, 5]]).T
        s = oapackage.oahelper.sortcols(a)
        self.assertTrue(np.all(s == [3, 2, 1, 0, 4]))


if __name__ == '__main__':
    unittest.main()
