""" Orthogonal Array package test functions
"""

import sys
import os
import numpy as np
import tempfile
import unittest
if sys.version_info >= (3, 4):
    import unittest.mock as mock
    from unittest.mock import patch
else:
    import mock
    from mock import patch

import oalib
import oapackage
import oapackage.scanf
import oapackage.Doptim
import oapackage.graphtools


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
        assert(len(v) == G.shape[0])

    def test_exampleArray(self):
        # test a selection of the example arrays
        al = oapackage.exampleArray(5)
        assert(al.md5() == '3885c883d3bee0c7546511255bb5c3ae')
        al = oapackage.exampleArray(20)
        assert(np.array(al).shape == (24, 3))
        al = oapackage.exampleArray(40)
        assert(np.array(al).shape == (14, 5))

    def test_scanf(self):
        r = oapackage.scanf.sscanf('1', '%d')
        assert(r[0] == 1)

    def test_oa2graph(self):
        al = oapackage.exampleArray(2, 0)
        adata = oapackage.arraylink2arraydata(al)
        g = oapackage.graphtools.oa2graph(al, adata)
        assert(g[0].shape == (34, 34))


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

    try:
        # not possible right now...
        if verbose:
            print('direct float')
        A = np.eye(3).astype(float)
        al = oapackage.array_link(A)
        if verbose:
            al.showarray()
    except:
        pass


def test_nauty(verbose=0):
    if verbose:
        print('test_nauty: test reduction to normal form')
    al = oapackage.exampleArray(1, verbose)
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
    al = oalib.exampleArray(5, 1)
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

    at = oalib.reductionDOP(al)
    check = at.apply(al) == al.reduceDOP()
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

#%%


class TestOAfiles(unittest.TestCase):
    """ Test functionality related to orthogonal array files """

    def test_misc_file_operations(self):
        array_filename = tempfile.mktemp(suffix='.oa')
        lst = [oapackage.exampleArray(4, 1)]
        oapackage.writearrayfile(array_filename, lst)
        assert(oapackage.oahelper.oaIsBinary(array_filename) is False)
        oapackage.writearrayfile(array_filename, oapackage.exampleArray(4, 1), oapackage.ABINARY)
        assert(oapackage.oahelper.oaIsBinary(array_filename))

        oapackage.oahelper.oainfo(array_filename)

        _ = oapackage.oahelper.compressOAfile(array_filename)

    def test_findfilesR(self):
        _ = oapackage.oahelper.findfilesR(tempfile.tempdir, '.*oa')

    def test_checkArrayFile(self):
        a = tempfile.mktemp(suffix='.oa')
        self.assertFalse(oapackage.oahelper.checkArrayFile(a))
        self.assertTrue(oapackage.oahelper.checkArrayFile(a, -1))

    def test_finddirectories(self):
        _ = oapackage.oahelper.finddirectories(os.getcwd())


class TestOAhelper(unittest.TestCase):
    """ Test functionality contained in oahelper module """

    # def test_tilefigs(self):
    #   oapackage.oahelper.tilefigs([], geometry=[2,2])

    def setUp(self):
        self.test_array = oapackage.exampleArray(1, 1)

    def test_array2latex(self):

        latex_str = oapackage.oahelper.array2latex(np.array(self.test_array))
        self.assertEqual(latex_str[0:15], r'\begin{tabular}')
        latex_str = oapackage.oahelper.array2latex(np.array(self.test_array), mode='psmallmatrix')
        self.assertEqual(latex_str[0:15], r'\begin{psmallma')
        latex_str = oapackage.oahelper.array2latex(np.array(self.test_array), mode='pmatrix')
        self.assertEqual(latex_str[0:15], r'\begin{pmatrix}')

    def test_gwlp2str(self):
        self.assertEqual(oapackage.oahelper.gwlp2str([1, 2, 3]), '')
        self.assertEqual(oapackage.oahelper.gwlp2str([1, 0, .1]), '1.00,0.00,0.10')

    def test_argsort(self):
        idx = oapackage.oahelper.argsort([1, 2, 3])
        assert(idx == [0, 1, 2])
        idx = oapackage.oahelper.argsort([2, 2, 1])
        assert(idx == [2, 0, 1])

    def test_plot2Dline(self):
        with mock.patch('matplotlib.pyplot.plot') as MockPlt:
            _ = oapackage.oahelper.plot2Dline([1, 0, 0])
            self.assertTrue(MockPlt.called)

    def test_deprecated(self):
        def func():
            return 'hi'
        deprecated_function = oapackage.oahelper.deprecated(func)
        with self.assertWarns(Warning):
            _ = deprecated_function()

    def test_formatC(self):
        c_code = oapackage.oahelper.formatC(self.test_array)
        self.assertEqual(
            c_code, '\tarray_link al ( 16,5, 0 );\n\tint tmp[] = {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,1,0,1,1,1,0,1,1,1,0,0,0,1,0,0,1,0,1,0,1,1,1,0,1,1,0,0,1,0,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0};')

    def test_runExtend(self):
        N = 24
        k = 5
        t = 2
        r = oapackage.oahelper.runExtend(N, 4, t=t, l=2, verbose=3)
        rx = oapackage.oahelper.runExtend(N, k, t=t, l=2, initsols=r, verbose=3)
        self.assertTrue(len(r) == 10)
        self.assertTrue(len(rx) == 63)

    def test_joinArrayLists(self):
        l1 = [oapackage.exampleArray(2)]
        l2 = [oapackage.exampleArray(2), oapackage.exampleArray(2)]
        l = oapackage.oahelper.joinArrayLists([l1, l2])
        assert(len(l) == len(l1) + len(l2))

    def test_checkFiles(self):
        self.assertTrue(oapackage.oahelper.checkFiles([], -1))
        self.assertTrue(oapackage.oahelper.checkFiles(['nosuchfile'], -1))

    def test_checkFilesOA(self):
        r = oapackage.oahelper.checkFilesOA([], cache=1, verbose=0)
        self.assertTrue(r)

    def test_floatformat(self):
        s = oapackage.oahelper.floatformat(3.14, mind=2, maxd=4)
        self.assertEqual(s, '3.14')
        s = oapackage.oahelper.floatformat(3.14, mind=1, maxd=2)
        self.assertEqual(s, '3.1')

    def test_safemin(self):
        r = oapackage.oahelper.safemin(np.array([1, -2, 3]), default=0)
        self.assertEqual(r, -2)
        r = oapackage.oahelper.safemin(np.array([]), default=-2)
        self.assertEqual(r, -2)

    def test_choose(self):
        test_cases = [((3, 2), 3), ((10, 1), 10), ((5, 2), 10), ((1, 0), 1), ((-1, 0), 1)]
        for args, expected in test_cases:
            print(args)
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


class TestDoptimize(unittest.TestCase):

    def setUp(self):

        self.arrayclass = oapackage.arraydata_t(2, 16, 2, 6)
        self.dds = np.random.rand(20, 3)
        self.dds2 = np.array([[1, 1, 1], [1, 2, 1], [1, 2, 3], [2, 0, 1]])

        self.guitest = True
        try:
            import matplotlib.pyplot
        except:
            self.guitest = False
        print('guitest %s'  % self.guitest)
            
    def test_custom_optim(self):
        def optimfunc(x): return x[0] + x[1] + x[2]
        scores, dds, sols, n = oapackage.Doptim.Doptimize(self.arrayclass, nrestarts=2, optimfunc=optimfunc, verbose=1,
                                                          maxtime=18, selectpareto=False, nout=None, method=oalib.DOPTIM_UPDATE, niter=1000, nabort=0, dverbose=0)
        self.assertEqual(len(scores), n)
        self.assertEqual(len(dds), n)
        self.assertEqual(len(sols), n)
        scores, dds, sols, n = oapackage.Doptim.Doptimize(self.arrayclass, nrestarts=2, optimfunc=None, verbose=1,
                                                          maxtime=6, selectpareto=False, nout=None, method=oalib.DOPTIM_UPDATE, niter=30, nabort=0, dverbose=0)

    def test_unittest(self):
        scores, dds, sols, _ = oapackage.Doptim.Doptimize(self.arrayclass, nrestarts=10, optimfunc=[
                                                          1, 0, 0], verbose=1, maxtime=9, selectpareto=False, nout=None, method=oalib.DOPTIM_UPDATE, niter=1000, nabort=0, dverbose=0)

        result = oapackage.Doptim.selectDn(scores, dds, sols, nout=1, sortfull=True)

    def test_optimDeffPython(self):
        al = oapackage.exampleArray(2)
        _, al = oapackage.Doptim.optimDeffPython(
            al, arrayclass=None, niter=1000, nabort=1500, verbose=1, alpha=[1, 0, 0], method=0)

        for method in [oapackage.oalib.DOPTIM_SWAP, oapackage.oalib.DOPTIM_FLIP, oapackage.oalib.DOPTIM_UPDATE]:
            r, al = oapackage.Doptim.optimDeffPython(
                al, arrayclass=None, niter=100, nabort=200, verbose=0, alpha=[1, 0, 0], method=method)

    def test_generateDscatter(self):
        if self.guitest:
            fig = 100
        else:
            fig=None
        r = oapackage.Doptim.generateDscatter(self.dds, second_index=0, first_index=1, lbls=None, verbose=1,
                                              ndata=3, nofig=True, fig=fig, scatterarea=80)

    def test_generateDpage(self):
        outputdir = tempfile.mkdtemp()
        allarrays = [oapackage.exampleArray(2), oapackage.exampleArray(2)]
        dds = np.array([A.Defficiencies() for A in allarrays])
        arrayclass = oapackage.arraylink2arraydata(allarrays[0])
        page = oapackage.Doptim.generateDpage(outputdir, arrayclass, dds, allarrays,
                                              fig=None, optimfunc=[1, 0, 0], nofig=True)
        if self.guitest:
            print('test_generateDpage: run gui test')
            # page = oapackage.Doptim.generateDpage(outputdir, arrayclass, dds, allarrays,
            #                                  fig=100, optimfunc=[1, 0, 0], nofig=True)
            try:
                matplotlib.pyplot.close(100)
            except:
                pass

    def test_runcommand(self):
        oapackage.oahelper.runcommand('dir', dryrun=1, verbose=1)

    def test_filterPareto(self):
        dds = self.dds2
        scores = np.arange(dds.shape[0])
        sols = [None] * scores.size
        s, _, _ = oapackage.Doptim.filterPareto(scores, dds, sols, verbose=0)
        self.assertEqual(list(s), [2, 3])

    def test_calcScore(self):
        dds = np.random.rand(10, 3)
        scores = oapackage.Doptim.calcScore(dds, optimfunc=[1, 2, 3])
        assert(scores.shape == (dds.shape[0], ))

    def test_array2Dtable(self):
        sols = [oapackage.exampleArray(9, 0)]
        _ = oapackage.Doptim.array2Dtable(sols, verbose=1, titlestr=None)


if __name__ == '__main__':
    """ Test code """
    unittest.main()
    # t=TestDoptimize()
