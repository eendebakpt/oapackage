""" Orthogonal Array package test functions
"""

import numpy as np
import tempfile
import unittest

import oalib
import oapackage
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


def test_reduceGraphNauty():
    al = oapackage.exampleArray(3)
    v = oapackage.reduceGraphNauty(al)


def test_exampleArray():
    # test a selection of the example arrays
    al = oapackage.exampleArray(5)
    assert(al.md5() == '3885c883d3bee0c7546511255bb5c3ae')
    al = oapackage.exampleArray(20)
    assert(np.array(al).shape == (24, 3))
    al = oapackage.exampleArray(40)
    assert(np.array(al).shape == (14, 5))


def test_scanf():
    import oapackage
    r = oapackage.scanf.sscanf('1', '%d')
    assert(r[0] == 1)


def test_oa2graph():
    al = oapackage.exampleArray(2, 0)
    adata = oapackage.arraylink2arraydata(al)
    g = oapackage.graphtools.oa2graph(al, adata)
    assert(g[0].shape == (34, 34))

    # t = graph2arrayTransformation(pp, arrayclass, verbose=0):


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

    if 0:
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
        a = tempfile.mktemp(suffix='.oa')
        lst = [oapackage.exampleArray(4, 1)]
        oapackage.writearrayfile(a, lst)
        assert(oapackage.oahelper.oaIsBinary(a) is False)
        oapackage.writearrayfile(a, oapackage.exampleArray(4, 1), oapackage.ABINARY)
        assert(oapackage.oahelper.oaIsBinary(a))

        oapackage.oahelper.oainfo(a)

    def test_findfilesR(self):
        _ = oapackage.oahelper.findfilesR(tempfile.tempdir, '.*oa')

    def test_checkArrayFile(self):
        a = tempfile.mktemp(suffix='.oa')
        self.assertFalse(oapackage.oahelper.checkArrayFile(a))


class TestOAhelper(unittest.TestCase):
    """ Test functionality contained in oahelper module """

    # def test_tilefigs(self):
    #   oapackage.oahelper.tilefigs([], geometry=[2,2])

    def test_argsort():
        idx = oapackage.oahelper.argsort([1, 2, 3])
        assert(idx == [0, 1, 2])
        idx = oapackage.oahelper.argsort([2, 2, 1])
        assert(idx == [2, 0, 1])

    def test_runExtend(self):
        N = 24
        k = 5
        t = 2
        l = 2
        r = oapackage.oahelper.runExtend(N, 4, t=t, l=2, verbose=3)
        rx = oapackage.oahelper.runExtend(N, k, t=t, l=2, initsols=r, verbose=3)
        self.assertTrue(len(r) == 10)
        self.assertTrue(len(rx) == 63)

    def test_joinArrayLists(self):
        l1 = [oapackage.exampleArray(2)]
        l2 = [oapackage.exampleArray(2), oapackage.exampleArray(2)]
        l = oapackage.oahelper.joinArrayLists([l1, l2])
        assert(len(l) == len(l1) + len(l2))

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

    def test_create_pareto_element(self):
        values = [1, 2, 3]
        p = oapackage.oahelper.create_pareto_element(values, pareto=None)

    def test_designStandardError(self):
        al = oapackage.exampleArray(14, 0)
        v = oapackage.oahelper.designStandardError(al)

    def test_fac(self):
        self.assertEqual(oapackage.oahelper.fac(4), 24)

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

    def test_custom_optim(self):
        def optimfunc(x): return x[0] + x[1] + x[2]
        scores, dds, sols, n = oapackage.Doptim.Doptimize(self.arrayclass, nrestarts=2, optimfunc=optimfunc, verbose=1,
                                                          maxtime=18, selectpareto=False, nout=None, method=oalib.DOPTIM_UPDATE, niter=1000, nabort=0, dverbose=0)

        scores, dds, sols, n = oapackage.Doptim.Doptimize(self.arrayclass, nrestarts=2, optimfunc=None, verbose=1,
                                                          maxtime=6, selectpareto=False, nout=None, method=oalib.DOPTIM_UPDATE, niter=30, nabort=0, dverbose=0)

    def test_unittest(self):
        scores, dds, sols, n = oapackage.Doptim.Doptimize(self.arrayclass, nrestarts=10, optimfunc=[
                                                          1, 0, 0], verbose=1, maxtime=9, selectpareto=False, nout=None, method=oalib.DOPTIM_UPDATE, niter=1000, nabort=0, dverbose=0)

        r = oapackage.Doptim.selectDn(scores, dds, sols, nout=1, sortfull=True)

    def test_optimDeffPython(self):
        al = oapackage.exampleArray(2)
        r, al = oapackage.Doptim.optimDeffPython(
            al, arrayclass=None, niter=1000, nabort=1500, verbose=1, alpha=[1, 0, 0], method=0)

        for method in [oapackage.oalib.DOPTIM_SWAP, oapackage.oalib.DOPTIM_FLIP, oapackage.oalib.DOPTIM_UPDATE]:
            r, al = oapackage.Doptim.optimDeffPython(
                al, arrayclass=None, niter=100, nabort=200, verbose=0, alpha=[1, 0, 0], method=method)

    def test_generateDscatter(self):
        r = oapackage.Doptim.generateDscatter(self.dds, si=0, fi=1, lbls=None,
                                              ndata=3, nofig=True, fig=20, scatterarea=80)

    def test_generateDpage(self):
        outputdir = tempfile.mkdtemp()
        allarrays = [oapackage.exampleArray(2), oapackage.exampleArray(2)]
        dds = np.array([A.Defficiencies() for A in allarrays])
        arrayclass = oapackage.arraylink2arraydata(allarrays[0])
        p = oapackage.Doptim.generateDpage(outputdir, arrayclass, dds, allarrays,
                                           fig=None, optimfunc=[1, 0, 0], nofig=True)

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
        t = oapackage.Doptim.array2Dtable(sols, verbose=1, titlestr=None)


class TestCppLibrary(unittest.TestCase):

    def test_miscunittest(self):
        miscunittest()


if __name__ == '__main__':
    """ Test code """
    unittest.main()
    # t=TestDoptimize()
