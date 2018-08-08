""" Orthogonal Array package

"""

import numpy as np
import oapackage
import oalib
import oapackage.Doptim
import oapackage.graphtools
import unittest


def autodoctest():
    """ Test the module using autodoc
    Example:
      >>> import oapackage
      >>> arrayclass=oapackage.arraydata_t(2, 40, 0, 7)
      >>> print(arrayclass)
      arrayclass: N 40, k 7, strength 0, s {2,2,2,2,2,2,2}, order 0

    """
    return


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


class TestOahelper(unittest.TestCase):

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


class TestDoptimize(unittest.TestCase):

    def setUp(self):

        self.arrayclass = oapackage.arraydata_t(2, 16, 2, 6)
        self.dds = np.random.rand(20, 3)

    def test_unittest(self):
        scores, dds, sols, n = oapackage.Doptim.Doptimize(self.arrayclass, nrestarts=10, optimfunc=[
                                                          1, 0, 0], verbose=0, maxtime=18, selectpareto=False, nout=None, method=oalib.DOPTIM_UPDATE, niter=1000, nabort=0, dverbose=0)

        r = oapackage.Doptim.selectDn(scores, dds, sols, nout=1, sortfull=True)

    def test_optimDeffPython(self):
        al = oapackage.exampleArray(2)
        r, al = oapackage.Doptim.optimDeffPython(
            al, arrayclass=None, niter=10000, nabort=2500, verbose=0, alpha=[1, 0, 0], method=0)

    def test_generateDscatter(self):
        r = oapackage.Doptim.generateDscatter(self.dds, si=0, fi=1, lbls=None,
                                              ndata=3, nofig=True, fig=20, scatterarea=80)
        # plt.close(r['ax'])


class TestCppLibrary(unittest.TestCase):

    def test_miscunittest(self):
        miscunittest()


if __name__ == '__main__':
    """ Test code """
    unittest.main()
    # t=TestDoptimize()
