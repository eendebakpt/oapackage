""" Orthogonal Array package test functions
"""

import sys
import numpy as np
import tempfile
import logging
import unittest
if sys.version_info >= (3, 4):
    import unittest.mock as mock
    import io
    from unittest.mock import patch
    python3 = True
else:
    try:
        import mock
    except ImportError as ex:
        logging.exception(ex)
        raise Exception('to perform tests with python2 install the mock package (see https://pypi.org/project/mock/)')
    python3 = False
    patch = None

import oalib
import oapackage
import oapackage.scanf
import oapackage.Doptim


def only_python3(function):
    if python3:
        def only_python3_function(*args, **kwargs):
            return function(*args, **kwargs)
    else:
        def only_python3_function(*args, **kwargs):
            return None
    return only_python3_function


class TestDoptimize(unittest.TestCase):

    def setUp(self):

        self.arrayclass = oapackage.arraydata_t(2, 16, 0, 6)
        self.dds = np.random.rand(20, 3)
        self.dds2 = np.array([[1, 1, 1], [1, 2, 1], [1, 2, 3], [2, 0, 1]])

        self.guitest = True
        try:
            import matplotlib.pyplot
        except:
            self.guitest = False

    @only_python3
    def test_custom_optim(self):
        def optimfunc(x): return x[0] + x[1] + x[2]
        with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            scores, dds, sols, n = oapackage.Doptim.Doptimize(self.arrayclass, nrestarts=2, optimfunc=optimfunc, verbose=1,
                                                              maxtime=18, selectpareto=False, nout=None, method=oalib.DOPTIM_UPDATE, niter=1000, nabort=0, dverbose=0)
        self.assertEqual(len(scores), n)
        self.assertEqual(len(dds), n)
        self.assertEqual(len(sols), n)
        scores, dds, sols, n = oapackage.Doptim.Doptimize(self.arrayclass, nrestarts=2, optimfunc=None, verbose=0,
                                                          maxtime=6, selectpareto=False, nout=None, method=oalib.DOPTIM_UPDATE, niter=30, nabort=0, dverbose=0)

    @only_python3
    def test_Doptimize_nonzero_strength(self):
        arrayclass = oapackage.arraydata_t(2, 16, 2, 6)
        with self.assertWarns(UserWarning):
            scores, dds, sols, _ = oapackage.Doptim.Doptimize(
                arrayclass, nrestarts=1, verbose=0, maxtime=9, selectpareto=False, nout=None, niter=1000, dverbose=0)

    def test_Doptimize_selectDn(self):
        scores, dds, sols, _ = oapackage.Doptim.Doptimize(self.arrayclass, nrestarts=10, optimfunc=[
                                                          1, 0, 0], verbose=0, maxtime=9, selectpareto=False, nout=None, method=oalib.DOPTIM_UPDATE, niter=1000, nabort=0, dverbose=0)

        result = oapackage.Doptim.selectDn(scores, dds, sols, nout=1, sortfull=True)
        self.assertTrue(len(result[2]) == 1)

    @only_python3
    def test_optimDeffPython(self):
        al = oapackage.exampleArray(2, 0)
        with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            _, al = oapackage.Doptim.optimDeffPython(
                al, arrayclass=None, niter=1000, nabort=1500, verbose=1, alpha=[1, 0, 0], method=0)
            stdout = mock_stdout.getvalue()
            self.assertTrue(stdout.startswith('optimDeff: initial Deff 0.0000'))

        for method in [oapackage.oalib.DOPTIM_SWAP, oapackage.oalib.DOPTIM_FLIP, oapackage.oalib.DOPTIM_UPDATE]:
            r, al = oapackage.Doptim.optimDeffPython(
                al, arrayclass=None, niter=100, nabort=200, verbose=0, alpha=[1, 0, 0], method=method)

    @only_python3
    def test_generateDscatter(self):
        if self.guitest:
            fig = 100
            with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
                r = oapackage.Doptim.generateDscatter(self.dds, second_index=0, first_index=1, lbls=None, verbose=1,
                                                      ndata=3, nofig=True, fig=fig, scatterarea=80)
        else:
            pass

    @only_python3
    def test_generateDpage(self):
        outputdir = tempfile.mkdtemp()
        allarrays = [oapackage.exampleArray(2, 0), oapackage.exampleArray(2, 0)]
        dds = np.array([A.Defficiencies() for A in allarrays])
        arrayclass = oapackage.arraylink2arraydata(allarrays[0])

        with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            page = oapackage.Doptim.generateDpage(outputdir, arrayclass, dds, allarrays,
                                                  fig=None, optimfunc=[1, 0, 0], nofig=True)
        if self.guitest and 0:
            print('test_generateDpage: run gui test')
            page = oapackage.Doptim.generateDpage(outputdir, arrayclass, dds, allarrays,
                                                  fig=100, optimfunc=[1, 0, 0], nofig=True)

            self.assertIsInstance(page, oapackage.markup.page)
            try:
                import matplotlib
                matplotlib.pyplot.close(100)
            except ImportError:
                pass

    def test_filterPareto(self):
        dds = self.dds2
        scores = np.arange(dds.shape[0])
        sols = [None] * scores.size
        s, _, _ = oapackage.Doptim.filterPareto(scores, dds, sols, verbose=0)
        self.assertEqual(list(s), [2, 3])

    def test_calcScore(self):
        dds = np.random.rand(10, 3)
        scores = oapackage.Doptim.calcScore(dds, optimfunc=[1, 2, 3])
        self.assertEqual(scores.shape, (dds.shape[0], ))

    def test_array2Dtable(self):
        sols = [oapackage.exampleArray(9, 0)]
        _ = oapackage.Doptim.array2Dtable(sols, verbose=1, titlestr=None)


if __name__ == '__main__':
    """ Test code """
    unittest.main()
    # t=TestDoptimize()
