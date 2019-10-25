""" Orthogonal Array package test functions
"""

import sys
import os
import logging
import unittest
import tempfile
import numpy as np

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


def is_sorted(l):
    return all(a <= b for a, b in zip(l, l[1:]))


def only_python3(function):
    if python3:
        def only_python3_function(*args, **kwargs):
            return function(*args, **kwargs)
    else:
        def only_python3_function(*args, **kwargs):
            return None
    return only_python3_function


import oapackage


class TestCompareMethods(unittest.TestCase):

    def test_compareLMC(self):
        al1 = oapackage.exampleArray(24, 0)
        al2 = oapackage.exampleArray(25, 0)
        al3 = oapackage.exampleArray(26, 0)

        self.assertEqual(oapackage.compareLMC(al1, al2), -1)
        self.assertEqual(oapackage.compareLMC(al2, al1), 1)
        self.assertEqual(oapackage.compareLMC(al1, al3), 1)

        al1 = oapackage.exampleArray(2, 0).selectFirstColumns(5)
        al2 = oapackage.exampleArray(2, 0).selectLastColumns(5)
        self.assertEqual(oapackage.compareLMC(al1, al2), -1)

    def test_compareJ54(self):
        al1 = oapackage.exampleArray(24, 0)
        al2 = oapackage.exampleArray(25, 0)
        al3 = oapackage.exampleArray(26, 0)

        self.assertEqual(oapackage.compareJ54(al1, al2), -1)
        self.assertEqual(oapackage.compareJ54(al2, al1), 1)

        self.assertEqual(oapackage.compareJ54(al1, al3), 1)
        self.assertEqual(oapackage.compareJ54(al2, al3), 1)


class TestMinimalFormCheck(unittest.TestCase):

    def test_LMCcheck(self):
        array = oapackage.exampleArray(1)

        lmc_type = oapackage.LMCcheck(array)
        self.assertEqual(lmc_type, oapackage.LMC_MORE)

        array2 = array.selectColumns([0, 1, 2, 4, 3])
        lmc_type = oapackage.LMCcheck(array2)
        self.assertEqual(lmc_type, oapackage.LMC_LESS)

        array2 = array.selectColumns([2, 1, 0, 3, 4])
        lmc_type = oapackage.LMCcheck(array2)
        self.assertEqual(lmc_type, oapackage.LMC_LESS)

        arrayclass = oapackage.arraylink2arraydata(array)
        reduction = oapackage.LMCreduction_t(arrayclass)
        oaextend = oapackage.OAextend(arrayclass)
        lmc_type = oapackage.LMCcheck(array2.array, arrayclass, oaextend, reduction)
        self.assertEqual(lmc_type, oapackage.LMC_LESS)

    def test_LMCcheckOriginal(self):
        array = oapackage.exampleArray(1)
        lmc_type = oapackage.LMCcheckOriginal(array)
        self.assertEqual(lmc_type, oapackage.LMC_MORE)

        array2 = array.selectColumns([0, 1, 2, 4, 3])
        lmc_type = oapackage.LMCcheck(array2)
        self.assertEqual(lmc_type, oapackage.LMC_LESS)

        array = oapackage.exampleArray(5, 0)
        with self.assertRaises(RuntimeError):
            lmc_type = oapackage.LMCcheckOriginal(array)


class TestReductions(unittest.TestCase):

    def test_LMC(self):
        al = oapackage.array_link(2, 2, 0)
        al[1, 1] = -1
        self.assertRaises(RuntimeError, oapackage.reduceLMCform, al)

        al = oapackage.exampleArray(8, 0)
        alr = oapackage.reduceLMCform(al)
        self.assertTrue(alr == al)

        array0 = oapackage.exampleArray(1, 0).selectFirstColumns(3)
        array = array0.randomperm()
        reduced_array = oapackage.reduceLMCform(array)
        self.assertEqual(array0, reduced_array)

    @only_python3
    def test_DOP(self):
        al = oapackage.exampleArray(1, 0)
        transformation = oapackage.reductionDOP(al)
        self.assertTrue(transformation.isIdentity())

        dof_values = oapackage.projectionDOFvalues(al)
        pvalues = [pvalue.raw_values() for pvalue in dof_values]
        self.assertTrue(is_sorted(pvalues))

        al = oapackage.exampleArray(7, 0)
        transformation = oapackage.reductionDOP(al)
        self.assertTrue(transformation.isIdentity())

        al = oapackage.exampleArray(8, 0)
        transformation = oapackage.reductionDOP(al)

        with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            transformation.show()
            stdout = mock_stdout.getvalue()
            self.assertTrue(stdout.startswith('array transformation: N 40'))
            self.assertIn('column permutation: {2,0,3,1,4,5,6}', stdout)
        alr = oapackage.reduceDOPform(al)
        self.assertTrue(transformation.apply(al) == alr)


class TestModelmatrix(unittest.TestCase):

    def test_modelmatrix_main_effects(self):
        # test model matrix main effects are helmert contrasts
        array = oapackage.array_link(np.array([[0, 1, 2, 3]]).T)

        M = oapackage.array2modelmatrix(array, 'm')
        hc = oapackage.oahelper.helmert_contrasts(4)
        np.testing.assert_array_almost_equal(hc, M[:, 1:])

    def test_modelmatrix(self):
        al = oapackage.exampleArray(1, 0)

        sizes = oapackage.array2modelmatrix_sizes(al)
        k = al.n_columns
        self.assertEqual(sizes, (1, 1 + k, int(1 + k + k * (k - 1) / 2), int(1 + k + k * (k + 1) / 2)))

        model_matrix = oapackage.array2modelmatrix(al, "main", verbose=0)
        np.testing.assert_array_equal(2 * np.array(al) - 1, model_matrix[:, 1:])

        conf_design = oapackage.exampleArray(41, 0)
        sizes = oapackage.array2modelmatrix_sizes(conf_design)
        k = conf_design.n_columns
        self.assertEqual(sizes, (1, 1 + k, int(1 + k + k * (k - 1) / 2), int(1 + k + k * (k + 1) / 2)))
        model_matrix = oapackage.array2modelmatrix(conf_design, "i", 0)
        self.assertTrue(model_matrix.shape[1] == sizes[2])
        model_matrix = oapackage.array2modelmatrix(conf_design, "q", 0)

        last_column = np.array(conf_design)[:, -1]

        np.testing.assert_array_equal(model_matrix[:, -1], last_column * last_column)

        self.assertTrue(model_matrix.shape[1] == sizes[3])

    def test_array2eigenModelMatrixMixed(self):
        array = oapackage.exampleArray(0, 0)
        r = oapackage.array2eigenModelMatrixMixed(array.selectFirstColumns(1), verbose=0)
        main_effects = r[0]
        np.testing.assert_array_equal(main_effects, np.array([[-1., -1., -1., -1., 1., 1., 1., 1.]]).T)

        array = oapackage.exampleArray(10, 0).selectFirstColumns(2)
        r = oapackage.array2eigenModelMatrixMixed(array, verbose=0)
        main_effects = r[0]
        interaction_model = r[1]
        np.testing.assert_array_almost_equal(main_effects, np.array([[-1.22474487, -0.70710678, -1.22474487, -0.70710678],
                                                                     [-1.22474487, -0.70710678, 1.22474487, -0.70710678],
                                                                     [-1.22474487, -0.70710678, 0., 1.41421356],
                                                                     [1.22474487, -0.70710678, 1.22474487, -0.70710678],
                                                                     [1.22474487, -0.70710678, 1.22474487, -0.70710678],
                                                                     [1.22474487, -0.70710678, 0., 1.41421356],
                                                                     [0., 1.41421356, -1.22474487, -0.70710678],
                                                                     [0., 1.41421356, -1.22474487, -0.70710678],
                                                                     [0., 1.41421356, 0., 1.41421356]]))
        np.testing.assert_array_almost_equal(interaction_model, np.array([[1.5, 0.8660254, 0.8660254, 0.5],
                                                                          [-1.5, 0.8660254, -0.8660254, 0.5],
                                                                          [-0., -1.73205081, -0., -1.],
                                                                          [1.5, -0.8660254, -0.8660254, 0.5],
                                                                          [1.5, -0.8660254, -0.8660254, 0.5],
                                                                          [0., 1.73205081, -0., -1.],
                                                                          [-0., -0., -1.73205081, -1.],
                                                                          [-0., -0., -1.73205081, -1.],
                                                                          [0., 0., 0., 2.]]))

    @only_python3
    def test_modelmatrix_verbosity(self):
        conf_design = oapackage.exampleArray(41, 0)
        with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            oapackage.array2modelmatrix(conf_design, "i", 1)
            stdout = mock_stdout.getvalue()
            self.assertIn('array2modelmatrix: type conference, model_type_idx 2', stdout)


class TestArrayLink(unittest.TestCase):

    @only_python3
    def test_selectFirstColumns(self):
        al = oapackage.exampleArray(41, 0)
        al = al.selectFirstColumns(3)
        self.assertEqual(al.n_columns, 3)

        al = oapackage.array_link()
        with self.assertRaises(RuntimeError):
            al = al.selectFirstColumns(1)

    def test_strength(self):
        example_strength_pairs = [(2, 2), (3, 3), (4, 2), (6, 2), (7, 2)]
        for idx, strength in example_strength_pairs:
            array = oapackage.exampleArray(idx, 0)
            self.assertEqual(array.strength(), strength)

    def test_basic_array_link_functionality(self):
        al2a = oapackage.array_link(2, 2, 0)
        al2b = oapackage.array_link(2, 2, 0)
        al2b.setconstant(1)

        al3 = oapackage.array_link(3, 2, 0)
        self.assertTrue(al2a != al2b)
        self.assertFalse(al2a == al2b)
        self.assertTrue(al2a.equalsize(al2b))
        self.assertFalse(al2a.equalsize(al3))
        self.assertTrue(al2a < al2b)
        self.assertFalse(al2b < al2a)

    def test_array_link_dimensions(self):
        for example_idx in [0, 2, 6, 10]:
            al = oapackage.exampleArray(example_idx, 0)
            self.assertEqual(al.size, al.n_rows * al.n_columns)
            self.assertEqual(al.shape, (al.n_rows, al.n_columns))
        with self.assertRaises(AttributeError):
            al.shape = 1

    def test_is_ortogonal_array(self):
        al = oapackage.exampleArray(1, 0)
        self.assertTrue(al.is_orthogonal_array())
        al[1, 2] = -1
        self.assertFalse(al.is_orthogonal_array())
        al = oapackage.exampleArray(41, 0)
        self.assertFalse(al.is_orthogonal_array())

    def test_array_class_functions(self):
        al = oapackage.exampleArray(1, 0)
        self.assertTrue(al.is2level())
        self.assertFalse(al.is_mixed_level())

        al = oapackage.exampleArray(11, 0)
        self.assertTrue(al.is2level())
        self.assertFalse(al.is_mixed_level())

        al = oapackage.array_link(2, 2, 0)
        al[1, 1] = -1
        self.assertFalse(al.is2level())
        self.assertFalse(al.is_mixed_level())


class TestArraydata_t(unittest.TestCase):

    def test_factor_levels(self):
        array = oapackage.exampleArray(5, 0)
        arrayclass = oapackage.arraylink2arraydata(array)
        factor_levels = arrayclass.factor_levels()
        self.assertEqual(factor_levels, (4, 3, 2, 2, 2))

    @only_python3
    def test_arraydata_t_oaindex(self):
        for ii in range(1, 4):
            arrayclass = oapackage.arraydata_t([2, 2, 2], 4 * ii, 2, 3)
            self.assertEqual(arrayclass.oaindex, ii)

        with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            arrayclass = oapackage.arraydata_t([4, 3, 3], 20, 2, 3)
            std_output = mock_stdout.getvalue()
            self.assertIn('arraydata_t: warning: no orthogonal arrays exist with the specified strength', std_output)
            self.assertEqual(arrayclass.oaindex, 0)

        with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            arrayclass = oapackage.arraydata_t([2, 3, 4], 20, 2, 3)
            std_output = mock_stdout.getvalue()
            self.assertIn('the factor levels of the structure are not sorted, this can lead to undefined behaviour', std_output)

        with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            arrayclass = oapackage.arraydata_t([6, 5], 10, 1, 2)
            self.assertEqual(arrayclass.oaindex, 0)
            std_output = mock_stdout.getvalue()


class TestJcharacteristics(unittest.TestCase):

    @only_python3
    def test_jstruct_conference(self):
        conf_design = oapackage.exampleArray(30, 0)
        js = oapackage.jstructconference_t(conf_design, 4)
        self.assertEqual(js.Jvalues(), (4, 0))

        with self.assertRaises(RuntimeError):
            js = oapackage.jstructconference_t(conf_design, 3)

        array = oapackage.array_link(10, 4, 0)
        with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            with self.assertRaises(Exception):
                js = oapackage.jstructconference_t(array, 4)
            std_output = mock_stdout.getvalue()
            self.assertIn('array is not conference design', std_output)

    def test_Jcharacteristics(self):
        al = oapackage.exampleArray(30, 0)
        jx = al.Jcharacteristics(4)
        self.assertEqual(jx, (0,))
        jx = al.Jcharacteristics(2)
        self.assertEqual(jx, (0, 0, 0, 0, 0, 0))

        al = oapackage.exampleArray(48, 0).selectFirstColumns(6)
        jx = al.Jcharacteristics(4)
        self.assertEqual(jx, (4, 4, 4, 0, 8, -8, 8, 4, -4, -4, 8, 4, -4, 8, -4))


class TestConferenceDesigns(unittest.TestCase):

    def test_conf2dsd(self):
        al = oapackage.exampleArray(42, 0)
        dsd = oapackage.conference2DSD(al)
        self.assertEqual(dsd.n_rows, 2 * al.n_rows + 1)

        conf = np.array(al)
        dsd = np.array(dsd)

        np.testing.assert_array_equal(conf, dsd[0:al.n_rows, :])
        np.testing.assert_array_equal(conf, -dsd[al.n_rows:2 * al.n_rows, :])
        np.testing.assert_array_equal(0 * conf[0], dsd[-1, :])

    def test_isConferenceFoldover(self):
        al = oapackage.exampleArray(37, 0)
        self.assertTrue(oapackage.isConferenceFoldover(al))
        al[2, 0] = -1
        self.assertFalse(oapackage.isConferenceFoldover(al))
        al = oapackage.exampleArray(7, 0)
        self.assertFalse(oapackage.isConferenceFoldover(al))
        al = oapackage.exampleArray(45, 0)
        self.assertFalse(oapackage.isConferenceFoldover(al))

    def test_double_conference_foldover_permutation(self):
        al = oapackage.exampleArray(37, 0)
        expected = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 11, 16, 19, 17, 18, 15, 14, 13, 12]
        permutation = oapackage.double_conference_foldover_permutation(al)
        self.assertEqual(list(permutation), expected)

        N = int(al.n_rows / 2)
        folded = np.array(oapackage.array_link(np.array(al)[permutation, :]))
        np.testing.assert_array_equal(folded[0:N, :], -folded[N:, :])

        expected[5], expected[3] = expected[3], expected[5]
        expected[-1], expected[4] = expected[4], expected[-1]
        expected[-3], expected[7] = expected[7], expected[-3]
        al = oapackage.array_link(np.array(al)[expected, :])
        permutation = oapackage.double_conference_foldover_permutation(al)
        folded = np.array(oapackage.array_link(np.array(al)[permutation, :]))
        N = int(al.n_rows / 2)
        np.testing.assert_array_equal(folded[0:N, :], -folded[N:, :])

        al = oapackage.exampleArray(45, 0)
        permutation = oapackage.double_conference_foldover_permutation(al)
        self.assertEqual(permutation[0], -1)

        al = oapackage.exampleArray(5, 0)
        with self.assertRaises(RuntimeError):
            permutation = oapackage.double_conference_foldover_permutation(al)

        al = oapackage.array_link(3, 4, 0)
        with self.assertRaises(RuntimeError):
            permutation = oapackage.double_conference_foldover_permutation(al)


class TestArrayFiles(unittest.TestCase):

    def test_write_latex_format(self):
        lst = [oapackage.exampleArray(2, 0)]
        filename = tempfile.mktemp(suffix='.tex')
        oapackage.writearrayfile(filename, oapackage.arraylist_t(lst), oapackage.ALATEX)
        with open(filename, 'rt') as fid:
            latex = fid.read()
        self.assertIsInstance(latex, str)


class TestCppLibrary(unittest.TestCase):

    def setUp(self):
        if getattr(self, 'assertRaisesRegex', None) is None:
            self.assertRaisesRegex = self.assertRaisesRegexp

    def test_root_form(self):
        array = oapackage.exampleArray(1, 0)
        self.assertTrue(oapackage.is_root_form(array, 2))
        self.assertFalse(oapackage.is_root_form(array, 5))

    @only_python3
    def test_runExtend_increasing_factor_levels(self):
        """ We test the usage of an increasing factor levels raises are warning, but nevertheless gives the correct results """
        N = 18
        k = 9
        t = 2
        l = [2, 3]
        rr = []
        with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            oapackage.oahelper.runExtend(N, k, t, l, verbose=1, nums=rr)
            self.assertIn(
                'warning: the factor levels of the structure are not sorted, this can lead to undefined behaviour', mock_stdout.getvalue())

        self.assertEqual(rr, [3, 15, 48, 19, 12, 3, 0])

    @only_python3
    def test_projectionDOFvalues(self):
        array = oapackage.exampleArray(5, 0)
        arrayclass = oapackage.arraylink2arraydata(array)
        dof_values = oapackage.projectionDOFvalues(array)
        sg = oapackage.symmetry_group(dof_values, False)
        with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            sg.show(1)
        self.assertIn('symmetry group: 5 elements, 5 subgroups: 1 1 1 1 1', mock_stdout.getvalue())

        for column in range(array.n_columns):
            dof_element = list(dof_values[column].raw_values())
            self.assertEqual(dof_element[0], -arrayclass.factor_levels()[column])

    @only_python3
    def test_exampleArray(self):
        with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            al = oapackage.exampleArray(5, 1)
            self.assertEqual(mock_stdout.getvalue(), 'exampleArray 5: array 0 in OA(24, 2, 4 3 2^a)\n')
            self.assertEqual(al.md5(), '3885c883d3bee0c7546511255bb5c3ae')

        with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            al = oapackage.exampleArray(-1, 1)
            lines = mock_stdout.getvalue().strip().split('\n')
            self.assertTrue(np.all([lines[ii].startswith('exampleArray %d:' % ii) for ii in range(len(lines))]))

        al = oapackage.exampleArray(51, 0)
        self.assertEqual(al.md5(), '662c1d9c51475b42539620385fa22338')

        with self.assertRaises(Exception):
            al = oapackage.exampleArray(1000, 0)

    def test_mvalue_t(self):
        input_vector = [1., 2., 2.]
        m = oapackage.mvalue_t_double(input_vector)
        self.assertEqual(m.size(), len(input_vector))
        self.assertEqual(list(m.values), input_vector)
        self.assertIsInstance(m.__array_interface__, dict)

        input_vector = [1, 2, -1]
        mvalue_long = oapackage.mvalue_t_long(input_vector)
        self.assertEqual(list(mvalue_long.values), input_vector)

    def test_splits(self):
        self.assertEqual(oapackage.splitTag([10, 12]), '10.12')
        self.assertEqual(oapackage.splitFile([]), '')
        self.assertEqual(oapackage.splitDir([1, 2]), 'sp0-split-1' + os.path.sep + 'sp1-split-2' + os.path.sep)

    @only_python3
    def test_array_transformation_t(self):
        at = oapackage.array_transformation_t()
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            at.show()
            std_output = mock_stdout.getvalue()
            self.assertEqual(std_output, 'array transformation: no class defined\n')

        al = oapackage.exampleArray(0, 0)
        arrayclass = oapackage.arraylink2arraydata(al)
        at = oapackage.array_transformation_t(arrayclass)
        at.setcolperm([1, 0])
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            _ = at.show()
            std_output = mock_stdout.getvalue()
            self.assertEqual(
                std_output, 'array transformation: N 8\ncolumn permutation: {1,0}\nlevel perms:\n{0,1}\n{0,1}\nrow permutation: {0,1,2,3,4,5,6,7}\n')
        self.assertEqual(at.colperm(), (1, 0))

    def test_arraylink2arraydata(self):
        al = oapackage.exampleArray(0, 0)
        adata = oapackage.arraylink2arraydata(al)
        self.assertEqual(str(adata), r'arrayclass: N 8, k 2, strength 2, s {2,2}, order 0')
        al = oapackage.array_link(4, 4, 0)
        al.setconstant(0)
        adata = oapackage.arraylink2arraydata(al)
        al.setconstant(-1)
        with self.assertRaises(RuntimeError):
            _ = oapackage.arraylink2arraydata(al)

        for ii in [0, 4, 8, 10]:
            al = oapackage.exampleArray(ii, 0)
            arrayclass = oapackage.arraylink2arraydata(al, strength=-1)
            self.assertEqual(arrayclass.strength, al.strength())

    def test_exception_handling(self):
        with self.assertRaises(RuntimeError):
            oapackage.mycheck_handler("file", "function", 10, 0, "hi")
        oapackage.mycheck_handler("file", "function", 10, 1, "hi")

        with self.assertRaises(RuntimeError):
            oapackage.throw_runtime_exception("should throw error")

        al = oapackage.oalib.exampleArray(18, 0)
        with self.assertRaisesRegex(RuntimeError, "array cannot have negative elements"):
            _ = oapackage.array2eigenModelMatrixMixed(al, 2)

    def test_array2eigenModelMatrixMixed(self):
        array = oapackage.exampleArray(5, 0)
        model = oapackage.array2eigenModelMatrixMixed(array)
        expected_main_effects = np.frombuffer(b'\xcc;\x7ff\x9e\xa0\xf6\xbf=,\x0cp\xbd \xea\xbf\x1c3\x90E\xa7y\xe2\xbf.!\t\x14\x8e\x98\xf3\xbf\xcc;\x7ff\x9e\xa0\xe6\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\xcc;\x7ff\x9e\xa0\xf6\xbf=,\x0cp\xbd \xea\xbf\x1c3\x90E\xa7y\xe2\xbf.!\t\x14\x8e\x98\xf3\xbf\xcc;\x7ff\x9e\xa0\xe6\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\xcc;\x7ff\x9e\xa0\xf6\xbf=,\x0cp\xbd \xea\xbf\x1c3\x90E\xa7y\xe2\xbf.!\t\x14\x8e\x98\xf3?\xcc;\x7ff\x9e\xa0\xe6\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\xcc;\x7ff\x9e\xa0\xf6\xbf=,\x0cp\xbd \xea\xbf\x1c3\x90E\xa7y\xe2\xbf.!\t\x14\x8e\x98\xf3?\xcc;\x7ff\x9e\xa0\xe6\xbf\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0?\xcc;\x7ff\x9e\xa0\xf6\xbf=,\x0cp\xbd \xea\xbf\x1c3\x90E\xa7y\xe2\xbf\x00\x00\x00\x00\x00\x00\x00\x00\xcc;\x7ff\x9e\xa0\xf6?\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0?\xcc;\x7ff\x9e\xa0\xf6\xbf=,\x0cp\xbd \xea\xbf\x1c3\x90E\xa7y\xe2\xbf\x00\x00\x00\x00\x00\x00\x00\x00\xcc;\x7ff\x9e\xa0\xf6?\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0?\xcc;\x7ff\x9e\xa0\xf6?=,\x0cp\xbd \xea\xbf\x1c3\x90E\xa7y\xe2\xbf.!\t\x14\x8e\x98\xf3\xbf\xcc;\x7ff\x9e\xa0\xe6\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0?\xcc;\x7ff\x9e\xa0\xf6?=,\x0cp\xbd \xea\xbf\x1c3\x90E\xa7y\xe2\xbf.!\t\x14\x8e\x98\xf3\xbf\xcc;\x7ff\x9e\xa0\xe6\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0?\xcc;\x7ff\x9e\xa0\xf6?=,\x0cp\xbd \xea\xbf\x1c3\x90E\xa7y\xe2\xbf.!\t\x14\x8e\x98\xf3?\xcc;\x7ff\x9e\xa0\xe6\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0?\xcc;\x7ff\x9e\xa0\xf6?=,\x0cp\xbd \xea\xbf\x1c3\x90E\xa7y\xe2\xbf.!\t\x14\x8e\x98\xf3?\xcc;\x7ff\x9e\xa0\xe6\xbf\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0\xbf\xcc;\x7ff\x9e\xa0\xf6?=,\x0cp\xbd \xea\xbf\x1c3\x90E\xa7y\xe2\xbf\x00\x00\x00\x00\x00\x00\x00\x00\xcc;\x7ff\x9e\xa0\xf6?\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\xcc;\x7ff\x9e\xa0\xf6?=,\x0cp\xbd \xea\xbf\x1c3\x90E\xa7y\xe2\xbf\x00\x00\x00\x00\x00\x00\x00\x00\xcc;\x7ff\x9e\xa0\xf6?\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\x00\x00=,\x0cp\xbd \xfa?\x1c3\x90E\xa7y\xe2\xbf.!\t\x14\x8e\x98\xf3\xbf\xcc;\x7ff\x9e\xa0\xe6\xbf\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\x00\x00=,\x0cp\xbd \xfa?\x1c3\x90E\xa7y\xe2\xbf.!\t\x14\x8e\x98\xf3\xbf\xcc;\x7ff\x9e\xa0\xe6\xbf\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\x00\x00=,\x0cp\xbd \xfa?\x1c3\x90E\xa7y\xe2\xbf.!\t\x14\x8e\x98\xf3?\xcc;\x7ff\x9e\xa0\xe6\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\x00\x00=,\x0cp\xbd \xfa?\x1c3\x90E\xa7y\xe2\xbf.!\t\x14\x8e\x98\xf3?\xcc;\x7ff\x9e\xa0\xe6\xbf\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\x00\x00=,\x0cp\xbd \xfa?\x1c3\x90E\xa7y\xe2\xbf\x00\x00\x00\x00\x00\x00\x00\x00\xcc;\x7ff\x9e\xa0\xf6?\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\x00\x00=,\x0cp\xbd \xfa?\x1c3\x90E\xa7y\xe2\xbf\x00\x00\x00\x00\x00\x00\x00\x00\xcc;\x7ff\x9e\xa0\xf6?\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xaaLX\xe8z\xb6\xfb?.!\t\x14\x8e\x98\xf3\xbf\xcc;\x7ff\x9e\xa0\xe6\xbf\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xaaLX\xe8z\xb6\xfb?.!\t\x14\x8e\x98\xf3\xbf\xcc;\x7ff\x9e\xa0\xe6\xbf\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xaaLX\xe8z\xb6\xfb?.!\t\x14\x8e\x98\xf3?\xcc;\x7ff\x9e\xa0\xe6\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xaaLX\xe8z\xb6\xfb?.!\t\x14\x8e\x98\xf3?\xcc;\x7ff\x9e\xa0\xe6\xbf\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xaaLX\xe8z\xb6\xfb?\x00\x00\x00\x00\x00\x00\x00\x00\xcc;\x7ff\x9e\xa0\xf6?\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xaaLX\xe8z\xb6\xfb?\x00\x00\x00\x00\x00\x00\x00\x00\xcc;\x7ff\x9e\xa0\xf6?\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\xf0?').reshape(24, 8)
        main_effects, interaction_effects = model[0], model[1]
        np.testing.assert_array_equal(main_effects, expected_main_effects)
        self.assertEqual(interaction_effects.shape, (24, 24))

    def test_mycheck_handler(self):
        oapackage.mycheck_handler('a', 'b', 1, 1, 'test mycheck_handler')
        with self.assertRaises(RuntimeError):
            oapackage.mycheck_handler('a', 'b', 1, 0, 'mycheck_handler raise')

    @only_python3
    def test_arrayrankInfo(self):
        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            rank = oapackage.arrayrankInfo(oapackage.exampleArray(21))
            self.assertEqual(rank, 4)
            self.assertIn('FullPivLU: rank 4', mock_stdout.getvalue())

    def test_rankStructure(self):
        array = oapackage.exampleArray(45, 0)
        array2 = oapackage.exampleArray(46, 0)
        rank_structure = oapackage.rankStructure(array)
        rank_structure.verbose = 2
        if python3:
            with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
                rank_structure.info()

        rank_structure.verbose = 0
        rank_direct = rank_structure.rankxf(array)
        rank = rank_structure.rankxf(array)
        rank2 = rank_structure.rankxf(array2)
        self.assertEqual(rank_direct, 20)
        self.assertEqual(rank, 20)
        self.assertEqual(rank2, 20)

    def distance_distribution(self):
        al = oapackage.array_link(2, 2, 0)
        distance_distrib = oapackage.distance_distribution(al)
        self.assertEqual(distance_distrib, (2.0, 0.0, 0.0))

        al[1, 0] = 1
        distance_distrib = oapackage.distance_distribution(al)
        self.assertEqual(distance_distrib, (1.0, 1.0, 0.0))

        al = oapackage.exampleArray(2, 0)

        distance_distrib = oapackage.distance_distribution(al)
        self.assertEqual(distance_distrib, (1.25, 0.75, 1.5, 6.5, 5.25, 0.75, 0.0))

    def test_Defficiencies(self):
        array = oapackage.exampleArray(0, 0)
        efficiencies = array.Defficiencies()
        self.assertEqual(efficiencies, (1.0, 1.0, 1.0))

        self.assertAlmostEqual(array.Eefficiency(), 1)
        self.assertAlmostEqual(array.VIFefficiency(), 1)

        array = oapackage.exampleArray(11, 0)
        efficiencies = array.Defficiencies()
        self.assertEqual(efficiencies, (0.8879176205539137, 0.8058954593581515, 0.9936928808878267))

        self.assertAlmostEqual(array.Eefficiency(), 0.3602369388406959)
        self.assertAlmostEqual(array.VIFefficiency(), 1.2648199507220503)

        for idx in [0, 4, 11, 16]:
            array = oapackage.exampleArray(idx, 0)
            efficiencies = array.Defficiencies()
            D = array.Defficiency()
            Ds = array.DsEfficiency()
            self.assertAlmostEqual(efficiencies[0], D)
            self.assertAlmostEqual(efficiencies[1], Ds)

    @only_python3
    def test_projection_efficiencies(self):
        al = oapackage.exampleArray(11, 0)
        with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            d = oapackage.projDeff(al, 3, 1)
            D = al.selectFirstColumns(3).Defficiency()
            std_output = mock_stdout.getvalue()
            self.assertIn('projDeff: k 8, kp 3: start with 56 combinations', std_output)
        self.assertEqual(D, d[0])
        np.testing.assert_almost_equal(np.mean(d), 0.99064112542249538329031111061340197921)

        pec_seq = oapackage.PECsequence(al)
        np.testing.assert_equal(pec_seq, (1.0,) * len(pec_seq))
        pic_seq = oapackage.PICsequence(al)
        np.testing.assert_almost_equal(pic_seq, (0.9985780064264659, 0.9965697009006985, 0.9906411254224957,
                                                 0.9797170906488152, 0.9635206782887167, 0.9421350381959234, 0.9162739059686846, 0.8879176205539139))

    def test_PICsequence_length(self):
        array = oapackage.exampleArray(8, 0)
        pic = oapackage.PICsequence(array)
        self.assertEqual(len(pic), array.n_columns)

    def test_arraylink(self):
        al = oapackage.exampleArray(0, 0)
        self.assertEqual(al.at(0, 1), 0)
        self.assertEqual(al.at(0), 0)
        with self.assertRaises(IndexError):
            al.at(-1, 0)
        with self.assertRaises(IndexError):
            al.at(-1)
        with self.assertRaises(IndexError):
            al.at(0, al.n_columns)

    def test_arraylink_slicing(self):
        numpy_array = np.arange(0, 6 * 10).reshape((6, 10))

        al = oapackage.makearraylink(numpy_array)
        self.assertTrue(al[0] == numpy_array.flatten()[0])
        self.assertTrue(al[0, 1] == numpy_array[0, 1])
        self.assertTrue(al[4, 2] == numpy_array[4, 2])
        np.testing.assert_equal(al[0:4, 1:5], np.array(al)[0:4, 1:5])
        np.testing.assert_equal(al[0:1, 0:10:2], np.array(al)[0:1, 0:10:2])
        np.testing.assert_equal(al[3, 3::], np.array(al)[3:4, 3::])
        np.testing.assert_equal(al[2, :8:2], np.array(al)[2:3, :8:2])
        np.testing.assert_equal(al[2:3, :8:2], np.array(al)[2:3, :8:2])

        with self.assertRaises(IndexError):
            _ = al[-1, 1]

    def test_conference_generation(self):

        al = oapackage.exampleArray(42, 0)
        lst = [al]
        conference_type = oapackage.conference_t(al.n_rows, al.n_rows, 0)

        extensions = oapackage.extend_conference_restricted(lst, conference_type, verbose=1)
        self.assertEqual(len(extensions), 10)

        self.assertEqual(extensions[0].md5(), 'f759e75d3ce6adda5489fed4c528a6fb')
        self.assertEqual(extensions[-1].md5(), 'ec13ed02c70dbc9975a99e70e23154b3')


if __name__ == '__main__':
    """ Test code """
    unittest.main()
