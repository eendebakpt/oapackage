import unittest

import oapackage


class TestGWLP(unittest.TestCase):

    def test_gwlp(self):
        array = oapackage.exampleArray(1, 0)
        array = array.selectFirstColumns(3)
        gwlp = oapackage.GWLP(array)
        gwlp_mixed = oapackage.GWLPmixed(array, 0)

        self.assertEqual(gwlp,  (1.0, 0.0, 0.0, 0.25))
        self.assertEqual(gwlp, gwlp_mixed)

    def test_gwlp_mixed(self):
        array = oapackage.exampleArray(56, 0)
        gwlp = oapackage.GWLPmixed(array, 0)
        self.assertEqual(gwlp,  (1.0, 0.0, 0.0, 48.0, 185.0, 524.0, 1381.0, 2792.0, 4327.0, 5556.0,
                         5961.0, 5152.0, 3539.0, 1972.0, 911.0, 328.0, 76.0, 12.0, 3.0, 0.0))

    def test_gwlp_mixed2(self):
        al = oapackage.exampleArray(4, 0)
        self.assertEqual(al.GWLP(), (1.0, 0.0, 0.0, 3.5, 2.5, 0.5, 0.5, 0.0))

    def test_regression_gwlp_non_mixed(self):
        al = oapackage.exampleArray(57, 0)
        gwlp = al.GWLP()
        self.assertEqual(gwlp, (1.0, 0.0, 0.0, 13.0, 13.5, 9.0, 4.0))
        self.assertEqual(oapackage.GWLPmixed(al), (1.0, 0.0, 0.0, 13.0, 13.5, 9.0, 4.0))


if __name__ == '__main__':
    """ Test code """
    unittest.main()
