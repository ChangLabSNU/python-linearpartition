from unittest import TestCase
from linearpartition import partition

class TestLinearPartition(TestCase):

    seq1 = 'UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA'
    probmtx1 = [
        (0, 71, 0.5397), (1, 70, 0.9983), (2, 69, 0.9994), (3, 68, 0.9998),
        (4, 67, 0.9982), (5, 66, 0.9993), (5, 67, 0.0002), (6, 65, 0.9061),
        (6, 66, 0.0003), (7, 64, 0.9293), (8, 63, 0.9976), (9, 62, 0.9999),
        (10, 61, 0.9999), (11, 60, 0.9977), (12, 59, 0.9955), (13, 58, 0.9979),
        (14, 57, 0.9983), (15, 56, 0.9899), (17, 54, 0.0001), (17, 55, 0.9982),
        (18, 54, 0.9995), (19, 53, 0.9998), (20, 52, 0.9993), (21, 50, 0.0011),
        (21, 51, 0.9947), (23, 49, 0.9764), (24, 48, 0.9997), (25, 47, 0.9989),
        (26, 46, 0.9986), (27, 45, 0.9975), (29, 43, 0.0271), (29, 44, 0.9506),
        (30, 41, 0.0214), (30, 43, 0.9369), (30, 44, 0.0004), (30, 46, 0.0000),
        (31, 38, 0.0007), (31, 40, 0.0248), (31, 42, 0.7542), (31, 45, 0.0001),
        (32, 41, 0.4639), (32, 44, 0.0001), (33, 41, 0.4888), (34, 38, 0.0054),
        (34, 39, 0.0122), (34, 40, 0.7858), (35, 39, 0.5950), (37, 44, 0.0001),
        (38, 43, 0.0001)]
    freeenergy1 = -36.0677
    mea_structure = '((((((((((((((((.(((((.(((((.(((.(((...))))))))))).)))))))))))))))))))))'

    def test_1(self):
        pred = partition(self.seq1)
        for expected, actual in zip(self.probmtx1, sorted(pred['bpp'].tolist())):
            self.assertEqual(expected[:2], actual[:2])
            self.assertAlmostEqual(expected[2], actual[2], places=4)
        self.assertAlmostEqual(pred['free_energy'], self.freeenergy1, places=4)
        self.assertEqual(pred['structure'], self.mea_structure)

if __name__ == '__main__':
    import unittest
    unittest.main()
