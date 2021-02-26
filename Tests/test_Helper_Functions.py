from unittest import TestCase
from functions.Helper_Functions import gccontent, comaprePAM, permPAM


class Test(TestCase):
    def test_gccontent(self):
        self.assertEqual(0.8260869565217391, gccontent("CCTCCTGCTGCGCGCGCGCTCCC"))
        self.assertEqual(1, gccontent("CCG"))
        self.assertEqual(0.5, gccontent("gCGtiChFYg"))

    def test_comaprePAM(self):
        self.assertTrue(comaprePAM("N", "G"))
        self.assertTrue(comaprePAM("W", "T"))
        self.assertTrue(comaprePAM("W", "A"))
        self.assertTrue(comaprePAM("S", "C"))
        self.assertTrue(comaprePAM("S", "G"))
        self.assertTrue(comaprePAM("M", "A"))
        self.assertTrue(comaprePAM("M", "C"))
        self.assertTrue(comaprePAM("K", "G"))
        self.assertTrue(comaprePAM("K", "T"))
        self.assertTrue(comaprePAM("R", "A"))
        self.assertTrue(comaprePAM("R", "G"))
        self.assertTrue(comaprePAM("Y", "C"))
        self.assertTrue(comaprePAM("Y", "T"))
        self.assertTrue(comaprePAM("B", "H"))
        self.assertTrue(comaprePAM("D", "B"))
        self.assertTrue(comaprePAM("H", "S"))
        self.assertTrue(comaprePAM("V", "J"))
        self.assertFalse(comaprePAM("D", "C"))
        self.assertFalse(comaprePAM("H", "G"))
        self.assertFalse(comaprePAM("O", "S"))

    def test_permPAM(self):
        self.assertEqual(["AA", "AC", "AT", "AG"], permPAM("AN"))
        self.assertEqual(["AC", "TC", "AG", "TG"], permPAM("WS"))
        self.assertNotEqual(["T"], permPAM("A"))










