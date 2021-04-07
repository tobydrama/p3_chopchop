from unittest import TestCase
from functions.eval import gccontent, comparePAM, permPAM
from functions.set_default_modes import getAllowedFivePrime


class Test(TestCase):
    def test_gccontent(self):
        self.assertEqual(0.8260869565217391, gccontent("CCTCCTGCTGCGCGCGCGCTCCC"))
        self.assertEqual(1, gccontent("CCG"))
        self.assertEqual(0.5, gccontent("gCGtiChFYg"))

    def test_comaprePAM(self):
        self.assertTrue(comparePAM("N", "G"))
        self.assertTrue(comparePAM("W", "T"))
        self.assertTrue(comparePAM("W", "A"))
        self.assertTrue(comparePAM("S", "C"))
        self.assertTrue(comparePAM("S", "G"))
        self.assertTrue(comparePAM("M", "A"))
        self.assertTrue(comparePAM("M", "C"))
        self.assertTrue(comparePAM("K", "G"))
        self.assertTrue(comparePAM("K", "T"))
        self.assertTrue(comparePAM("R", "A"))
        self.assertTrue(comparePAM("R", "G"))
        self.assertTrue(comparePAM("Y", "C"))
        self.assertTrue(comparePAM("Y", "T"))
        self.assertTrue(comparePAM("B", "H"))
        self.assertTrue(comparePAM("D", "B"))
        self.assertTrue(comparePAM("H", "S"))
        self.assertTrue(comparePAM("V", "J"))
        self.assertFalse(comparePAM("D", "C"))
        self.assertFalse(comparePAM("H", "G"))
        self.assertFalse(comparePAM("O", "S"))

    def test_permPAM(self):
        self.assertEqual(["AA", "AC", "AT", "AG"], permPAM("AN"))
        self.assertEqual(["AC", "TC", "AG", "TG"], permPAM("WS"))
        self.assertNotEqual(["T"], permPAM("A"))

    def test_getAllowedFivePrime(self):
        self.assertEqual(("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG",
                          "TT"), getAllowedFivePrime("NN"))
