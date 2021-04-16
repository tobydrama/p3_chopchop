from unittest import TestCase
from functions.eval import gc_content, compare_PAM, perm_PAM
from functions.set_default_modes import get_allowed_five_prime


class Test(TestCase):
    def test_gccontent(self):
        self.assertEqual(0.8260869565217391, gc_content("CCTCCTGCTGCGCGCGCGCTCCC"))
        self.assertEqual(1, gc_content("CCG"))
        self.assertEqual(0.5, gc_content("gCGtiChFYg"))

    def test_comaprePAM(self):
        self.assertTrue(compare_PAM("N", "G"))
        self.assertTrue(compare_PAM("W", "T"))
        self.assertTrue(compare_PAM("W", "A"))
        self.assertTrue(compare_PAM("S", "C"))
        self.assertTrue(compare_PAM("S", "G"))
        self.assertTrue(compare_PAM("M", "A"))
        self.assertTrue(compare_PAM("M", "C"))
        self.assertTrue(compare_PAM("K", "G"))
        self.assertTrue(compare_PAM("K", "T"))
        self.assertTrue(compare_PAM("R", "A"))
        self.assertTrue(compare_PAM("R", "G"))
        self.assertTrue(compare_PAM("Y", "C"))
        self.assertTrue(compare_PAM("Y", "T"))
        self.assertTrue(compare_PAM("B", "H"))
        self.assertTrue(compare_PAM("D", "B"))
        self.assertTrue(compare_PAM("H", "S"))
        self.assertTrue(compare_PAM("V", "J"))
        self.assertFalse(compare_PAM("D", "C"))
        self.assertFalse(compare_PAM("H", "G"))
        self.assertFalse(compare_PAM("O", "S"))

    def test_permPAM(self):
        self.assertEqual(["AA", "AC", "AT", "AG"], perm_PAM("AN"))
        self.assertEqual(["AC", "TC", "AG", "TG"], perm_PAM("WS"))
        self.assertNotEqual(["T"], perm_PAM("A"))

    def test_getAllowedFivePrime(self):
        self.assertEqual(("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG",
                          "TT"), get_allowed_five_prime("NN"))
