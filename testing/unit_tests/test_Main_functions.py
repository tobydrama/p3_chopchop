from unittest import TestCase

from classes.ProgramMode import ProgramMode
from constants import CPF1_DEFAULT, TALEN_DEFAULT, CRISPR_DEFAULT, NICKASE_DEFAULT
from functions.main_functions import mode_select
from functions.set_default_modes import get_allowed_five_prime, get_mismatch_vectors, get_cpf1_mismatch_vectors


class Test(TestCase):
    def test_getAllowedFivePrime(self):
        self.assertEqual(("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG",
                          "TT"), get_allowed_five_prime("NN"))
        self.assertEqual({'AA': True, 'CA': True, 'GA': True, 'TA': True}, get_allowed_five_prime("NA"))
        self.assertEqual({'AA': True, 'AC': True, 'AG': True, 'AT': True}, get_allowed_five_prime("AN"))
        self.assertEqual({'AA': True}, get_allowed_five_prime("AA"))

    def test_getMismatchVectors(self):
        expected = ([True, True, True, True, True, True, True, False, False],
                    [True, True, True, True, True, True, False, False, False])
        self.assertEqual(expected, get_mismatch_vectors("NGG", 9, False))

        expected = ([True, True, True, True, True, True, True, True, True, True], [False])
        self.assertEqual(expected, get_mismatch_vectors("N", 1, True))
        self.assertEqual(([True], [False]), get_mismatch_vectors("N", 1, False))
        self.assertEqual(([False], [False]), get_mismatch_vectors("A", 1, False))
        self.assertEqual(([False], [False]), get_mismatch_vectors("A", 1, None))

    def test_getCpf1MismatchVectors(self):
        self.assertEqual(([False], [False]), get_cpf1_mismatch_vectors("A", 1))
        expected = ([True, False, False, True, True, True, True, True, True, True, True, True, True, True, True, True,
                     True, True, True, True, True, True, True, True, True, True, True],
                    [False, False, False, True, True, True, True, True, True, True, True, True, True, True, True, True,
                     True, True, True, True, True, True, True, True, True, True, True])
        self.assertEqual(expected, get_cpf1_mismatch_vectors("NGG", 27))
        self.assertEqual(([True, False, False], [False, False, False]), get_cpf1_mismatch_vectors("NGG", 1))

    def test_mode_select(self):
        self.assertEqual("?", mode_select("?", 1, None))
        self.assertNotEqual("?", mode_select(False, 1, None))

        self.assertEqual(CRISPR_DEFAULT["GUIDE_SIZE"], mode_select(None, "GUIDE_SIZE", ProgramMode.CRISPR))
        self.assertEqual(CRISPR_DEFAULT["PAM"], mode_select(None, "PAM", ProgramMode.CRISPR))
        self.assertEqual(CRISPR_DEFAULT["MAX_OFFTARGETS"], mode_select(None, "MAX_OFFTARGETS", ProgramMode.CRISPR))
        self.assertEqual(CRISPR_DEFAULT["MAX_MISMATCHES"], mode_select(None, "MAX_MISMATCHES", ProgramMode.CRISPR))
        self.assertEqual(CRISPR_DEFAULT["SCORE_GC"], mode_select(None, "SCORE_GC", ProgramMode.CRISPR))
        self.assertEqual(CRISPR_DEFAULT["SCORE_FOLDING"], mode_select(None, "SCORE_FOLDING", ProgramMode.CRISPR))

        self.assertEqual(CPF1_DEFAULT["GUIDE_SIZE"], mode_select(None, "GUIDE_SIZE", ProgramMode.CPF1))
        self.assertEqual(CPF1_DEFAULT["PAM"], mode_select(None, "PAM", ProgramMode.CPF1))
        self.assertEqual(CPF1_DEFAULT["MAX_OFFTARGETS"], mode_select(None, "MAX_OFFTARGETS", ProgramMode.CPF1))
        self.assertEqual(CPF1_DEFAULT["MAX_MISMATCHES"], mode_select(None, "MAX_MISMATCHES", ProgramMode.CPF1))
        self.assertEqual(CPF1_DEFAULT["SCORE_GC"], mode_select(None, "SCORE_GC", ProgramMode.CPF1))
        self.assertEqual(CPF1_DEFAULT["SCORE_FOLDING"], mode_select(None, "SCORE_FOLDING", ProgramMode.CPF1))

        self.assertEqual(TALEN_DEFAULT["GUIDE_SIZE"], mode_select(None, "GUIDE_SIZE", ProgramMode.TALENS))
        self.assertEqual(TALEN_DEFAULT["PAM"], mode_select(None, "PAM", ProgramMode.TALENS))
        self.assertEqual(TALEN_DEFAULT["MAX_OFFTARGETS"], mode_select(None, "MAX_OFFTARGETS", ProgramMode.TALENS))
        self.assertEqual(TALEN_DEFAULT["MAX_MISMATCHES"], mode_select(None, "MAX_MISMATCHES", ProgramMode.TALENS))
        self.assertEqual(TALEN_DEFAULT["SCORE_GC"], mode_select(None, "SCORE_GC", ProgramMode.TALENS))
        self.assertEqual(TALEN_DEFAULT["SCORE_FOLDING"], mode_select(None, "SCORE_FOLDING", ProgramMode.TALENS))

        self.assertEqual(NICKASE_DEFAULT["GUIDE_SIZE"], mode_select(None, "GUIDE_SIZE", ProgramMode.NICKASE))
        self.assertEqual(NICKASE_DEFAULT["PAM"], mode_select(None, "PAM", ProgramMode.NICKASE))
        self.assertEqual(NICKASE_DEFAULT["MAX_OFFTARGETS"], mode_select(None, "MAX_OFFTARGETS", ProgramMode.NICKASE))
        self.assertEqual(NICKASE_DEFAULT["MAX_MISMATCHES"], mode_select(None, "MAX_MISMATCHES", ProgramMode.NICKASE))
        self.assertEqual(NICKASE_DEFAULT["SCORE_GC"], mode_select(None, "SCORE_GC", ProgramMode.NICKASE))
        self.assertEqual(NICKASE_DEFAULT["SCORE_FOLDING"], mode_select(None, "SCORE_FOLDING", ProgramMode.NICKASE))
