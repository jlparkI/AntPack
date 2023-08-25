"""Tests basic functionality for the antibody_scoring module and
toolkit."""
import sys
import unittest
import numpy as np
from scipy.linalg import hadamard

class TestFastHadamardTransform(unittest.TestCase):


    def test_3d_array_transform(self):
        self.assertTrue(outcome_f)

        dim = (250, 1, 4096)
        outcome_f, outcome_d = run_fht_test(dim)
        self.assertTrue(outcome_d)
        self.assertTrue(outcome_f)




if __name__ == "__main__":
    unittest.main()
