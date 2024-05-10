"""Tests basic functionality for the LiabilitySearchTool class."""
import os
import unittest
from antpack import LiabilitySearchTool

class TestLiabilitySearchTool(unittest.TestCase):


    def test_liability_tool(self):
        """Check that sequences with liabilities are flagged as such
        while other sequences are passed."""
        search_tool = LiabilitySearchTool()

        dummy_sequence = "AATHSCSSC"
        result = search_tool.analyze_seq(dummy_sequence)
        self.assertTrue(result[0][1] == "Sequence contains invalid characters")

        generic_test_seq = "VHLQQSGAELMKPGASVKISCKASGYTFITYWIEWVKQRPGHGLEWIGDILPGSGSTNYNENFKGKATFTADSSSNTAYMQLSSLTSEDSAVYYCARSGYYGNSGFAYWGQGTLVTVSA"
        result = search_tool.analyze_seq(generic_test_seq)
        self.assertTrue(len(result)==3)
        self.assertTrue(result[0][1] == "Deamidation (elevated risk)")
        self.assertTrue(result[1][1] == "Methionine / Tryptophan oxidation moderate risk")
        self.assertTrue(result[2][1] == "Deamidation (low risk)")

        corrected_seq = list(generic_test_seq)
        corrected_seq[102:104] = "AL"
        self.assertTrue(len(search_tool.analyze_seq("".join(corrected_seq))) == 2)
        corrected_seq[31] = "A"
        self.assertTrue(len(search_tool.analyze_seq("".join(corrected_seq))) == 1)
        corrected_seq[56:58] = "AL"
        self.assertTrue(len(search_tool.analyze_seq("".join(corrected_seq))) == 0)

        liability_seq = corrected_seq.copy()
        liability_seq[3] = "C"
        result = search_tool.analyze_seq("".join(liability_seq))
        self.assertTrue(result[0][1] == "Unusual cysteine")

        liability_seq = corrected_seq.copy()
        liability_seq[-5:-2] = "NAS"
        result = search_tool.analyze_seq("".join(liability_seq))
        self.assertTrue(result[0][1] == "N-glycosylation")

        liability_seq = corrected_seq.copy()
        liability_seq[3:5] = "M"
        result = search_tool.analyze_seq("".join(liability_seq))
        self.assertTrue(len(result) == 0)


if __name__ == "__main__":
    unittest.main()
