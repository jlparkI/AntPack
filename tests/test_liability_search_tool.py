"""Tests basic functionality for the LiabilitySearchTool class."""
import os
import unittest
from antpack import LiabilitySearchTool, SingleChainAnnotator

class TestLiabilitySearchTool(unittest.TestCase):


    def test_liability_tool(self):
        """Check that sequences with liabilities are flagged as such
        while other sequences are passed."""
        search_tool = LiabilitySearchTool()
        annotation_tool = SingleChainAnnotator()

        generic_test_seq = "VHLQQSGAELMKPGASVKISCKASGYTFITYWIEWVKQRPGHGLEWIGDILPGSGSTNYNENFKGKATFTADSSSNTAYMQLSSLTSEDSAVYYCARSGYYGNSGFAYWGQGTLVTVSA"
        alignment = annotation_tool.analyze_seq(generic_test_seq)
        result = search_tool.analyze_seq(generic_test_seq, alignment,
                "imgt", "imgt")
        self.assertTrue(len(result)==3)
        self.assertTrue(result[0][1] == "Methionine / Tryptophan oxidation moderate risk")
        self.assertTrue(result[1][1] == "Deamidation (low risk)")
        self.assertTrue(result[2][1] == "Deamidation (elevated risk)")

        corrected_seq = list(generic_test_seq)
        corrected_seq[102:104] = "AL"
        corrected_seq = "".join(corrected_seq)
        alignment = annotation_tool.analyze_seq(corrected_seq)
        self.assertTrue(len(search_tool.analyze_seq(corrected_seq,
            alignment, "imgt", "imgt")) == 2)

        corrected_seq = list(corrected_seq)
        corrected_seq[31] = "A"
        corrected_seq = "".join(corrected_seq)
        alignment = annotation_tool.analyze_seq(corrected_seq)
        self.assertTrue(len(search_tool.analyze_seq(corrected_seq,
            alignment, "imgt", "imgt")) == 1)

        corrected_seq = list(corrected_seq)
        corrected_seq[56:58] = "AL"
        corrected_seq = "".join(corrected_seq)
        alignment = annotation_tool.analyze_seq(corrected_seq)
        self.assertTrue(len(search_tool.analyze_seq(corrected_seq,
            alignment, "imgt", "imgt")) == 0)

        corrected_seq = list(corrected_seq)
        liability_seq = corrected_seq.copy()
        liability_seq[3] = "C"
        liability_seq = "".join(liability_seq)
        alignment = annotation_tool.analyze_seq(liability_seq)
        result = search_tool.analyze_seq(liability_seq,
                alignment, "imgt", "imgt")
        self.assertTrue(result[0][1] == "Unusual cysteine")

        liability_seq = corrected_seq.copy()
        liability_seq[-5:-2] = "NAS"
        liability_seq = "".join(liability_seq)
        alignment = annotation_tool.analyze_seq(liability_seq)
        result = search_tool.analyze_seq(liability_seq,
                alignment, "imgt", "imgt")
        self.assertTrue(result[0][1] == "N-glycosylation risk")

        liability_seq = corrected_seq.copy()
        liability_seq[3:5] = "M"
        liability_seq = "".join(liability_seq)
        alignment = annotation_tool.analyze_seq(liability_seq)
        result = search_tool.analyze_seq(liability_seq,
                alignment, "imgt", "imgt")
        self.assertTrue(len(result) == 0)


if __name__ == "__main__":
    unittest.main()
