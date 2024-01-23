"""Tests basic functionality for the SequenceScoringTool class.
Assumes that pandas is installed."""
import os
import unittest
import numpy as np
import pandas as pd
from antpack import SequenceScoringTool

class TestSequenceScoringTool(unittest.TestCase):


    def test_error_checking(self):
        """Check that invalid data passed to the sequence scoring
        tool raises expected exceptions."""
        score_tool = SequenceScoringTool()

        self.assertTrue(np.isnan(score_tool.score_seq("WoW")))
        self.assertTrue(np.isnan(score_tool.batch_score_seqs(["WoW"])[0]))
        
        with self.assertRaises(ValueError):
            score_tool.per_position_probs("WoW")
        with self.assertRaises(ValueError):
            score_tool.get_closest_clusters("WoW")
        with self.assertRaises(ValueError):
            score_tool.suggest_mutations("WoW")

        with self.assertRaises(ValueError):
            score_tool.models["human"]["H"].score(np.zeros((1,173)))
        with self.assertRaises(ValueError):
            score_tool.models["human"]["H"].score(np.zeros((173), dtype=np.uint8))
        with self.assertRaises(ValueError):
            score_tool.models["human"]["H"].score(np.zeros((1,141), dtype=np.uint8))

        with self.assertRaises(ValueError):
            score_tool.models["human"]["L"].score(np.zeros((1,173)))
        with self.assertRaises(ValueError):
            score_tool.models["human"]["L"].score(np.zeros((173), dtype=np.uint8))
        with self.assertRaises(ValueError):
            score_tool.models["human"]["L"].score(np.zeros((1,173), dtype=np.uint8))


    def basic_checks(self):
        """Check that constants are what we expect."""
        score_tool = SequenceScoringTool()
        unadj_score_tool = SequenceScoringTool(False)

        self.assertTrue(unadj_score_tool.score_adjustments["H"] == 0)
        self.assertTrue(unadj_score_tool.score_adjustments["L"] == 0)

        self.assertTrue(score_tool.score_adjustments["H"] == -59.5617)
        self.assertTrue(score_tool.score_adjustments["L"] == -33.22445)


    def check_sequence_extraction(self):
        """Make sure the sequence prep feature is providing expected
        results."""


    def test_scoring_consistency(self):
        """Test that scores assigned by different procedures give what
        we expect."""
        score_tool = SequenceScoringTool(offer_classifier_option=True)
        unadj_score_tool = SequenceScoringTool(False)

        start_dir = os.path.abspath(os.path.dirname(__file__))

        raw_data = pd.read_csv(os.path.join(start_dir, "test_data", "imgt_comp_scoring.csv"))
        sequences = raw_data["sequences"].tolist()

        # First, make sure scores assigned by score tool are equivalent to historical.
        self.assertTrue(np.allclose(raw_data["term_del_scores"].tolist(),
            [score_tool.score_seq(s, mask_term_dels=True) for s in sequences]))
        self.assertTrue(np.allclose(raw_data["batched_scores"].values,
            score_tool.batch_score_seqs(sequences) ))
        self.assertTrue(np.allclose(raw_data["gapped_scores"].values,
            score_tool.batch_score_seqs(sequences, mask_heavy_gaps = True,
                mask_light_gaps = True) ))



if __name__ == "__main__":
    unittest.main()
