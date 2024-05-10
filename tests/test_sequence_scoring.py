"""Tests basic functionality for the SequenceScoringTool class.
Assumes that pandas is installed."""
import os
import unittest
import numpy as np
import pandas as pd
from antpack import SequenceScoringTool, SingleChainAnnotator, HumanizationTool

class TestSequenceScoringTool(unittest.TestCase):


    def test_error_checking(self):
        """Check that invalid data passed to the sequence scoring
        tool raises expected exceptions."""
        score_tool = SequenceScoringTool()
        humanization_tool = HumanizationTool()

        # Sequences containing unrecognized characters should raise
        # a value error for these two routines, and return np.nan
        # for the other two.
        with self.assertRaises(ValueError):
            humanization_tool.suggest_mutations("WoW")

        self.assertTrue(score_tool.get_diagnostic_info("WoW")[2]=="unknown")
        self.assertTrue(np.isnan(score_tool.score_seqs(["WoW"])[0]))

        with self.assertRaises(ValueError):
            score_tool.models["human"]["H"].score(np.zeros((1,173)))
        with self.assertRaises(ValueError):
            score_tool.models["human"]["H"].score(np.zeros((173), dtype=np.uint8))
        with self.assertRaises(ValueError):
            score_tool.models["human"]["H"].score(np.zeros((1,141), dtype=np.uint8))

        with self.assertRaises(ValueError):
            score_tool.models["human"]["L"].score(np.zeros((1,141)))
        with self.assertRaises(ValueError):
            score_tool.models["human"]["L"].score(np.zeros((141), dtype=np.uint8))
        with self.assertRaises(ValueError):
            score_tool.models["human"]["L"].score(np.zeros((1,173), dtype=np.uint8))

        proba_test = score_tool.models["human"]["H"].predict(np.zeros((1,173), dtype=np.uint8),
                return_raw_probs = True)
        self.assertTrue(len(proba_test.shape) == 2)
        self.assertTrue(proba_test.shape[1] == 1)
        self.assertTrue(proba_test.shape[0] == score_tool.models["human"]["H"].n_components)


    def test_sequence_extraction(self):
        """Make sure the sequence prep feature is providing expected
        results."""
        score_tool = SequenceScoringTool()

        start_dir = os.path.abspath(os.path.dirname(__file__))

        raw_data = pd.read_csv(os.path.join(start_dir, "test_data",
            "imgt_comp_scoring.csv.gz"))

        aligners = {"H":SingleChainAnnotator(["H"]),
            "K":SingleChainAnnotator(["K"]),
            "L":SingleChainAnnotator(["L"])  }

        for seq, chain_name in zip(raw_data["sequences"].tolist()[:10],
                raw_data["chain_types"].tolist()[:10]):
            p_extract, assigned_chain, _, _, _, _ = score_tool._full_prep_sequence(seq)
            if chain_name == "K":
                self.assertTrue(assigned_chain == "L")
            else:
                self.assertTrue(assigned_chain == chain_name)

            numbering = aligners[assigned_chain].analyze_seq(seq)[0]
            position_dict = score_tool.position_dict[assigned_chain]
            seq_extract = ["-" for i in range(len(position_dict))]

            for position, aa in zip(numbering, seq):
                if position == "-" or position not in position_dict:
                    continue
                seq_extract[position_dict[position]] = aa

            self.assertTrue(seq_extract == p_extract)




    def test_scoring_consistency(self):
        """Test that scores assigned by different procedures give what
        we expect."""
        score_tool = SequenceScoringTool(offer_classifier_option=True)
        adj_score_tool = SequenceScoringTool(offer_classifier_option=False,
                normalization="training_set_adjust")
        adj_score_tool.aligner = SingleChainAnnotator(compress_init_gaps=False)

        start_dir = os.path.abspath(os.path.dirname(__file__))

        raw_data = pd.read_csv(os.path.join(start_dir, "test_data",
            "imgt_comp_scoring.csv.gz"))

        heavy_chains = raw_data[raw_data["chain_types"]=="H"]
        light_chains = raw_data[raw_data["chain_types"].isin(["K", "L"])]

        heavy_arr = [score_tool._simple_prep_sequence(s)[1] for s
                in heavy_chains["sequences"].tolist()]
        light_arr = [score_tool._simple_prep_sequence(s)[1] for s
                in light_chains["sequences"].tolist()]

        heavy_arr = np.vstack(heavy_arr)
        light_arr = np.vstack(light_arr)

        # Generate scores using a variety of different procedures. We will compare
        # these to make sure everything makes sense.

        raw_heavy_score = score_tool.models["human"]["H"].score(heavy_arr, n_threads = 2)
        raw_light_score = score_tool.models["human"]["L"].score(light_arr, n_threads = 2)

        unadj_heavy_score = score_tool.score_seqs(heavy_chains["sequences"].tolist())
        unadj_light_score = score_tool.score_seqs(light_chains["sequences"].tolist())
        self.assertTrue(np.allclose(unadj_heavy_score, raw_heavy_score))
        self.assertTrue(np.allclose(unadj_light_score, raw_light_score))

        og_term_del_scores = raw_data["term_del_scores"].values
        term_del_scores = adj_score_tool.score_seqs(raw_data["sequences"].tolist(),
                mask_terminal_dels=True)
        self.assertTrue(np.allclose(term_del_scores, og_term_del_scores))



if __name__ == "__main__":
    unittest.main()
