"""Tests basic functionality for the SequenceScoringTool class.
Assumes that pandas is installed."""
import os
import unittest
import numpy as np
import pandas as pd
from antpack import SequenceScoringTool, SingleChainAnnotator

class TestSequenceScoringTool(unittest.TestCase):


    def test_error_checking(self):
        """Check that invalid data passed to the sequence scoring
        tool raises expected exceptions."""
        score_tool = SequenceScoringTool()

        # Sequences containing unrecognized characters should raise
        # a value error for these two routines, and return np.nan
        # for the other two.
        with self.assertRaises(ValueError):
            score_tool.get_closest_clusters("WoW")
        with self.assertRaises(ValueError):
            score_tool.suggest_mutations("WoW")

        self.assertTrue(np.isnan(score_tool.score_seq("WoW")))
        self.assertTrue(np.isnan(score_tool.batch_score_seqs(["WoW"])[0]))

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

        proba_test = score_tool.models["human"]["H"].predict_proba(np.zeros((1,173), dtype=np.uint8))
        self.assertTrue(len(proba_test.shape) == 2)
        self.assertTrue(proba_test.shape[1] == 1)
        self.assertTrue(proba_test.shape[0] == score_tool.models["human"]["H"].n_components)


    def test_constants(self):
        """Check that constants are what we expect."""
        score_tool = SequenceScoringTool()
        unadj_score_tool = SequenceScoringTool(False)

        self.assertTrue(unadj_score_tool.score_adjustments["H"] == 0)
        self.assertTrue(unadj_score_tool.score_adjustments["L"] == 0)

        self.assertTrue(score_tool.score_adjustments["H"] == -59.5617)
        self.assertTrue(score_tool.score_adjustments["L"] == -33.22445)


    def test_sequence_extraction(self):
        """Make sure the sequence prep feature is providing expected
        results."""
        score_tool = SequenceScoringTool()

        start_dir = os.path.abspath(os.path.dirname(__file__))

        raw_data = pd.read_csv(os.path.join(start_dir, "test_data", "imgt_comp_scoring.csv"))

        aligners = {"H":SingleChainAnnotator(["H"]),
            "K":SingleChainAnnotator(["K"]),
            "L":SingleChainAnnotator(["L"])  }

        for seq, chain_name in zip(raw_data["sequences"].tolist()[:10],
                raw_data["chain_types"].tolist()[:10]):
            p_extract, assigned_chain, _, _, _, _ = score_tool._prep_sequence(seq)
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
        unadj_score_tool = SequenceScoringTool(False)

        start_dir = os.path.abspath(os.path.dirname(__file__))

        raw_data = pd.read_csv(os.path.join(start_dir, "test_data", "imgt_comp_scoring.csv"))

        heavy_chains = raw_data[raw_data["chain_types"]=="H"]
        light_chains = raw_data[raw_data["chain_types"].isin(["K", "L"])]

        heavy_arr = [score_tool._prep_sequence(s)[2] for s in heavy_chains["sequences"].tolist()]
        light_arr = [score_tool._prep_sequence(s)[2] for s in light_chains["sequences"].tolist()]

        heavy_arr = np.vstack(heavy_arr)
        light_arr = np.vstack(light_arr)

        # Generate scores using a variety of different procedures. We will compare
        # these to make sure everything makes sense.

        raw_heavy_score = score_tool.models["human"]["H"].score(heavy_arr, n_threads = 2)
        raw_light_score = score_tool.models["human"]["L"].score(light_arr, n_threads = 2)

        unadj_heavy_score = unadj_score_tool.batch_score_seqs(heavy_chains["sequences"].tolist())
        unadj_light_score = unadj_score_tool.batch_score_seqs(light_chains["sequences"].tolist())
        self.assertTrue(np.allclose(unadj_heavy_score, raw_heavy_score))
        self.assertTrue(np.allclose(unadj_light_score, raw_light_score))

        adj_heavy_score = score_tool.batch_score_seqs(heavy_chains["sequences"].tolist())
        adj_light_score = score_tool.batch_score_seqs(light_chains["sequences"].tolist())
        self.assertTrue(np.allclose(adj_heavy_score, raw_heavy_score -
            score_tool.score_adjustments["H"]))
        self.assertTrue(np.allclose(adj_light_score, raw_light_score -
            score_tool.score_adjustments["L"]))

        self.assertTrue(np.allclose(adj_heavy_score, heavy_chains["batched_scores"].values))
        self.assertTrue(np.allclose(adj_light_score, light_chains["batched_scores"].values))

        gapped_heavy_score = score_tool.models["human"]["H"].gapped_score(heavy_arr, n_threads = 2)
        gapped_light_score = score_tool.models["human"]["L"].gapped_score(light_arr, n_threads = 2)
        self.assertTrue(np.allclose(heavy_chains["gapped_scores"].values, gapped_heavy_score -
            score_tool.score_adjustments["H"]))
        self.assertTrue(np.allclose(light_chains["gapped_scores"].values, gapped_light_score -
            score_tool.score_adjustments["L"]))

        self.assertTrue(np.allclose(raw_data["term_del_scores"].tolist(),
            [score_tool.score_seq(s, mask_term_dels=True) for s in raw_data["sequences"].tolist() ]))


    def test_masked_scoring(self):
        """Check that scoring using a mask is consistent with
        expectations. For simplicity, we just ask whether
        adding gaps and using the terminal mask scoring gives the
        same result as masking the terminal (it should, unless the
        deletion results in different numbering, which is possible --
        for this reason we manually call the categorical mix functions
        rather than using the score_tool wrapper."""
        test_seq = "EVHLQQSGAQELMKPGASVKISCKASGYTFITYWIEWVKQRPGHGLEWIGDILPGSGSTNYNENFKGKATFTADSSSNTAYMQLSSLTSEDSAVYYCARSGYYGNSGFAYWGQGTLVTVSA"
        test_annotator = SingleChainAnnotator(compress_init_gaps=False)
        score_tool = SequenceScoringTool(adjusted_scores = False,
                offer_classifier_option=False)
        score_tool.aligner = test_annotator
        numbering = test_annotator.analyze_seq(test_seq)[0]
        
        _, _, basic_arr, _, _, _ = score_tool._prep_sequence(test_seq)
        position_dict = score_tool.position_dict['H']

        mask1_pos = {'109', '110', '111', '112', '113', '114'}
        mask2_pos = {'2', '3', '4', '5', '6'}

        gapped_arr1, gapped_arr2 = basic_arr.copy(), basic_arr.copy()
        for mpos in list(mask1_pos):
            gapped_arr1[0, position_dict[mpos]] = 20
        for mpos in list(mask2_pos):
            gapped_arr2[0, position_dict[mpos]] = 20

        g1_score = float(score_tool.models["human"]['H'].gapped_score(
                    gapped_arr1, n_threads = 1)[0])
        g2_score = float(score_tool.models["human"]['H'].gapped_score(
                    gapped_arr2, n_threads = 1)[0])

        mask1 = [t not in mask1_pos for t in numbering]
        mask2 = [t not in mask2_pos for t in numbering]

        m1_score = score_tool.score_masked_sequences([test_seq],
                ["H"], [numbering], [mask1])
        m2_score = score_tool.score_masked_sequences([test_seq],
                ["H"], [numbering], [mask2])
        self.assertTrue(np.allclose(g1_score, m1_score))
        self.assertTrue(np.allclose(g2_score, m2_score))


if __name__ == "__main__":
    unittest.main()
