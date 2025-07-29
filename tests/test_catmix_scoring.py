"""Tests scoring outputs of the categorical mixture model by comparing them
to a slow but easy to cross-check Python routine."""
import os
import unittest
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from antpack import SequenceScoringTool
from antpack.antpack_cpp_ext import py_logsumexp_axis0


class TestCatmixScoring(unittest.TestCase):


    def test_catmix_scoring(self):
        """Check the workhorse functions that provides scores
        and cluster assignments (predict and score)
        for consistency with a simple Python version."""
        start_dir = os.path.abspath(os.path.dirname(__file__))
        score_tool = SequenceScoringTool()
        raw_data = pd.read_csv(os.path.join(start_dir, "test_data",
            "imgt_comp_scoring.csv.gz"))

        heavy_chains = raw_data[raw_data["chain_types"]=="H"]
        light_chains = raw_data[raw_data["chain_types"].isin(["K", "L"])]

        heavy_arr = [score_tool._simple_prep_sequence(s)[1] for s
                in heavy_chains["sequences"].tolist()[:10]]
        light_arr = [score_tool._simple_prep_sequence(s)[1] for s
                in light_chains["sequences"].tolist()[:10]]

        heavy_arr = np.vstack(heavy_arr)
        light_arr = np.vstack(light_arr)
        heavy_model = score_tool.models["human"]["H"]
        light_model = score_tool.models["human"]["L"]


        gt_heavy_score = np.zeros((heavy_model.get_specs()[0],
            heavy_arr.shape[0]))
        gt_light_score = np.zeros((light_model.get_specs()[0],
            light_arr.shape[0]))

        log_mu_mix, log_mix_weights = heavy_model.get_model_parameters()
        log_mu_mix = np.log(log_mu_mix)
        log_mix_weights = np.log(log_mix_weights)

        for i in range(heavy_arr.shape[0]):
            for k in range(log_mu_mix.shape[0]):
                for j in range(log_mu_mix.shape[1]):
                    gt_heavy_score[k,i] += log_mu_mix[k,j,heavy_arr[i,j]]

        # Check predict proba.
        gt_heavy_score += log_mix_weights[:,None]
        test_heavy_score = heavy_model.predict_proba(heavy_arr)
        self.assertTrue(np.allclose(gt_heavy_score, test_heavy_score))

        # Now check predict.
        test_heavy_preds = heavy_model.predict(heavy_arr)
        gt_heavy_preds = np.argmax(gt_heavy_score, axis=0)
        self.assertTrue(np.allclose(gt_heavy_preds, test_heavy_preds))

        # Now check score.
        gt_heavy_score = logsumexp(gt_heavy_score, axis=0)
        test_heavy_score = heavy_model.score(heavy_arr)
        self.assertTrue(np.allclose(gt_heavy_score, test_heavy_score))



        # Do it again doing light chains (different set of test data and
        # parameters).
        log_mu_mix, log_mix_weights = light_model.get_model_parameters()
        log_mu_mix = np.log(log_mu_mix)
        log_mix_weights = np.log(log_mix_weights)

        for i in range(light_arr.shape[0]):
            for k in range(log_mu_mix.shape[0]):
                for j in range(log_mu_mix.shape[1]):
                    gt_light_score[k,i] += log_mu_mix[k,j,light_arr[i,j]]

        # Check predict_proba.
        gt_light_score += log_mix_weights[:,None]
        test_light_score = light_model.predict_proba(light_arr)
        self.assertTrue(np.allclose(gt_light_score, test_light_score))

        # Now check predict.
        gt_light_preds = np.argmax(gt_light_score, axis=0)
        test_light_preds = light_model.predict(light_arr)
        self.assertTrue(np.allclose(gt_light_preds, test_light_preds))

        # Now check score.
        gt_light_score = logsumexp(gt_light_score, axis=0)
        test_light_score = light_model.score(light_arr)
        self.assertTrue(np.allclose(gt_light_score, test_light_score))




    def test_terminal_deletions(self):
        """Tests the categorical mixture masked scoring function by
        comparing it with a simple python implementation."""
        start_dir = os.path.abspath(os.path.dirname(__file__))
        score_tool = SequenceScoringTool()
        raw_data = pd.read_csv(os.path.join(start_dir, "test_data",
            "imgt_comp_scoring.csv.gz"))

        heavy_chains = raw_data[raw_data["chain_types"]=="H"].sequences.tolist()
        light_chains = raw_data[raw_data["chain_types"].isin(["K", "L"])].sequences.tolist()

        heavy_chains[0] = heavy_chains[0][10:]
        heavy_chains[1] = heavy_chains[1][:-10]
        light_chains[0] = light_chains[0][10:]
        light_chains[1] = light_chains[1][:-25]

        heavy_arr = [score_tool._simple_prep_sequence(s)[1] for s
                in heavy_chains[:10]]
        light_arr = [score_tool._simple_prep_sequence(s)[1] for s
                in light_chains[:10]]
        heavy_arr = np.vstack(heavy_arr)
        light_arr = np.vstack(light_arr)

        dummy_heavy_arr = heavy_arr.copy()
        dummy_light_arr = light_arr.copy()
        score_tool.models["human"]["H"].mask_terminal_deletions(heavy_arr, dummy_heavy_arr)
        score_tool.models["human"]["H"].mask_terminal_deletions(light_arr, dummy_light_arr)

        for i in range(light_arr.shape[0]):
            for k in range(light_arr.shape[1]):
                if light_arr[i,k] != 20:
                    break
                light_arr[i,k] = 255
            for k in range(light_arr.shape[1]-1, 0, -1):
                if light_arr[i,k] != 20:
                    break
                light_arr[i,k] = 255

        for i in range(heavy_arr.shape[0]):
            for k in range(heavy_arr.shape[1]):
                if heavy_arr[i,k] != 20:
                    break
                heavy_arr[i,k] = 255
            for k in range(heavy_arr.shape[1]-1, 0, -1):
                if heavy_arr[i,k] != 20:
                    break
                heavy_arr[i,k] = 255

        self.assertTrue(np.allclose(dummy_heavy_arr, heavy_arr))
        self.assertTrue(np.allclose(dummy_light_arr, light_arr))



    def test_mask_scoring(self):
        """Ensure that the terminal deletion masking function yields
        expected results."""
        start_dir = os.path.abspath(os.path.dirname(__file__))
        score_tool = SequenceScoringTool()
        raw_data = pd.read_csv(os.path.join(start_dir, "test_data",
            "imgt_comp_scoring.csv.gz"))

        heavy_chains = raw_data[raw_data["chain_types"]=="H"]
        light_chains = raw_data[raw_data["chain_types"].isin(["K", "L"])]

        heavy_arr = [score_tool._simple_prep_sequence(s)[1] for s
                in heavy_chains["sequences"].tolist()[:10]]
        light_arr = [score_tool._simple_prep_sequence(s)[1] for s
                in light_chains["sequences"].tolist()[:10]]
        heavy_arr = np.vstack(heavy_arr)
        light_arr = np.vstack(light_arr)

        heavy_masks, light_masks = [], []

        heavy_masks = [(i < 10 or i > 20) for i in range(heavy_arr.shape[1])]
        light_masks = [(i < 10 or i > 20) for i in range(light_arr.shape[1])]

        heavy_model = score_tool.models["human"]["H"]
        light_model = score_tool.models["human"]["L"]

        gt_heavy_score = np.zeros((heavy_model.get_specs()[0], heavy_arr.shape[0]))
        gt_light_score = np.zeros((light_model.get_specs()[0], light_arr.shape[0]))

        log_mu_mix, log_mix_weights = heavy_model.get_model_parameters()
        log_mu_mix = np.log(log_mu_mix)
        log_mix_weights = np.log(log_mix_weights)

        for i in range(heavy_arr.shape[0]):
            for k in range(log_mu_mix.shape[0]):
                for j in range(log_mu_mix.shape[1]):
                    if not heavy_masks[j]:
                        continue
                    gt_heavy_score[k,i] += log_mu_mix[k,j,heavy_arr[i,j]]

        heavy_masks = np.array(heavy_masks, dtype=bool)
        test_heavy_score = heavy_model.score(heavy_arr, heavy_masks)

        gt_heavy_score += log_mix_weights[:,None]
        gt_heavy_score = logsumexp(gt_heavy_score, axis=0)
        self.assertTrue(np.allclose(gt_heavy_score, test_heavy_score))


        log_mu_mix, log_mix_weights = light_model.get_model_parameters()
        log_mu_mix = np.log(log_mu_mix)
        log_mix_weights = np.log(log_mix_weights)

        for i in range(light_arr.shape[0]):
            for k in range(log_mu_mix.shape[0]):
                for j in range(log_mu_mix.shape[1]):
                    if not light_masks[j]:
                        continue
                    gt_light_score[k,i] += log_mu_mix[k,j,light_arr[i,j]]

        light_masks = np.array(light_masks, dtype=bool)
        test_light_score = light_model.score(light_arr, light_masks)

        gt_light_score += log_mix_weights[:,None]
        gt_light_score = logsumexp(gt_light_score, axis=0)
        self.assertTrue(np.allclose(gt_light_score, test_light_score))



if __name__ == "__main__":
    unittest.main()
