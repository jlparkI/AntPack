"""Tests scoring outputs of the c extension by comparing them
to a slow but easy to cross-check Python routine."""
import os
import unittest
import numpy as np
import pandas as pd
from antpack.utilities.special_functions import logsumexp
from antpack import SequenceScoringTool
from antpack_cpp_ext import mask_terminal_deletions


class TestCExtScoring(unittest.TestCase):


    def test_cext_scoring(self):
        """Check the workhorse function that provides scores
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

        raw_heavy_score = heavy_model.score(heavy_arr, n_threads = 2)
        raw_light_score = light_model.score(light_arr, n_threads = 2)

        comp_heavy_score = np.zeros((heavy_model.n_components, raw_heavy_score.shape[0]))
        comp_light_score = np.zeros((light_model.n_components, raw_light_score.shape[0]))

        for i in range(heavy_arr.shape[0]):
            for k in range(heavy_model.log_mu_mix.shape[0]):
                for j in range(heavy_model.log_mu_mix.shape[1]):
                    comp_heavy_score[k,i] += heavy_model.log_mu_mix[k,j,heavy_arr[i,j]]

        comp_heavy_score += heavy_model.log_mix_weights[:,None]
        comp_heavy_score = logsumexp(comp_heavy_score, axis=0)
        self.assertTrue(np.allclose(comp_heavy_score, raw_heavy_score))

        for i in range(light_arr.shape[0]):
            for k in range(light_model.log_mu_mix.shape[0]):
                for j in range(light_model.log_mu_mix.shape[1]):
                    comp_light_score[k,i] += light_model.log_mu_mix[k,j,light_arr[i,j]]

        comp_light_score += light_model.log_mix_weights[:,None]
        comp_light_score = logsumexp(comp_light_score, axis=0)
        self.assertTrue(np.allclose(comp_light_score, raw_light_score))


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
        mask_terminal_deletions(dummy_heavy_arr)
        mask_terminal_deletions(dummy_light_arr)

        for i in range(light_arr.shape[0]):
            for k in range(light_arr.shape[1]):
                if light_arr[i,k] != 20:
                    break
                light_arr[i,k] = 21
            for k in range(light_arr.shape[1]-1, 0, -1):
                if light_arr[i,k] != 20:
                    break
                light_arr[i,k] = 21

        for i in range(heavy_arr.shape[0]):
            for k in range(heavy_arr.shape[1]):
                if heavy_arr[i,k] != 20:
                    break
                heavy_arr[i,k] = 21
            for k in range(heavy_arr.shape[1]-1, 0, -1):
                if heavy_arr[i,k] != 20:
                    break
                heavy_arr[i,k] = 21

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

        comp_heavy_score = np.zeros((heavy_model.n_components, heavy_arr.shape[0]))
        comp_light_score = np.zeros((light_model.n_components, light_arr.shape[0]))

        for i in range(heavy_arr.shape[0]):
            for k in range(heavy_model.log_mu_mix.shape[0]):
                for j in range(heavy_model.log_mu_mix.shape[1]):
                    if not heavy_masks[j]:
                        continue
                    comp_heavy_score[k,i] += heavy_model.log_mu_mix[k,j,heavy_arr[i,j]]

        for i in range(light_arr.shape[0]):
            for k in range(light_model.log_mu_mix.shape[0]):
                for j in range(light_model.log_mu_mix.shape[1]):
                    if not light_masks[j]:
                        continue
                    comp_light_score[k,i] += light_model.log_mu_mix[k,j,light_arr[i,j]]

        heavy_masks = np.array(heavy_masks, dtype=bool)
        light_masks = np.array(light_masks, dtype=bool)

        raw_heavy_score = heavy_model.score(heavy_arr, heavy_masks, n_threads = 2)
        raw_light_score = light_model.score(light_arr, light_masks, n_threads = 2)

        comp_heavy_score += heavy_model.log_mix_weights[:,None]
        comp_heavy_score = logsumexp(comp_heavy_score, axis=0)
        self.assertTrue(np.allclose(comp_heavy_score, raw_heavy_score))

        comp_light_score += light_model.log_mix_weights[:,None]
        comp_light_score = logsumexp(comp_light_score, axis=0)
        self.assertTrue(np.allclose(comp_light_score, raw_light_score))



if __name__ == "__main__":
    unittest.main()
