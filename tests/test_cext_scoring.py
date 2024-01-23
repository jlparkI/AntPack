"""Tests scoring outputs of the c extension by comparing them
to a slow but easy to cross-check Python routine."""
import os
import unittest
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from antpack import SequenceScoringTool, SingleChainAnnotator


class TestCExtScoring(unittest.TestCase):


    def test_cext_scoring(self):
        """Check the workhorse function that provides scores
        for consistency with a simple Python version."""
        start_dir = os.path.abspath(os.path.dirname(__file__))
        score_tool = SequenceScoringTool()
        raw_data = pd.read_csv(os.path.join(start_dir, "test_data", "imgt_comp_scoring.csv"))

        heavy_chains = raw_data[raw_data["chain_types"]=="H"]
        light_chains = raw_data[raw_data["chain_types"].isin(["K", "L"])]

        heavy_arr = [score_tool._prep_sequence(s)[2] for s in heavy_chains["sequences"].tolist()[:10]]
        light_arr = [score_tool._prep_sequence(s)[2] for s in light_chains["sequences"].tolist()[:10]]

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


if __name__ == "__main__":
    unittest.main()
