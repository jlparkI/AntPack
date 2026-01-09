"""Tests scoring outputs of the categorical mixture model by comparing them
to a slow but easy to cross-check Python routine."""
import os
import pytest
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from antpack.scoring_tools.scoring_constants import allowed_imgt_pos as ahip
from antpack import (SequenceScoringTool, SingleChainAnnotator,
        SequenceTemplateAligner)
from ..conftest import (standard_aa_list, get_test_data_filepath)



def test_catmix_scoring(scoring_comp_sequences):
    """Check the workhorse functions that provides scores
    and cluster assignments (predict and score)
    for consistency with a simple Python version."""
    heavy_seqs, heavy_arr, light_seqs, light_arr = scoring_comp_sequences
    score_tool = SequenceScoringTool()
    heavy_model = score_tool.models["human"]["H"]
    light_model = score_tool.models["human"]["L"]

    gt_heavy_score = np.zeros((heavy_model.get_model_specs()[0],
        heavy_arr.shape[0]))
    gt_light_score = np.zeros((light_model.get_model_specs()[0],
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
    test_heavy_score = heavy_model.predict_proba(heavy_seqs)
    assert np.allclose(gt_heavy_score, test_heavy_score)

    # Now check predict.
    test_heavy_preds = heavy_model.predict(heavy_seqs)
    gt_heavy_preds = np.argmax(gt_heavy_score, axis=0)
    assert np.allclose(gt_heavy_preds, test_heavy_preds)

    # Now check score.
    gt_heavy_score = logsumexp(gt_heavy_score, axis=0)
    test_heavy_score = heavy_model.score(heavy_seqs)
    assert np.allclose(gt_heavy_score, test_heavy_score)



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
    test_light_score = light_model.predict_proba(light_seqs)
    assert np.allclose(gt_light_score, test_light_score)

    # Now check predict.
    gt_light_preds = np.argmax(gt_light_score, axis=0)
    test_light_preds = light_model.predict(light_seqs)
    assert np.allclose(gt_light_preds, test_light_preds)

    # Now check score.
    gt_light_score = logsumexp(gt_light_score, axis=0)
    test_light_score = light_model.score(light_seqs)
    assert np.allclose(gt_light_score, test_light_score)




def test_terminal_deletions(scoring_comp_sequences):
    """Tests the categorical mixture masked scoring function by
    comparing it with a simple python implementation."""
    _, heavy_arr, _, light_arr = scoring_comp_sequences
    score_tool = SequenceScoringTool()
    dummy_heavy_arr = heavy_arr.copy()
    dummy_light_arr = light_arr.copy()
    score_tool.models["human"]["H"].em_cat_mixture_model.mask_terminal_deletions(
            heavy_arr, dummy_heavy_arr)
    score_tool.models["human"]["H"].em_cat_mixture_model.mask_terminal_deletions(
            light_arr, dummy_light_arr)

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

    assert np.allclose(dummy_heavy_arr, heavy_arr)
    assert np.allclose(dummy_light_arr, light_arr)



def test_mask_scoring(scoring_comp_sequences):
    """Ensure that the terminal deletion masking function yields
    expected results."""
    score_tool = SequenceScoringTool()
    heavy_seqs, heavy_arr, light_seqs, light_arr = scoring_comp_sequences

    heavy_masks, light_masks = [], []

    heavy_masks = [(i < 10 or i > 20) for i in range(heavy_arr.shape[1])]
    light_masks = [(i < 10 or i > 20) for i in range(light_arr.shape[1])]

    heavy_model = score_tool.models["human"]["H"]
    light_model = score_tool.models["human"]["L"]

    gt_heavy_score = np.zeros((heavy_model.get_model_specs()[0],
        heavy_arr.shape[0]))
    gt_light_score = np.zeros((light_model.get_model_specs()[0],
        light_arr.shape[0]))

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
    test_heavy_score = heavy_model.score(
            heavy_seqs, heavy_masks)

    gt_heavy_score += log_mix_weights[:,None]
    gt_heavy_score = logsumexp(gt_heavy_score, axis=0)
    assert np.allclose(gt_heavy_score, test_heavy_score)


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
    test_light_score = light_model.score(
            light_seqs, light_masks)

    gt_light_score += log_mix_weights[:,None]
    gt_light_score = logsumexp(gt_light_score, axis=0)
    assert np.allclose(gt_light_score, test_light_score)



@pytest.fixture
def scoring_comp_sequences(get_test_data_filepath, standard_aa_list):
    """Loads sequences for scoring tests, converts them to an
    MSA and encodes them for use in testing functions."""
    raw_data = pd.read_csv(os.path.join(get_test_data_filepath,
            "imgt_comp_scoring.csv.gz"))
    def prep_sequences(input_seqs):
        sca = SingleChainAnnotator()
        annotations = sca.analyze_seqs(input_seqs)
        if annotations[0][2] == "H":
            chain_type = "H"
            allowed_positions = ahip.heavy_allowed_positions
        else:
            chain_type = "L"
            allowed_positions = ahip.light_allowed_positions

        sta = SequenceTemplateAligner(allowed_positions, chain_type,
                "imgt")

        seqs = [sta.align_sequence(s, n[0], False) for s,n in zip(
            input_seqs, annotations)]

        aa_codes = {k:i for i,k in enumerate(standard_aa_list)}
        encoded_data = np.empty((len(seqs), len(seqs[0])),
                dtype=np.uint8)

        for i, seq in enumerate(seqs):
            for j, letter in enumerate(seq):
                encoded_data[i,j] = aa_codes[letter]

        return seqs, encoded_data

    heavy_seqs, heavy_arr = prep_sequences(
                raw_data[raw_data["chain_types"]=="H"].sequences.tolist()[:10])
    light_seqs, light_arr = prep_sequences(
                raw_data[raw_data["chain_types"]!="H"].sequences.tolist()[:10])
    return heavy_seqs, heavy_arr, light_seqs, light_arr
