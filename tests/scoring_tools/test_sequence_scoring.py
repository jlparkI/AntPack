"""Tests basic functionality for the SequenceScoringTool class.
Assumes that pandas is installed."""
import os
import unittest
import pytest
import numpy as np
import pandas as pd
from antpack import SequenceScoringTool, SingleChainAnnotator
from antpack.scoring_tools.scoring_constants import allowed_imgt_pos as ahip
from antpack.scoring_tools.scoring_constants import scoring_constants as useful_constants
from antpack.utilities.model_loader_utils import load_model

from ..conftest import (get_test_data_filepath, get_test_base_filepath)


def test_get_standard_positions():
    """Ensure that the standard positions are returned
    as expected."""
    score_tool = SequenceScoringTool()
    assert (score_tool.get_standard_positions("H")==
        ahip.heavy_allowed_positions)
    assert (score_tool.get_standard_positions("L")==
        ahip.light_allowed_positions)



def test_get_standard_mask():
    """Ensures that mask construction proceeds as expected."""
    score_tool = SequenceScoringTool()
    aligner = SingleChainAnnotator()

    for chain_type in ["H", "L"]:
        std_nmbr = score_tool.get_standard_positions(chain_type)
        for nmbr_scheme in ["kabat", "imgt"]:
            gt_region_labels = aligner.assign_cdr_labels(std_nmbr,
                    chain_type, nmbr_scheme)

            for region in ["fmwk", "cdr"]:
                gt_mask = [False for reg in gt_region_labels]
                for i, token in enumerate(gt_region_labels):
                    if token.startswith(region):
                        gt_mask[i] = True

                gt_mask = np.array(gt_mask)
                assert (np.allclose(gt_mask,
                        score_tool.get_standard_mask(
                            chain_type, region, nmbr_scheme)))


def test_convert_to_array(raw_scoring_data):
    """Checks that the conversion of a sequence to
    an array is working correctly."""
    score_tool = SequenceScoringTool()
    raw_data = raw_scoring_data

    for seq in raw_data["sequences"].tolist()[:10]:
        alternate_array = score_tool.convert_sequence_to_array(seq)[1]
        assigned_chain, aligned_seq = \
                score_tool._prep_sequence(seq)
        gt_array = np.zeros((1, len(aligned_seq)), dtype=np.uint8)
        for i, letter in enumerate(aligned_seq):
            gt_array[0,i] = useful_constants.aa_list.index(letter)

        assert np.allclose(gt_array, alternate_array)



def test_calc_aa_probs_closest_cluster(raw_scoring_data):
    """Ensure the calc per aa probs function is working
    correctly and also get_closest_clusters."""
    score_tool = SequenceScoringTool()
    raw_data = raw_scoring_data

    for seq in raw_data["sequences"].tolist()[:10]:
        chain_type, seq_as_array = score_tool.convert_sequence_to_array(seq)
        closest_cluster = score_tool.get_closest_clusters(seq)[0][0]

        alternate_chain_type, test_mu, test_likely_aas = \
                score_tool.calc_per_aa_probs(seq, closest_cluster)
        assert alternate_chain_type==chain_type

        cat_model = score_tool.models["human"][chain_type]
        assigned_cluster = np.zeros((1), dtype=np.int64)
        cat_model.em_cat_mixture_model.predict_cpp(
                seq_as_array, assigned_cluster, True)
        assigned_cluster = int(assigned_cluster[0])
        assert assigned_cluster==closest_cluster

        model_mu = cat_model.get_model_parameters()[0]
        best_mu = model_mu[assigned_cluster,...]
        seq_as_array = seq_as_array.flatten()
        selected_mu = []
        gt_likely_aas = []
        for i, token in enumerate(seq_as_array.tolist()):
            if seq_as_array[i] >= 20:
                continue
            selected_mu.append(best_mu[i,token])
            gt_likely_aas.append(useful_constants.aa_list[
                np.argmax(best_mu[i,:])])

        selected_mu = np.log(np.array(selected_mu).clip(min=1e-16))
        assert test_likely_aas==gt_likely_aas
        assert np.allclose(np.array(selected_mu),
            test_mu)



def test_error_checking(raw_scoring_data):
    """Check that invalid data passed to the sequence scoring
    tool raises expected exceptions."""
    score_tool = SequenceScoringTool()
    raw_data = raw_scoring_data

    # Sequences containing unrecognized characters should raise
    # a value error for these two routines, and return np.nan
    # for the other two.
    with pytest.raises(ValueError):
        score_tool.suggest_humanizing_mutations("WoW")

    assert score_tool.get_diagnostic_info("WoW")[2]=="unknown"
    assert np.isnan(score_tool.score_seqs(["WoW"])[0])

    aligned_seq = score_tool._prep_sequence(
            raw_data["sequences"].tolist()[0])[1]

    proba_test = score_tool.models["human"]["H"].predict_proba(
            [aligned_seq])
    assert len(proba_test.shape) == 2
    assert proba_test.shape[1] == 1
    assert (proba_test.shape[0] ==
            score_tool.models["human"]["H"].get_model_specs()[0])


def test_sequence_extraction(raw_scoring_data):
    """Make sure the sequence prep feature is providing expected
    results."""
    score_tool = SequenceScoringTool()
    raw_data = raw_scoring_data
    aligners = {"H":SingleChainAnnotator(["H"]),
        "K":SingleChainAnnotator(["K"]),
        "L":SingleChainAnnotator(["L"])  }

    chain_position_dicts = {"H":
            {k:i for i, k in enumerate(ahip.heavy_allowed_positions)},
            "L":
            {k:i for i, k in enumerate(ahip.light_allowed_positions)}
            }

    for seq, chain_name in zip(raw_data["sequences"].tolist()[:10],
            raw_data["chain_types"].tolist()[:10]):
        assigned_chain, aligned_seq = score_tool._prep_sequence(seq)
        if chain_name == "K":
            assert assigned_chain == "L"
        else:
            assert assigned_chain == chain_name

        position_dict = chain_position_dicts[assigned_chain]

        numbering = aligners[assigned_chain].analyze_seq(seq)[0]
        seq_extract = ["-" for i in range(len(position_dict))]

        for position, aa in zip(numbering, seq):
            if position == "-" or position not in position_dict:
                continue
            seq_extract[position_dict[position]] = aa
        seq_extract = ''.join(seq_extract)
        assert seq_extract == aligned_seq




def test_scoring_consistency(raw_scoring_data, get_test_base_filepath):
    """Test that scores assigned by different procedures give what
    we expect."""
    score_tool = SequenceScoringTool(offer_classifier_option=True)
    adj_score_tool = SequenceScoringTool(offer_classifier_option=False,
            normalization="training_set_adjust")
    # This is an ugly hack, but. In prior versions of antpack, weights
    # were clipped and not allowed to go below 1e-14 for numerical
    # stability. This is higher than necessary however and 1e-16
    # (current value) works fine. For compatibility with the older
    # version against which the test was designed, we load weights
    # with the 1e-14 lower bound.
    # TODO: Update these with scores calculated without the old lower
    # bound and get rid of this ugly hack.
    project_dir = os.path.join(get_test_base_filepath, "src", "antpack",
            "scoring_tools")
    score_tool.models = {"human":{"H":load_model(project_dir, "heavy",
                                max_threads=2,
                                clip_mu_at_old_value=True),
                    "L":load_model(project_dir, "light",
                                max_threads=2,
                                clip_mu_at_old_value=True)} }
    adj_score_tool.models = {"human":{"H":load_model(project_dir, "heavy",
                                max_threads=2,
                                clip_mu_at_old_value=True),
                    "L":load_model(project_dir, "light",
                                max_threads=2,
                                clip_mu_at_old_value=True)} }

    raw_data = raw_scoring_data
    heavy_chains = raw_data[raw_data["chain_types"]=="H"]
    light_chains = raw_data[raw_data["chain_types"].isin(["K", "L"])]

    heavy_seqs = [score_tool._prep_sequence(s)[1] for s
            in heavy_chains["sequences"].tolist()]
    light_seqs = [score_tool._prep_sequence(s)[1] for s
            in light_chains["sequences"].tolist()]

    # Generate scores using a variety of different procedures. We will compare
    # these to make sure everything makes sense.

    raw_heavy_score = score_tool.models["human"]["H"].score(
            heavy_seqs)
    raw_light_score = score_tool.models["human"]["L"].score(
            light_seqs)

    unadj_heavy_score = score_tool.score_seqs(heavy_chains["sequences"].tolist())
    unadj_light_score = score_tool.score_seqs(light_chains["sequences"].tolist())
    assert np.allclose(unadj_heavy_score, raw_heavy_score)
    assert np.allclose(unadj_light_score, raw_light_score)

    og_term_del_scores = raw_data["term_del_scores"].values
    term_del_scores = adj_score_tool.score_seqs(raw_data["sequences"].tolist(),
            mask_terminal_dels=True)
    assert np.allclose(term_del_scores, og_term_del_scores)

    og_scores = raw_data["batched_scores"].values
    batched_scores = adj_score_tool.score_seqs(raw_data["sequences"].tolist())
    assert np.allclose(batched_scores, og_scores)



@pytest.fixture
def raw_scoring_data(get_test_data_filepath):
    """Loads the raw scoring data used by all tests."""
    return pd.read_csv(os.path.join(get_test_data_filepath,
        "imgt_comp_scoring.csv.gz"))
