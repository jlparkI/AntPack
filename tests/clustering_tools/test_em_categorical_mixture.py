"""Tests fitting procedures for the categorical mixture."""
import os
import random
import pytest
import numpy as np
from scipy.special import logsumexp
from antpack.clustering_tools import EMCategoricalMixture
from antpack.antpack_cpp_ext import (py_logsumexp_axis0,
        py_test_em_online, py_test_em_offline, py_random_initialize_weights)




def test_illegal_inputs(load_mab_emclust_test_data,
        load_non_mab_test_data):
    """Test that illegal inputs raise appropriate
    exceptions."""
    seqs, _ = load_non_mab_test_data
    base_model = build_default_model_non_mab_data(load_params = False)

    with pytest.raises(RuntimeError):
        base_model.predict(seqs)

    with pytest.raises(RuntimeError):
        base_model.predict_proba(seqs)

    with pytest.raises(RuntimeError):
        base_model.score(seqs)

    with pytest.raises(RuntimeError):
        base_model.get_model_parameters()

    n_components, sequence_length, max_allowed_aas = \
            base_model.get_model_specs()
    bad_mixweights = np.zeros((n_components+1))
    bad_mu_mix = np.zeros((n_components, sequence_length-2,
        max_allowed_aas))

    with pytest.raises(RuntimeError):
        base_model.load_params(bad_mu_mix, bad_mixweights)

    base_model = build_default_model_non_mab_data(
            load_params = True)

    seqs[0] += "AAA"

    with pytest.raises(RuntimeError):
        base_model.predict(seqs)

    seqs[1] = list(seqs[1])
    seqs[1][10] = "*"
    seqs[1] = ''.join(seqs[1])

    with pytest.raises(RuntimeError):
        base_model.predict(seqs)

    seqs, numbering, _, _ = load_mab_emclust_test_data

    with pytest.raises(RuntimeError):
        base_model.predict(seqs)

    cat_mix = EMCategoricalMixture(n_components=5,
            numbering=numbering + ['130'], max_threads=2,
            verbose=False)

    with pytest.raises(RuntimeError):
        cat_mix.fit(seqs)


    with pytest.raises(RuntimeError):
        cat_mix = EMCategoricalMixture(n_components=5,
            numbering=numbering, max_threads=2,
            verbose=False, region="gollum")

    cat_mix = EMCategoricalMixture(n_components=5,
            numbering=numbering, max_threads=2,
            verbose=False, region="cdr")
    seqs[0] += "A"

    with pytest.raises(RuntimeError):
        cat_mix.fit(seqs)


def test_pruning(load_non_mab_test_data):
    """Test that weight pruning after fitting does
    what we expect."""
    seqs, _ = load_non_mab_test_data
    base_model = build_default_model_non_mab_data(
            verbose=False)
    alternate_base_model = build_default_model_non_mab_data(
            verbose=False)

    base_model.fit(seqs, max_iter = 150, n_restarts=3,
            prune_after_fitting = True)
    alternate_base_model.fit(seqs, max_iter = 150, n_restarts=3,
            prune_after_fitting = False)
    mu_mix1, mixweights1 = base_model.get_model_parameters()
    mu_mix2, mixweights2 = alternate_base_model.get_model_parameters()

    idx = np.where(mixweights2>1e-16)[0]
    assert np.allclose(mu_mix2[idx,...], mu_mix1)
    assert np.allclose(mixweights2[idx], mixweights1)

    n_components, _, _ = base_model.get_model_specs()
    assert n_components==mu_mix1.shape[0]



def test_cluster_profiles(load_non_mab_test_data):
    """Test that cluster profile construction does
    what we expect."""
    random.seed(123)
    seqs, encoded_data = load_non_mab_test_data
    retained_seqs, retained_idx = [], []

    for i, seq in enumerate(seqs):
        if random.randint(0,3) == 0:
            seqs[i] = ""
            continue
        retained_seqs.append(seq)
        retained_idx.append(i)

    base_model = build_default_model_non_mab_data(
            verbose=False)

    base_model.fit(retained_seqs, max_iter = 150,
            n_restarts=3, prune_after_fitting = True)
    cprofiles = base_model.initialize_cluster_profiles()
    base_model.update_cluster_profiles(seqs, cprofiles)

    specs = base_model.get_model_specs()
    gt_cprofiles = np.zeros((specs[0], specs[1], specs[2]),
            dtype=np.int64)
    retained_idx = np.array(retained_idx)
    preds = np.zeros((len(seqs)), dtype=np.int64)
    preds[retained_idx] = base_model.predict(retained_seqs)

    for i in range(preds.shape[0]):
        if seqs[i] == "":
            continue
        for j in range(encoded_data.shape[1]):
            gt_cprofiles[preds[i], j,
                    encoded_data[i,j]] += 1

    assert np.allclose(cprofiles, gt_cprofiles)




def test_logsumexp():
    """Tests the logsumexp function."""
    rng = np.random.default_rng(123)
    r1 = rng.uniform(size=(250,1000))
    t1 = np.zeros((1000))
    py_logsumexp_axis0(r1, t1)
    assert np.allclose(t1, logsumexp(r1, axis=0))


def test_weight_initialization():
    """Ensure that weight initialization is as expected."""
    mix_weights, cluster_params = get_initial_params_non_mab_data()
    assert np.allclose(mix_weights.sum(), 1.)
    assert np.allclose(cluster_params.sum(axis=2),
        np.ones((10, 408)))


def test_em_online(load_non_mab_test_data):
    """Tests responsibility calculations conducted by the cpp
    extension for data in memory."""
    _, encoded_data = load_non_mab_test_data
    mix_weights, mu_init = get_initial_params_non_mab_data()
    gt_weights, gt_lower_bound, gt_rik_counts, \
            gt_net_resp, _ = ground_truth_em_calcs(encoded_data,
                    mix_weights, mu_init)

    lower_bound, net_resp = py_test_em_online(encoded_data, mix_weights,
            mu_init, 1)
    assert np.allclose(gt_weights, mix_weights)
    assert np.allclose(gt_rik_counts, mu_init)
    assert np.allclose(lower_bound, gt_lower_bound)
    assert np.allclose(gt_net_resp, net_resp)


def test_em_offline(tmp_path, load_non_mab_test_data):
    """Tests responsibility calculations conducted by the cpp
    extension for data on disk."""
    _, encoded_data = load_non_mab_test_data
    fpaths = [os.path.join(tmp_path, "TEST_DISCARD1.npy"),
              os.path.join(tmp_path, "TEST_DISCARD2.npy")]
    np.save(fpaths[0], encoded_data[:10,:])
    np.save(fpaths[1], encoded_data[10:,:])

    mix_weights, mu_init = get_initial_params_non_mab_data()
    gt_weights, gt_lower_bound, gt_rik_counts, \
            gt_net_resp, _ = ground_truth_em_calcs(encoded_data,
                    mix_weights, mu_init)

    lower_bound, net_resp = py_test_em_offline(fpaths, mix_weights,
            mu_init, 1)

    assert np.allclose(gt_weights, mix_weights)
    assert np.allclose(gt_rik_counts, mu_init)
    assert np.allclose(lower_bound, gt_lower_bound)
    assert np.allclose(gt_net_resp, net_resp)


def test_fitting_procedure(tmp_path, load_non_mab_test_data):
    """Runs some quick and dirty sanity checks on the
    fitting procedure for in-memory and on-disk data."""
    seqs, encoded_data = load_non_mab_test_data
    init_model = build_default_model_non_mab_data(load_params = True,
            verbose=False)
    base_model = build_default_model_non_mab_data(load_params = False,
            verbose=False)
    threaded_model = build_default_model_non_mab_data(load_params=False,
            verbose=False, max_threads=3)
    offline_threaded_model = build_default_model_non_mab_data(
            load_params=False, verbose=False, max_threads=3)

    init_bic = init_model.BIC(seqs)
    init_aic = init_model.AIC(seqs)

    base_model.fit(seqs, max_iter = 150, n_restarts=3,
            prune_after_fitting = False)
    final_bic = base_model.BIC(seqs)
    final_aic = base_model.AIC(seqs)

    assert final_bic < init_bic
    assert final_aic < init_aic

    fpaths = [os.path.join(tmp_path, "TEST_DISCARD1.npy"),
              os.path.join(tmp_path, "TEST_DISCARD2.npy")]
    np.save(fpaths[0], encoded_data[:51,:])
    np.save(fpaths[1], encoded_data[51:,:])

    base_model.fit(filepaths=fpaths, max_iter = 150,
            n_restarts=3, prune_after_fitting = False)
    final_offline_bic = base_model.BIC(filepaths=fpaths)
    final_offline_aic = base_model.AIC(filepaths=fpaths)

    assert final_bic < init_bic
    assert final_aic < init_aic
    assert np.allclose(final_bic, final_offline_bic)
    assert np.allclose(final_aic, final_offline_aic)

    threaded_model.fit(seqs, max_iter = 150, n_restarts=3,
            prune_after_fitting = False)
    final_threaded_bic = threaded_model.BIC(seqs)
    final_threaded_aic = threaded_model.AIC(seqs)
    assert np.allclose(final_threaded_bic,
        final_offline_bic)
    assert np.allclose(final_threaded_aic,
        final_offline_aic)

    offline_threaded_model.fit(filepaths=fpaths,
            max_iter = 150, n_restarts=3,
            prune_after_fitting = False)
    offline_threaded_bic = offline_threaded_model.BIC(
            filepaths=fpaths)
    offline_threaded_aic = offline_threaded_model.AIC(
            filepaths=fpaths)
    assert np.allclose(offline_threaded_bic,
        final_offline_bic)
    assert np.allclose(offline_threaded_aic,
        final_offline_aic)


    params = [model.get_model_parameters() for model in
            [base_model, threaded_model, offline_threaded_model]]
    for i, param_set1 in enumerate(params):
        for param_set2 in params[i+1:]:
            for p1, p2 in zip(param_set1, param_set2):
                assert np.allclose(p1, p2)



@pytest.mark.parametrize("region_code", ["all", "cdr", "cdr1", "fmwk", "cdr"])
def test_region_extraction(tmp_path, load_mab_emclust_test_data,
        region_code):
    """Check that encoding to a fasta file works properly and that
    extracting regions of interest works correctly."""
    seqs, position_codes, encoded_data, cdr_labels = load_mab_emclust_test_data
    with open(os.path.join(tmp_path, "TEMP_FASTA_FILE.fasta"), "w+",
            encoding="utf-8") as fhandle:
        for i, seq in enumerate(seqs):
            fhandle.write(f">{i}\n{seq}\n\n\n")

    if region_code == "all":
        cdr_idx = np.arange(encoded_data.shape[1])
    else:
        cdr_idx = np.array([i for i, token in enumerate(cdr_labels) if
            token.startswith(region_code)])
    selected_data = encoded_data[:,cdr_idx]

    offline_cat_mix = EMCategoricalMixture(n_components=10,
        numbering=position_codes, region=region_code,
        max_threads=5, verbose=False)
    online_cat_mix = EMCategoricalMixture(n_components=10,
        numbering=position_codes, region=region_code,
        max_threads=1, verbose=False)

    fpaths = offline_cat_mix.encode_fasta(
            os.path.join(tmp_path, "TEMP_FASTA_FILE.fasta"),
            tmp_path, chunk_size=100)

    np_data = np.vstack([np.load(v) for v in fpaths])
    assert np.allclose(np_data, selected_data)

    assert np.allclose(selected_data,
        online_cat_mix._prep_xdata(seqs))

    offline_cat_mix.fit(filepaths=fpaths)
    online_cat_mix.fit(seqs)

    mu_params1, mix_params1 = offline_cat_mix.get_model_parameters()
    mu_params2, mix_params2 = online_cat_mix.get_model_parameters()
    assert np.allclose(mu_params1, mu_params2)
    assert np.allclose(mix_params1, mix_params2)






def get_initial_params_non_mab_data(seq_length = 408, n_clusters = 10):
    """Gets initial parameters specialized for the non-mab
    test data."""
    mixture_weights = np.zeros((n_clusters))
    cluster_params = np.zeros((n_clusters, seq_length, 21))
    py_random_initialize_weights(mixture_weights, cluster_params,
            123)
    return mixture_weights, cluster_params


def build_default_model_non_mab_data(nclusters = 51, seq_length = 408,
        load_params = True, verbose=False, max_threads=1):
    """Builds an initial default model for the non mab data
    test set."""
    cat_mix = EMCategoricalMixture(n_components = nclusters,
            sequence_length = seq_length,
            max_threads=max_threads, verbose=verbose)
    if load_params:
        mix_weights, mu_mix = get_initial_params_non_mab_data(seq_length,
            nclusters)
        cat_mix.load_params(mu_mix, mix_weights)
    return cat_mix


def ground_truth_em_calcs(test_data, mix_weights, mu_in):
    """Generates a 'ground-truth' to compare
    against the em calculation routines."""
    log_mixweights = np.log(mix_weights.clip(min=1e-14))[:,None]
    mu_params = mu_in.copy()
    mu_params[mu_params < 1e-14] = 1e-14
    mu_params = np.log(mu_params)

    resp = np.zeros((mu_in.shape[0], test_data.shape[0]))
    lnorm = np.zeros((test_data.shape[0]))
    rik_counts = np.zeros(mu_in.shape)

    for k in range(mu_in.shape[0]):
        for i in range(test_data.shape[0]):
            resp_value = 0
            for j in range(test_data.shape[1]):
                resp_value += mu_params[k,j,test_data[i,j]]
            resp[k,i] = resp_value

    resp += log_mixweights
    lnorm[:] = logsumexp(resp, axis=0)
    with np.errstate(under="ignore"):
        resp[:] = np.exp(resp - lnorm[None,:])
    lower_bound = lnorm.sum()

    new_weights = resp.sum(axis=1)
    net_resp = new_weights.sum()

    for k in range(mu_in.shape[0]):
        for i in range(test_data.shape[0]):
            for j in range(test_data.shape[1]):
                rik_counts[k,j,test_data[i,j]] += resp[k,i]

    return new_weights, lower_bound, rik_counts, net_resp, test_data.shape[0]
