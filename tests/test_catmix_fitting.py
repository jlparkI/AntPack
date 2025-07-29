"""Tests fitting procedures for the categorical mixture."""
import os
import unittest
import numpy as np
from scipy.special import logsumexp
from antpack.clustering_tools import EMCategoricalMixture
from antpack.antpack_cpp_ext import (py_logsumexp_axis0,
        py_test_em_online, py_test_em_offline, py_random_initialize_weights)


class TestCatmixFitting(unittest.TestCase):


    def test_illegal_inputs(self):
        """Test that illegal inputs raise appropriate
        exceptions."""
        xdata, _ = load_test_data()
        base_model = build_default_model(load_params = False)

        with self.assertRaises(RuntimeError):
            base_model.predict(xdata)

        with self.assertRaises(RuntimeError):
            base_model.predict_proba(xdata)

        with self.assertRaises(RuntimeError):
            base_model.score(xdata)

        with self.assertRaises(RuntimeError):
            base_model.get_model_parameters()

        n_components, sequence_length, max_allowed_aas = \
                base_model.get_specs()
        bad_mixweights = np.zeros((n_components+1))
        bad_mu_mix = np.zeros((n_components, sequence_length-2,
            max_allowed_aas))

        with self.assertRaises(RuntimeError):
            base_model.load_params(bad_mu_mix, bad_mixweights)

        base_model = build_default_model(load_params = True)

        with self.assertRaises(RuntimeError):
            base_model.predict(np.ascontiguousarray(xdata[:,:-1]))

        with self.assertRaises(RuntimeError):
            base_model.score(np.ascontiguousarray(xdata[:,:-1]))


    def test_pruning(self):
        """Test that weight pruning after fitting does
        what we expect."""
        xdata, _ = load_test_data()
        base_model = build_default_model(verbose=False)
        alternate_base_model = build_default_model(verbose=False)

        base_model.fit(xdata, max_iter = 150, n_restarts=3,
                prune_after_fitting = True)
        alternate_base_model.fit(xdata, max_iter = 150, n_restarts=3,
                prune_after_fitting = False)
        mu_mix1, mixweights1 = base_model.get_model_parameters()
        mu_mix2, mixweights2 = alternate_base_model.get_model_parameters()

        idx = np.where(mixweights2>1e-16)[0]
        self.assertTrue(np.allclose(mu_mix2[idx,...], mu_mix1))
        self.assertTrue(np.allclose(mixweights2[idx], mixweights1))

        n_components, _, _ = base_model.get_specs()
        self.assertTrue(n_components==mu_mix1.shape[0])




    def test_logsumexp(self):
        """Tests the logsumexp function."""
        rng = np.random.default_rng(123)
        r1 = rng.uniform(size=(250,1000))
        t1 = np.zeros((1000))
        py_logsumexp_axis0(r1, t1)
        self.assertTrue(np.allclose(t1, logsumexp(r1, axis=0)))


    def test_weight_initialization(self):
        """Ensure that weight initialization is as expected."""
        mix_weights, cluster_params = get_initial_params()
        self.assertTrue(np.allclose(mix_weights.sum(), 1.))
        self.assertTrue(np.allclose(cluster_params.sum(axis=2),
            np.ones((10, 408)))
            )


    def test_em_online(self):
        """Tests responsibility calculations conducted by the cpp
        extension for data in memory."""
        xdata, _ = load_test_data()
        mix_weights, mu_init = get_initial_params()
        gt_weights, gt_lower_bound, gt_rik_counts, \
                gt_net_resp, _ = ground_truth_em_calcs(xdata, mix_weights, mu_init)

        lower_bound, net_resp = py_test_em_online(xdata, mix_weights,
                mu_init, 1)
        self.assertTrue(np.allclose(gt_weights, mix_weights))
        self.assertTrue(np.allclose(gt_rik_counts, mu_init))
        self.assertTrue(np.allclose(lower_bound, gt_lower_bound))
        self.assertTrue(np.allclose(gt_net_resp, net_resp))


    def test_em_offline(self):
        """Tests responsibility calculations conducted by the cpp
        extension for data on disk."""
        xdata, xfiles = load_test_data()
        mix_weights, mu_init = get_initial_params()
        gt_weights, gt_lower_bound, gt_rik_counts, \
                gt_net_resp, _ = ground_truth_em_calcs(xdata, mix_weights, mu_init)

        lower_bound, net_resp = py_test_em_offline(xfiles, mix_weights,
                mu_init, 1)
        self.assertTrue(np.allclose(gt_weights, mix_weights))
        self.assertTrue(np.allclose(gt_rik_counts, mu_init))
        self.assertTrue(np.allclose(lower_bound, gt_lower_bound))
        self.assertTrue(np.allclose(gt_net_resp, net_resp))


    def test_fitting_procedure(self):
        """Runs some quick and dirty sanity checks on the
        fitting procedure for in-memory and on-disk data."""
        xdata, xfiles = load_test_data()
        init_model = build_default_model(load_params = True, verbose=True)
        base_model = build_default_model(load_params = False, verbose=True)
        threaded_model = build_default_model(load_params=False,
                verbose=True, max_threads=3)
        offline_threaded_model = build_default_model(load_params=False,
                verbose=True, max_threads=3)

        init_bic = init_model.BIC(xdata)
        init_aic = init_model.AIC(xdata)

        base_model.fit(xdata, max_iter = 150, n_restarts=3,
                prune_after_fitting = False)
        final_bic = base_model.BIC(xdata)
        final_aic = base_model.AIC(xdata)

        self.assertTrue(final_bic < init_bic)
        self.assertTrue(final_aic < init_aic)

        base_model.fit(xfiles, max_iter = 150, n_restarts=3,
                prune_after_fitting = False)
        final_offline_bic = base_model.BIC(xfiles)
        final_offline_aic = base_model.AIC(xfiles)


        self.assertTrue(final_bic < init_bic)
        self.assertTrue(final_aic < init_aic)
        self.assertTrue(np.allclose(final_bic, final_offline_bic))
        self.assertTrue(np.allclose(final_aic, final_offline_aic))

        threaded_model.fit(xdata, max_iter = 150, n_restarts=3,
                prune_after_fitting = False)
        final_threaded_bic = threaded_model.BIC(xfiles)
        final_threaded_aic = threaded_model.AIC(xfiles)
        self.assertTrue(np.allclose(final_threaded_bic,
            final_offline_bic))
        self.assertTrue(np.allclose(final_threaded_aic,
            final_offline_aic))

        offline_threaded_model.fit(xfiles, max_iter = 150, n_restarts=3,
                prune_after_fitting = False)
        offline_threaded_bic = offline_threaded_model.BIC(xfiles)
        offline_threaded_aic = offline_threaded_model.AIC(xfiles)
        self.assertTrue(np.allclose(offline_threaded_bic,
            final_offline_bic))
        self.assertTrue(np.allclose(offline_threaded_aic,
            final_offline_aic))


        params = [model.get_model_parameters() for model in
                [base_model, threaded_model, offline_threaded_model]]
        for i, param_set1 in enumerate(params):
            for param_set2 in params[i+1:]:
                for p1, p2 in zip(param_set1, param_set2):
                    self.assertTrue(np.allclose(p1, p2))




def load_test_data():
    """Loads some saved test data."""
    current_dir = os.getcwd()
    start_dir = os.path.abspath(os.path.dirname(__file__))
    data_path = os.path.join(start_dir, "test_data",
            "mixture_model_test_data", "encoded_test_data.npy")
    xdata = np.load(data_path)
    xdata = np.vstack([xdata, xdata])
    xfiles = 2 * [os.path.abspath(data_path)]
    os.chdir(current_dir)
    return xdata, xfiles


def get_initial_params(seq_length = 408, n_clusters = 10):
    """Gets initial parameters specialized for the test set."""
    mixture_weights = np.zeros((n_clusters))
    cluster_params = np.zeros((n_clusters, seq_length, 21))
    py_random_initialize_weights(mixture_weights, cluster_params,
            123)
    return mixture_weights, cluster_params


def build_default_model(nclusters = 51, seq_length = 408,
        load_params = True, verbose=False, max_threads=1):
    """Builds an initial default model for the test set."""
    cat_mix = EMCategoricalMixture(n_components = nclusters,
            sequence_length = seq_length,
            max_threads=max_threads, verbose=verbose)
    if load_params:
        mix_weights, mu_mix = get_initial_params(seq_length,
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






if __name__ == "__main__":
    unittest.main()
