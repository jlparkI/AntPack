"""Wrapper for C++ code implementing a categorical mixture
fitted via EM."""
import os
import numpy as np
from antpack.antpack_cpp_ext import EMCategoricalMixtureCpp


class EMCategoricalMixture(EMCategoricalMixtureCpp):
    """A categorical mixture model that is fitted using
    the EM algorithm."""

    def __init__(self, n_components, sequence_length,
            max_threads):
        """Constructor.

        Args:
            n_components (int): The number of clusters to use.
            sequence_length (int): The length of the input
                sequences.
            max_threads (int): The maximum number of threads
                to use.
        """
        super().__init__(n_components, 21, sequence_length,
                max_threads)


    def get_model_parameters(self):
        """Returns the cluster parameters followed by
        the mixture weights as two numpy arrays."""
        n_components, sequence_length, n_aas = self.get_specs()
        mu_mix = np.zeros((n_components, sequence_length, n_aas))
        mix_weights = np.zeros((n_components))
        self.get_model_parameters_cpp(mu_mix, mix_weights)
        return mu_mix, mix_weights


    def predict(self, xdata, mask = None, mask_terminal_dels = False,
            mask_gaps = False, use_mixweights = True):
        """Determine the most probable cluster for each datapoint
        in a numpy array. Note that you should also check the
        overall probability of each datapoint. If a datapoint is
        very different from your training set, it will have
        very low overall probability, but this function will
        still assign it to the most likely cluster -- whichever
        that is -- by default.

        Args:
            xdata (np.ndarray): A 2d numpy array of type np.uint8.
            mask (np.ndarray): Either None or a numpy array of type bool
                and shape (xdata.shape[1]). If not None, indicated
                positions are masked, i.e. are not taken into account
                when calculating the score.
            mask_gaps (bool): If True, all non-filled IMGT positions in the sequence
                are ignored when calculating the score. This is useful when your
                sequence has unusual deletions and you would like to ignore these.
            mask_terminal_dels (bool): If True, ignore N- and C-terminal
                deletions when assigning to a cluster.
            use_mixweights (bool): If True, take mixture weights into
                account; otherwise, find the closest cluster (even if
                it is a low-probability cluster).

        Returns:
            preds (np.ndarray): An array of shape (xdata.shape[0])
                containing a number from 0 to self.n_components - 1
                indicating the predicted cluster for each datapoint.

        Raises:
            ValueError: Raised if unexpected inputs are supplied.
        """
        cluster_assignments = np.zeros((xdata.shape[0]))

        if mask is not None or mask_terminal_dels or mask_gaps:
            xmasked = xdata.copy()
            if mask is not None:
                self.apply_mask_to_input_data(mask, xmasked)
            if mask_terminal_dels:
                self.mask_terminal_deletions(xdata, xmasked)
            if mask_gaps:
                self.mask_gaps(xdata, xmasked)
            self.predict_cpp(xmasked, cluster_assignments,
                    use_mixweights, True)

        else:
            self.predict_cpp(xdata, cluster_assignments, use_mixweights)

        return cluster_assignments


    def predict_proba(self, xdata, mask = None, mask_terminal_dels = False,
            mask_gaps = False, use_mixweights = True):
        """Determine the probability of each datapoint for each cluster
        in the input array.

        Args:
            xdata (np.ndarray): A 2d numpy array of type np.uint8
                containing encoded sequence data.
            mask (np.ndarray): Either None or a numpy array of type bool
                and shape (xdata.shape[1]). If not None, indicated
                positions are masked, i.e. are not taken into account
                when calculating the score.
            mask_gaps (bool): If True, all non-filled IMGT positions in the sequence
                are ignored when calculating the score. This is useful when your
                sequence has unusual deletions and you would like to ignore these.
            mask_terminal_dels (bool): If True, ignore N- and C-terminal
                deletions when assigning to a cluster.
            use_mixweights (bool): If True, take mixture weights into
                account; otherwise, find the closest cluster (even if
                it is a low-probability cluster).

        Returns:
            proba (np.ndarray): An array of shape (ncomponents,
                xdata.shape[0]) indicating the probability of each
                cluster for each datapoint.

        Raises:
            ValueError: Raised if unexpected inputs are supplied.
        """
        proba = np.zeros((self.get_specs()[0], xdata.shape[0]))

        if mask is not None or mask_terminal_dels or mask_gaps:
            xmasked = xdata.copy()
            if mask is not None:
                self.apply_mask_to_input_data(mask, xmasked)
            if mask_terminal_dels:
                self.mask_terminal_deletions(xdata, xmasked)
            if mask_gaps:
                self.mask_gaps(xdata, xmasked)
            self.predict_proba_cpp(xmasked, proba,
                    use_mixweights, True)

        else:
            self.predict_proba_cpp(xdata, proba, use_mixweights)

        return proba


    def score(self, xdata, mask = None, mask_terminal_dels = False,
            mask_gaps = False, normalize_scores = False):
        """Generate the overall log-likelihood of individual datapoints.
        This is very useful to determine if a new datapoint is very
        different from the training set. If the log-likelihood of
        the new datapoints is much lower than the training set
        distribution of log-likelihoods, it is fairly unlikely
        that the new datapoint is a sample from the distribution
        represented by the model, and you should not try to
        assign it to a cluster.

        Args:
            xdata (np.ndarray): An array with the input data,
                of type np.uint8.
            mask (np.ndarray): Either None or a numpy array of type bool
                and shape (xdata.shape[1]). If not None, indicated
                positions are masked, i.e. are not taken into account
                when calculating the score.
            mask_terminal_dels (bool): If True, ignore N- and C-terminal
                deletions when calculating a score.
            mask_gaps (bool): If True, all non-filled IMGT positions in the sequence
                are ignored when calculating the score. This is useful when your
                sequence has unusual deletions and you would like to ignore these.
            normalize_scores (bool): If True, normalize the score by dividing by
                the number of non-masked residues in the input.

        Returns:
            loglik (np.ndarray): A float64 array of shape (x.shape[0])
                where each element is the log-likelihood of that
                datapoint given the model.

        Raises:
            ValueError: Raised if unexpected inputs are supplied.
        """
        loglik = np.zeros((xdata.shape[0]))
        if mask is not None or mask_terminal_dels or mask_gaps:
            xmasked = xdata.copy()
            if mask is not None:
                self.apply_mask_to_input_data(mask, xmasked)
            if mask_terminal_dels:
                self.mask_terminal_deletions(xdata, xmasked)
            if mask_gaps:
                self.mask_gaps(xdata, xmasked)
            self.score_cpp(xmasked, loglik,
                    normalize_scores, True)

        else:
            self.score_cpp(xdata, loglik, normalize_scores)

        return loglik
