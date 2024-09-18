"""Implements the CategoricalMixture class, for all operations involved
in generating mixture model predictions for new datapoints. The tools
necessary to fit the model (via EM) are in a separate repo.
Fitting right now is strictly via EM although we may evaluate
variational methods in future."""
import numpy as np
from ..utilities.special_functions import logsumexp
from .scoring_constants import catmix_constants as constants
from antpack.antpack_cpp_ext import getProbsCExt, getProbsCExt_masked, mask_terminal_deletions



class CategoricalMixture:
    """A CategoricalMixture model, with all the methods necessary
    to score and do inference.

    Attributes:
        mu_mix (np.ndarray): Array of type np.float64, shape (self.n_components,
            self.sequence_length, self.num_possible_items). The probability of
            each possible item for each point in sequence length for each cluster.
            Initialized to None, only converted to array once model is fitted.
        log_mu_mix (np.ndarray): Array of type np.float64, shape (self.n_components,
            self.sequence_length, self.num_possible_items). The log-probability of
            each possible item for each point in sequence length for each cluster.
            Initialized to None, only converted to array once model is fitted.
        mix_weights (np.ndarray): Array of type np.float64, shape (self.n_components).
            The weight for each distribution in the mixture.
        log_mix_weights (np.ndarray): Array of type np.float64, shape (self.n_components).
            The log of the mixture weight for each distribution in the mixture.
        n_components (int): The number of components.
        num_possible_items (int): The number of possible choices
            at each position in the sequence.
        sequence_length (int): The length of the sequences that the
            model will be fitted to / can analyze.
            
    """

    def __init__(self, n_components, num_possible_items = 21,
                sequence_length = 158):
        """Class constructor.

        Args:
            n_components (int): The number of mixture components
                (i.e. number of clusters).
            num_possible_items (int): The number of possible choices
                at each position in the sequence. For a protein sequence,
                for example, this might be the number of amino acid symbols
                that are possible; for a sequence of letters, this might
                be the number of letters in the alphabet; for shopping
                data, this might be the number of unique items the customer
                might purchase. Currently limited to the range from 1 - 255,
                this restriction will be lifted in a future version.
            sequence_length (int): The length of the sequences that the
                model will be fitted to / can analyze.

        Raises:
            ValueError: A ValueError is raised if unacceptable arguments are
                supplied.
        """
        if n_components <= 0:
            raise ValueError("n_components must be > 0.")
        if num_possible_items > 255 or num_possible_items <= 0:
            raise ValueError("Currently num_possible_items is limited to "
                    "values from 1 to 255, inclusive.")
        if sequence_length <= 0:
            raise ValueError("Sequence length must be positive.")

        self.mix_weights = None
        self.log_mix_weights = None
        self.mu_mix = None
        self.log_mu_mix = None
        self.n_components = n_components
        self.num_possible_items = num_possible_items
        self.sequence_length = sequence_length


    def load_params(self, mu_mix, mix_weights):
        """Checks input params to ensure they are compatible with selected
        settings, and if so, sets the model parameters to the inputs.

        Args:
            mu_mix (np.ndarray): Array of type np.float64, shape (self.n_components,
                self.sequence_length, self.num_possible_items). The probability of
                each possible item for each point in sequence length for each cluster.
            mix_weights (np.ndarray): Array of type np.float64, shape (self.n_components).
                The weight for each distribution in the mixture.

        Raises:
            ValueError: A ValueError is raised if the inputs are inappropriate.
        """
        if not isinstance(mu_mix, np.ndarray) or not isinstance(mix_weights, np.ndarray):
            raise ValueError("mu_mix and mix_weights should both be numpy arrays.")
        if len(mu_mix.shape) != 3 or len(mix_weights.shape) != 1:
            raise ValueError("mu_mix and mix_weights must be 3d and 1d arrays respectively.")
        if mu_mix.shape[0] != self.n_components or mu_mix.shape[1] != self.sequence_length \
                or mu_mix.shape[2] != self.num_possible_items:
            raise ValueError("mu_mix has an inappropriate shape.")
        if mix_weights.shape[0] != self.n_components:
            raise ValueError("mix_weights has an inappropriate shape.")

        if not mu_mix.flags["C_CONTIGUOUS"] or not mix_weights.flags["C_CONTIGUOUS"]:
            raise ValueError("Non-C-contiguous arrays supplied.")

        self.mix_weights = mix_weights
        self.log_mix_weights = mix_weights.copy()
        self.log_mix_weights[self.log_mix_weights < constants.MINIMUM_PROB_VAL] = constants.MINIMUM_PROB_VAL
        self.log_mix_weights = np.log(self.log_mix_weights)

        self.mu_mix = mu_mix
        self.log_mu_mix = mu_mix.copy()
        self.log_mu_mix[self.log_mu_mix < constants.MINIMUM_PROB_VAL] = constants.MINIMUM_PROB_VAL
        self.log_mu_mix[:] = np.log(self.log_mu_mix)



    def _check_input_array(self, xdata, n_threads):
        """Checks an input array to make sure that all have
        the correct format to ensure that no problems will be encountered.

        Args:
            x (np.ndarray): A numpy array with the input data.
            n_threads (int): The number of threads requested.

        Raises:
            ValueError: A ValueError is raised if unacceptable
                input data is supplied.
        """
        if not isinstance(xdata, np.ndarray):
            raise ValueError("Unexpected input supplied.")
        if np.max(xdata) >= self.num_possible_items or np.min(xdata) < 0:
            raise ValueError("Values in input data are out of range.")
        if xdata.dtype != "uint8":
            raise ValueError("Unexpected datatype for input data.")
        if len(xdata.shape) != 2:
            raise ValueError("Unexpected shape for input data.")
        if xdata.shape[1] != self.sequence_length or xdata.shape[0] < 1:
            raise ValueError("Unexpected shape for input data.")
        if not xdata.flags["C_CONTIGUOUS"]:
            raise ValueError("Input data is not C-contiguous.")

        if n_threads <= 0:
            raise ValueError("n_threads should be a positive integer.")
        if not isinstance(n_threads, int):
            raise ValueError("n_threads should be a positive integer.")


    def prep_masked_array(self, xdata, mask = None, mask_terminal_dels = False,
            mask_gaps = False):
        """Prepares a masked version of an input array if any masking options are
        supplied.

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

        Returns:
            xmasked (np.ndarray): A 2d numpy array of type np.uint8.

        Raises:
            ValueError: A ValueError is raised if bad inputs are supplied.
        """
        xmasked = xdata.copy()

        if mask is not None:
            if not isinstance(mask, np.ndarray):
                raise ValueError("Mask must be a numpy array.")
            if mask.shape[0] != xdata.shape[1]:
                raise ValueError("Mask shape must be consistent with xdata shape.")
            if len(mask.shape) != 1:
                raise ValueError("Mask must be a 1d array.")

            #There is no amino acid 21, so we use this as a convenient "please ignore"
            #indicator.
            xmasked[:,~mask] = 21

        if mask_terminal_dels:
            mask_terminal_deletions(xmasked)

        if mask_gaps:
            xmasked[xmasked==20] = 21

        return xmasked



    def predict(self, xdata, mask = None, mask_terminal_dels = False,
            mask_gaps = False, use_mixweights = True, return_raw_probs = False,
            n_threads = 1):
        """Determine the most probable cluster for each datapoint
        in a numpy array. Note that you should also check the
        overall probability of each datapoint. If a datapoint is
        very different from your training set, it will have
        very low overall probability, but this function will
        still assign it to the most likely cluster -- whichever
        that is -- by default.

        Args:
            xdata (np.ndarray): A 2d numpy array of type np.uint8.
            n_threads (int): The number of threads to use.
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
            return_raw_probs (bool): If True, return the raw probabilities
                instead of cluster assignments.

        Returns:
            preds (np.ndarray): EITHER an array of shape (xdata.shape[0])
                containing a number from 0 to self.n_components - 1
                indicating the predicted cluster for each datapoint if
                return_raw_probs is False, OR an array of shape (ncomponents,
                xdata.shape[0]) if return_raw_probs is True, indicating
                the probability of each cluster for each datapoint.

        Raises:
            ValueError: Raised if unexpected inputs are supplied.
        """
        self._check_input_array(xdata, n_threads)
        if self.log_mu_mix is None or self.log_mix_weights is None:
            raise ValueError("Model not fitted yet.")

        resp = np.zeros((self.log_mu_mix.shape[0], xdata.shape[0]))

        if mask is not None or mask_terminal_dels or mask_gaps:
            xmasked = self.prep_masked_array(xdata, mask, mask_terminal_dels, mask_gaps)
            getProbsCExt_masked(xmasked, self.log_mu_mix, resp, n_threads)

        else:
            getProbsCExt(xdata, self.log_mu_mix, resp, n_threads)

        if use_mixweights:
            resp += self.log_mix_weights[:,None]
        if return_raw_probs:
            return resp
        cluster_assignments = resp.argmax(axis=0).astype(np.uint32)
        return cluster_assignments



    def score(self, xdata, mask = None, mask_terminal_dels = False,
            mask_gaps = False, normalize_scores = False, n_threads = 1):
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
            n_threads (int): the number of threads to use.

        Returns:
            loglik (np.ndarray): A float64 array of shape (x.shape[0])
                where each element is the log-likelihood of that
                datapoint given the model.

        Raises:
            ValueError: Raised if unexpected inputs are supplied.
        """
        self._check_input_array(xdata, n_threads)
        if self.mu_mix is None or self.mix_weights is None:
            raise ValueError("Model not fitted yet.")

        resp = np.zeros((self.log_mu_mix.shape[0], xdata.shape[0]))
        if mask is not None or mask_terminal_dels or mask_gaps:
            xmasked = self.prep_masked_array(xdata, mask, mask_terminal_dels, mask_gaps)
            getProbsCExt_masked(xmasked, self.log_mu_mix, resp, n_threads)
            resp += self.log_mix_weights[:,None]
            resp = logsumexp(resp, axis=0)
            if normalize_scores:
                resp /= (xmasked < 21).sum(axis=1)

        else:
            getProbsCExt(xdata, self.log_mu_mix, resp, n_threads)
            resp += self.log_mix_weights[:,None]
            resp = logsumexp(resp, axis=0)
            if normalize_scores:
                resp /= xdata.shape[1]

        return resp
