"""Wrapper for C++ code implementing a categorical mixture
fitted via EM."""
import os
import numpy as np
from antpack.antpack_cpp_ext import EMCategoricalMixtureCpp


class EMCategoricalMixture(EMCategoricalMixtureCpp):
    """A categorical mixture model that is fitted using
    the EM algorithm."""

    def __init__(self, n_components:int, sequence_length:int=0,
            numbering:list=None, numbering_scheme:str="imgt",
            cdr_scheme:str="imgt", region:str="all",
            max_threads:int=2, verbose:bool=True):
        """Constructor.

        Args:
            n_components (int): The number of clusters to use.
            sequence_length (int): The length of the input
                sequences. This can be set to 0 if a list of
                position codes is supplied as the "numbering"
                argument. If no list is supplied for the
                "numbering" argument however this must be
                an integer greater than 0.
            numbering (list): A list of position codes for
                a numbering scheme. If your sequences are an MSA
                then this numbering list would be the list of position
                codes for the MSA (use SingleChainAnnotator's build_msa
                function to see an example). This argument does not have
                to be supplied -- you can leave this as an empty list.
                If it is an empty list however you must supply a number
                for sequence length. The advantage to supplying numbering
                is that you can use it to cluster just one region of your
                input sequences (e.g. the CDRs) if desired.
            numbering_scheme (str): A valid numbering scheme; one of
                'imgt', 'aho', 'martin', 'kabat'. If you supply a numbering
                list and a region (see below) this will be used to determine
                what region should be extracted; otherwise it is ignored.
            cdr_scheme (str): A valid scheme for assigning CDRs; one of
                'imgt', 'aho', 'martin', 'kabat', 'north'. If you supply a numbering
                list and a region (see below) this will be used to determine
                what region should be extracted; otherwise it is ignored.
                This value can be different from "numbering_scheme" to allow
                you to use e.g. Kabat or North CDR definitions with IMGT
                or Martin numbering.
            region (str): If "all" the full input sequences are
                clustered. Otherwise you can supply "fmwk1", "fmwk2",
                "fmwk3" or "fmwk4" to cluster a specific framework region,
                "cdr1", "cdr2", "cdr3" to cluster a specific cdr,
                "fmwk" to cluster all framework regions or
                "cdr" to cluster all cdrs, then the model will extract
                just that region during fitting and prediction. If this
                argument is NOT "all", you must supply the 'numbering',
                'numbering_scheme' and 'cdr_scheme' arguments.
            max_threads (int): The maximum number of threads
                to use. The tool will use up to this number of threads
                wherever it makes sense to do so.
            verbose (bool): If True, print loss on every fitting
                iteration.
        """
        super().__init__(n_components, 21, sequence_length,
                max_threads, verbose)


    def get_model_parameters(self):
        """Returns the cluster parameters followed by
        the mixture weights as two numpy arrays."""
        n_components, sequence_length, n_aas = self.get_specs()
        mu_mix = np.zeros((n_components, sequence_length, n_aas))
        mix_weights = np.zeros((n_components))
        self.get_model_parameters_cpp(mu_mix, mix_weights)
        return mu_mix, mix_weights


    def get_model_specs(self):
        """Returns n_components, sequence_length and
        the maximum number of allowed aas (generally 21)."""
        return self.get_specs()


    def fit(self, sequences:list=None, filepaths:list=None,
            max_iter:int=150, tol:float=1e-3,
            n_restarts:int=3, random_state:int=123,
            prune_after_fitting:bool=True):
        """Fits the mixture model to either a list of filepaths
        representing numpy arrays saved as npy files OR an
        input numpy array.

        Args:
            sequences (list): Either None or a list of sequences that
                are all the same length (i.e. an MSA). If None,
                you must supply a list of filepaths where data is
                stored (next argument).
            filepaths (list): Either None or a list of filepaths,
                each of which is a 2d numpy array of type uint8
                containing encoded sequence data. The easiest
                way to generate these is to call `encode_fasta`
                below, which will take an input fasta file (possibly
                compressed) and convert it to numpy arrays on disk.
                This is useful when your dataset is too large to
                load to memory. If this argument is None, sequences
                cannot also be None.
            max_iter (int): The maximum number of iterations to run
                per restart.
            tol (float): If the loss changes by less than this amount,
                assume convergence and end the restart.
            n_restarts (int): The number of times to re-run fitting
                with randomly reinitialized weights.
            random_state (int): The random seed.
            prune_after_fitting (bool): If True, empty clusters are
                removed at the end of fitting. This will change
                the n_components of the model. Call self.get_model_specs()
                to get the new n_components value.
        """
        if isinstance(sequences, list):
            self.fit_online(sequences, max_iter, tol,
                    n_restarts, random_state, prune_after_fitting)
        elif isinstance(filepaths, list):
            self.fit_offline(filepaths, max_iter, tol,
                    n_restarts, random_state, prune_after_fitting)
        else:
            raise RuntimeError("sequences and filepaths cannot both "
                    "be None.")


    def BIC(self, sequences, filepaths):
        """Calculates the BIC or Bayes information criterion for
        a fitted model.

        Args:
            sequences (list): A list of sequences that
                are all the same length (i.e. an MSA). If None,
                you must supply a list of filepaths where data is
                stored (next argument).
            filepaths (list): Either None or a list of filepaths,
                each of which is a 2d numpy array of type uint8
                containing encoded sequence data. The easiest
                way to generate these is to call `encode_fasta`
                below, which will take an input fasta file (possibly
                compressed) and convert it to numpy arrays on disk.
                This is useful when your dataset is too large to
                load to memory. If this argument is None, sequences
                cannot also be None.

        Returns:
            bic (float): The Bayes information criterion for the input
                dataset.
        """
        if isinstance(sequences, list):
            return self.BIC_online(sequences)
        if isinstance(filepaths, list):
            return self.BIC_offline(filepaths)
        raise RuntimeError("sequences and filepaths cannot both be None.")


    def AIC(self, sequences, filepaths):
        """Calculates the AIC or Akaike information criterion for
        a fitted model.

        Args:
            sequences (list): A list of sequences that
                are all the same length (i.e. an MSA). If None,
                you must supply a list of filepaths where data is
                stored (next argument).
            filepaths (list): Either None or a list of filepaths,
                each of which is a 2d numpy array of type uint8
                containing encoded sequence data. The easiest
                way to generate these is to call `encode_fasta`
                below, which will take an input fasta file (possibly
                compressed) and convert it to numpy arrays on disk.
                This is useful when your dataset is too large to
                load to memory. If this argument is None, sequences
                cannot also be None.

        Returns:
            aic (float): The AIC for the input dataset.
        """
        if isinstance(sequences, list):
            return self.AIC_online(sequences)
        if isinstance(filepaths, list):
            return self.AIC_offline(filepaths)
        raise RuntimeError("sequences and filepaths cannot both be None.")



    def predict(self, sequences, mask = None, mask_terminal_dels = False,
            mask_gaps = False, use_mixweights = True):
        """Determine the most probable cluster for each datapoint
        in a numpy array. Note that you should also check the
        overall probability of each datapoint. If a datapoint is
        very different from your training set, it will have
        very low overall probability, but this function will
        still assign it to the most likely cluster -- whichever
        that is -- by default.

        Args:
            sequences (list): A list of sequences that
                are all the same length (i.e. an MSA).
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
            RuntimeError: Raised if unexpected inputs are supplied.
        """
        xdata = np.zeros((len(sequences), self.get_extract_length() ))
        self.extract_and_encode(sequences, xdata)
        cluster_assignments = np.zeros((xdata.shape[0]), dtype=np.int64)

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


    def predict_proba(self, sequences, mask = None, mask_terminal_dels = False,
            mask_gaps = False, use_mixweights = True):
        """Determine the probability of each datapoint for each cluster
        in the input array.

        Args:
            sequences (list): A list of sequences that
                are all the same length (i.e. an MSA).
            mask (np.ndarray): Either None or a numpy array of type bool
                and with shape[0] == sequence length. If not None, indicated
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
            RuntimeError: Raised if unexpected inputs are supplied.
        """
        xdata = np.zeros((len(sequences), self.get_extract_length() ))
        self.extract_and_encode(sequences, xdata)
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


    def score(self, sequences, mask = None, mask_terminal_dels = False,
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
            RuntimeError: Raised if unexpected inputs are supplied.
        """
        xdata = np.zeros((len(sequences), self.get_extract_length() ))
        self.extract_and_encode(sequences, xdata)

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



    def encode_fasta(self, fasta_filepath, temporary_dir,
            chunk_size=10000):
        """This function takes all of the sequences in the fasta_filepath
        and saves them as 2d numpy arrays with chunk_size sequences in
        each array of type uint8 under
        temporary_dir, then returns the list of filepaths. You can
        supply this list of filepaths to `fit`, `BIC` and `AIC`, then
        you can `os.remove` the filepaths when you are done fitting.
        This is a convenient way to fit the model for datasets that
        are too large to load to memory. The sequences must all
        be the same length, i.e. should be an MSA.

        Expect this to take #sequences * sequence_length / 1e9 GB
        of disk space; so if you have 50 million sequences of length
        150, this will take 7.5 GB of disk space.

        This method is preferred to loading the sequences from the fasta
        file on each iteration because it allows us to greatly accelerate
        fitting (the model will load and process up to max_threads
        numpy files at a time).

        Args:
            fasta_filepath (str): The location of the fasta file. It
                can be a gzipped file or uncompressed.
            temporary_dir (str): A path to a (probably) temporary
                directory where the encoded sequences can be saved.
            chunk_size (int): The maximum number of sequences per numpy
                file. The model will try to load up to max_threads *
                chunk_size sequences at a time during fitting, so
                do not set this too high to avoid excess memory
                consumption. For a reasonable number of threads,
                10000 is likely fine.

        Returns:
            numpy_filepaths (list): A list of absolute filepaths
                where the temporary numpy files with encoded sequences
                for fitting are stored. You can directly supply this
                list to `fit`, `bic`, `aic`. Once you are done fitting
                these files can safely be deleted.

        Raises:
            RuntimeError: An error will be raised if invalid sequences
                (all different lengths, unexpected characters etc)
                are found in the fasta file.
        """
        return []
