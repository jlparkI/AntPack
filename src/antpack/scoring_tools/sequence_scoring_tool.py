"""Provides a tool for scoring individual sequences."""
import os
import numpy as np
from ..numbering_tools import SingleChainAnnotator
from .scoring_constants import scoring_constants as constants
from .scoring_constants import allowed_imgt_pos as ahip
from ..utilities.model_loader_utils import load_model
from ..utilities.data_compile_utils import get_position_dict, get_reverse_position_dict


class SequenceScoringTool():
    """Tool for scoring sequences."""

    def __init__(self, offer_classifier_option = False,
            normalization = "none"):
        """Class constructor.

        Args:
            offer_classifier_option (bool): If True, the object loads
                additional species-specific models (mouse, rhesus monkey,
                rat) and can use these to function as a classifier. If this
                is False, setting mode='classifier' when calling a scoring
                function will result in an error.
            normalization (str): One of "none", "training_set_adjust" and
                "normalize". If normalize, the score is divided by the
                number of non-masked residues. If "training_set_adjust",
                the median training set score is subtracted from the score.
                "none" is fine unless you want to combine scores for two
                different chains, compare scores across chains, or compare
                scores for different regions, in which case one of the
                other two options is recommended.
        """

        project_dir = os.path.abspath(os.path.dirname(__file__))

        self.position_dict = {"H":get_position_dict("heavy")[0],
                "L":get_position_dict("light")[0]}
        self.rev_position_dict = {"H":get_reverse_position_dict("heavy"),
                "L":get_reverse_position_dict("light")}

        self.chain_map = {"H":"H", "K":"L", "L":"L", "":"unknown"}

        if offer_classifier_option:
            self.models = {
                    "human":{"H":load_model(project_dir, "heavy"),
                        "L":load_model(project_dir, "light")},
                    "mouse":{"H":load_model(project_dir, "heavy", species="mouse"),
                        "L":load_model(project_dir, "light", species="mouse")},
                    "rhesus":{"H":load_model(project_dir, "heavy", species="rhesus"),
                        "L":load_model(project_dir, "light", species="rhesus")},
                    "rat":{"H":load_model(project_dir, "heavy", species="rat") }
                }
        else:
            self.models = {"human":{"H":load_model(project_dir, "heavy"),
                        "L":load_model(project_dir, "light")} }

        self.aa_list = constants.aa_list
        self.aa_dict = {aa:i for (i, aa) in enumerate(self.aa_list)}
        self.aligner = SingleChainAnnotator(chains=["H", "K", "L"],
                scheme = "imgt")

        self.normalize_scores = False
        self.score_adjustments = {"H":0.0, "L":0.0}
        if normalization == "training_set_adjust":
            self.score_adjustments = {"H":constants.HEAVY_MEDIAN_SCORE,
                    "L":constants.LIGHT_MEDIAN_SCORE}
        elif normalization == "normalize":
            self.normalize_scores = True



    def _simple_prep_sequence(self, seq:str):
        """Determines the chain type of a sequence and aligns it
        using the IMGT numbering. Generates / returns less info
        than full_prep_sequence, so useful if only standard
        scoring is desired.

        Args:
            seq (str): An input sequence

        Returns:
            chain_name (str): One of "L", "H". Note that "K" and "L"
                are combined since the same mixture model is used for
                "K" and "L".
            seq_array (np.ndarray): A numpy array of shape (1, L)
                where L is the total number of recognized IMGT positions
                of type uint8. This is the encoded version of seq extract.
        """
        numbering, _, chain_type, _ = self.aligner.analyze_seq(seq)
        chain_type = self.chain_map[chain_type]
        if chain_type == "unknown":
            output_seq = ["-" for a in seq]
            seq_arr = np.array([[self.aa_dict[letter] for letter
                    in output_seq]], dtype=np.uint8)
            return chain_type, seq_arr

        position_dict = self.position_dict[chain_type]
        output_seq = ["-" for i in range(len(position_dict))]

        for position, aa in zip(numbering, seq):
            if position == "-":
                continue
            if position not in position_dict:
                continue
            output_seq[position_dict[position]] = aa

        seq_arr = np.array([[self.aa_dict[letter] for letter
                in output_seq]], dtype=np.uint8)
        return chain_type, seq_arr




    def _full_prep_sequence(self, seq:str):
        """Determines the chain type of a sequence and aligns it
        using the IMGT numbering. Also generates a back mapping
        from array position to input that is useful for more
        complicated procedures (e.g. humanization) and lists of
        unusual positions / insertions that are useful for
        troubleshooting.

        Args:
            seq (str): An input sequence

        Returns:
            output_seq (str): The extracted and aligned sequence.
            chain_name (str): One of "L", "H". Note that "K" and "L"
                are combined since the same mixture model is used for
                "K" and "L".
            seq_array (np.ndarray): A numpy array of shape (1, L)
                where L is the total number of recognized IMGT positions
                of type uint8. This is the encoded version of seq extract.
            bad_gaps (list): A list of unexpected gaps at sites where
                a deletion is unusual in the IMGT numbering system.
            bad_positions (list): A list of unexpected insertions at
                sites where an insertion is not expected in IMGT.
            backmap (dict): A dictionary mapping from the position in the
                numbered sequence back to the original sequence.
        """
        numbering, _, chain_type, _ = self.aligner.analyze_seq(seq)
        chain_type = self.chain_map[chain_type]
        if chain_type == "unknown":
            output_seq = ["-" for a in seq]
            seq_arr = np.array([[self.aa_dict[letter] for letter
                    in output_seq]], dtype=np.uint8)
            return output_seq, chain_type, seq_arr, [], [], {}

        position_dict = self.position_dict[chain_type]
        rev_dict = self.rev_position_dict[chain_type]
        output_seq = ["-" for i in range(len(position_dict))]
        bad_positions, backmap, input_counter = [], {}, 0

        for position, aa in zip(numbering, seq):
            if position == "-":
                input_counter += 1
                continue
            if position not in position_dict:
                bad_positions.append(position)
                input_counter += 1
                continue
            output_seq[position_dict[position]] = aa
            if aa != "-":
                backmap[position_dict[position]] = input_counter
                input_counter += 1

        bad_gaps = [rev_dict[position] for position, aa in
                enumerate(output_seq) if aa == "-"]

        seq_arr = np.array([[self.aa_dict[letter] for letter
                in output_seq]], dtype=np.uint8)
        return output_seq, chain_type, seq_arr, bad_gaps, bad_positions, backmap


    def _build_mask_arr(self, mask_positions:list, mask_cdr3:bool,
            chain_type:str):
        """Constructs a mask array that can be passed to the wrapped
        categorical mix tool using a user-specified list of mask
        positions and/or a request to mask cdr3.

        Args:
            mask_positions (list): Either None or a (possibly empty) list of
                strings defining IMGT positions to be masked.
            mask_cdr3 (bool): Indicates whether cdr3 should be masked in
                addition to the user-specified mask (if any).
            chain_type (str): One of "H", "L".

        Returns:
            mask_arr (np.ndarray): A numpy array of type bool.
        """
        if chain_type not in ("H", "L"):
            raise ValueError("Unrecognized chain type passed to an internal function.")
        position_dict = self.position_dict[chain_type]
        scoring_mask = [True for i in range(len(position_dict))]

        if mask_positions is not None:
            for position in mask_positions:
                scoring_mask[position_dict[position]] = False

        if mask_cdr3:
            for position in ahip.imgt_cdrs[chain_type]["3"]:
                scoring_mask[position_dict[position]] = False

        return np.array(scoring_mask, dtype=bool)




    def get_standard_positions(self, chain_type):
        """Returns a list of the standard positions used by SAM when scoring a
        sequence. IMGT numbered positions outside this set are ignored when
        assigning a score. If you want to use a generic mask for an IMGT-
        defined region, call get_standard_mask. If you want to create a custom
        mask for everything except a specific region of interest, use
        this function to get a list of all positions."""
        if chain_type == "H":
            return ahip.heavy_allowed_positions
        if chain_type in ("K", "L"):
            return ahip.light_allowed_positions
        raise ValueError("Unrecognized chain type supplied.")



    def get_standard_mask(self, chain_type:str, region:str = "framework_1"):
        """Returns a mask for ALL positions EXCEPT a specified IMGT-
        defined region. You can then use this mask as input to
        score_seqs to see what the score would be if only that region
        were included.

        Args:
            chain_type (str): One of "H", "L".
            region (str): One of 'framework_1', 'framework_2', 'framework_3',
                'framework_4', 'cdr_1', 'cdr_2', 'cdr_3', 'cdr_4'. This
                function will construct a mask that excludes all other regions.

        Returns:
            mask (list): A list of excluded imgt positions that can be passed to
                one of the scoring functions.

        Raises:
            ValueError: A ValueError is raised if unexpected inputs are supplied.
        """
        if chain_type not in ["H", "L"]:
            raise ValueError("chain_type must be one of 'H', 'L'.")
        if region not in ("framework_1", "framework_2", "framework_3", "framework_4",
                "cdr_1", "cdr_2", "cdr_3"):
            raise ValueError("Unexpected region supplied.")

        if region.startswith("framework"):
            valid_pos = set(ahip.imgt_framework_regions[chain_type][region.split("_")[1]])
        else:
            valid_pos = set(ahip.imgt_cdrs[chain_type][region.split("_")[1]])

        mask = [p for p in self.position_dict[chain_type].keys() if
                p not in valid_pos]
        return mask




    def get_diagnostic_info(self, seq):
        """Gets diagnostic information for an input sequence that is useful
        if troubleshooting.

        Args:
            seq (str): The sequence. May be either heavy chain or
                light.

        Returns:
            gapped_imgt_positions (list): A list of gaps at sites that are
                sometimes filled in the IMGT numbering system. The presence
                of these is not a cause for concern but may be useful for
                diagnostic purposes.
            unusual_positions (list): A list of unexpected insertions at
                sites where an insertion is very unusual in the training set.
                Take note of these if any are found -- these are not taken into
                account when scoring the sequence, so if this is not an empty
                list, the score may be less reliable.
            chain_name (str): one of "H", "L" for heavy or light. Indicates the
                chain type to which the sequence was aligned.
        """
        _, chain_name, _, bad_gaps, bad_positions, _ = self._full_prep_sequence(seq)
        if chain_name not in ["H", "L"]:
            return bad_gaps, bad_positions, chain_name

        return bad_gaps, bad_positions, chain_name



    def score_seqs(self, seq_list, mask_cdr3:bool = False,
            custom_light_mask:list = None, custom_heavy_mask:list = None,
            mask_terminal_dels:bool = False, mask_gaps:bool = False,
            mode:str = "score"):
        """Scores a list of sequences in batches or assigns them
        to clusters. Can be used in conjunction with a user-supplied
        mask (for positions to ignore) and in conjunction with Substantially faster than
        single seq scoring but does not offer the option to retrieve
        diagnostic infoCan also be used to assign a large number of
        sequences to clusters as well.

        Args:
            seq_list (str): The list of input sequences. May contain both
                heavy and light.
            mask_cdr3 (bool): If True, ignore IMGT-defined CDR3 when assigning a
                score. CDR3 is not distinctive across species so this is often
                useful. Ignored if mode is 'assign', 'assign_no_weights'.
            custom_light_mask (list): Either None or a list of strings indicating
                IMGT positions to ignore. Use self.get_standard_positions and/or
                self.get_standard_mask to construct a mask. This can be useful
                if you just want to score a specific region, or if there is a large
                deletion that should be ignored.
            custom_heavy_mask (list): Either None or a list of strings indicating
                IMGT positions to ignore. Use self.get_standard_positions and/or
                self.get_standard_mask to construct a mask. This can be useful
                if you just want to score a specific region, or if there is a
                large deletion that should be ignored.
            mask_terminal_dels (bool): If True, N and C-terminal deletions are
                masked when calculating a score or assigning to a cluster. Useful
                if there are large unusual deletions at either end of the sequence that
                you would like to ignore when scoring.
            mask_gaps (bool): If True, all non-filled IMGT positions in the sequence
                are ignored when calculating the score. This is useful when your
                sequence has unusual deletions and you would like to ignore these.
            mode (str): One of 'score', 'assign', 'assign_no_weights', 'classifier'.
                If score, returns the human generative model score.
                If 'assign', provides the most likely cluster number
                for each input sequence. If 'assign_no_weights',
                assigns the closest cluster ignoring mixture weights,
                so that the closest cluster is assigned even if that
                cluster is a low-probability one. If 'classifier',
                assigns a score using the Bayes' rule classifier, which
                also takes into account some info regarding other species.
                'classifier' is not a good way to score sequences in general
                because it only works well for sequences of known origin,
                so it should only be used for testing.

        Returns:
            output_scores (np.ndarray): log( p(x) ) for all input sequences.
        """

        light_arr, heavy_arr = [], []
        light_idx, heavy_idx = [], []
        output_scores = np.zeros((len(seq_list)))

        heavy_mask, light_mask = None, None
        if custom_heavy_mask is not None or mask_cdr3:
            heavy_mask = self._build_mask_arr(custom_heavy_mask, mask_cdr3, "H")
        if custom_light_mask is not None or mask_cdr3:
            light_mask = self._build_mask_arr(custom_light_mask, mask_cdr3, "L")

        for i, seq in enumerate(seq_list):
            chain_name, arr = self._simple_prep_sequence(seq)

            if chain_name == "L":
                light_arr.append(arr)
                light_idx.append(i)
            elif chain_name == "H":
                heavy_arr.append(arr)
                heavy_idx.append(i)
            else:
                output_scores[i] = np.nan

        if len(heavy_arr) > 0:
            self._batch_score(heavy_arr, output_scores,
                    heavy_idx, "H", mode, heavy_mask,
                    mask_terminal_dels, mask_gaps)

        if len(light_arr) > 0:
            self._batch_score(light_arr, output_scores,
                    light_idx, "L", mode, light_mask,
                    mask_terminal_dels, mask_gaps)

        return output_scores



    def get_closest_clusters(self, seq:str, nclusters:int = 1, return_ids_only = True):
        """Gets the closest cluster(s) for a given sequence. These can be used
        for humanization, to determine which amino acids in the input sequence
        are most problematic, or to generate new sequences containing motifs of interest.

        Args:
            seq (str): The input sequence.
            nclusters (int): The number of clusters to retrieve.
            return_ids_only (bool): If True, only the cluster numbers are returned,
                not the actual clusters.

        Returns:
            cluster_idx (np.ndarray): The index number for each cluster in the
                original model that is returned. Sorted from highest probability
                cluster to lowest probability.
            chain_name (str): One of "H" or "L", indicating whether chain is
                heavy or light (K and L are both mapped to L).
            mu_mix (np.ndarray): An array of shape (nclusters, sequence_length,
                21), where 21 is the number of possible AAs. The clusters
                are sorted in order from most to least likely given the
                input sequence. Not returned if return_ids_only is True.
            mixweights (np.ndarray): An array of shape (nclusters) containing
                the mixture weights (probability of each cluster) associated with
                each cluster. Not returned if return_ids_only is True.

        Raises:
            ValueError: A ValueError is raised if the input sequence contains
                unrecognized characters or another serious issue is encountered.
        """
        _, chain_name, arr, _, _, _ = self._full_prep_sequence(seq)
        if chain_name not in ["H", "L"]:
            raise ValueError("The sequence provided does not recognizably "
                    "belong as a heavy or light chain.")

        # We can flatten here, because only one sequence is used as input.
        cluster_probs = self.models["human"][chain_name].predict(arr,
                return_raw_probs=True).flatten()
        best_clusters = np.argsort(cluster_probs)[-nclusters:]

        if return_ids_only:
            return best_clusters, chain_name

        mu_mix = self.models["human"][chain_name].mu_mix[best_clusters,...].copy()
        mixweights = self.models["human"][chain_name].mix_weights[best_clusters].copy()

        return best_clusters, chain_name, mu_mix, mixweights



    def retrieve_cluster(self, cluster_id, chain_type):
        """A convenience function to get the per-position probabilities
        associated with a particular cluster.

        Args:
            cluster_id (int): The id number of the cluster to retrieve. Can
                be generated by calling self.get_closest_clusters or
                self.batch_score_seqs with mode = "assign".
            chain_type (str): One of "H", "L".

        Returns:
            mu_mix (np.ndarray): An array of shape (1, sequence_length,
                21), where 21 is the number of possible AAs. The clusters
                are sorted in order from most to least likely given the
                input sequence.
            mixweights (float): The probability of this cluster in the mixture.
            aas (list): A list of amino acids in standard order. The last
                dimension of mu_mix corresponds to these aas in the order given.
        """
        mu_mix = self.models["human"][chain_type].mu_mix[cluster_id:cluster_id+1,...].copy()
        mixweights = self.models["human"][chain_type].mix_weights[cluster_id].copy()
        return mu_mix, mixweights, self.aa_list



    def convert_sequence_to_array(self, seq):
        """Converts an input sequence to a type uint8_t array where
        the integer at each position indicates the amino acid at that
        position. Can be used in conjunction with the cluster returned by
        retrieve_cluster or get_closest_clusters to determine which amino
        acids are contributing most (or least) to the humanness score.

        Args:
            seq (str): The sequence of interest.

        Returns:
            chain_name (str): The chain type; one of "H", "L".
            arr (np.ndarray): A numpy array of shape (1,M) where M is
                the sequence length after converting to a fixed length
                array.
        """
        return self._simple_prep_sequence(seq)



    def _batch_score(self, seq_array:list, output_scores:list,
            assigned_idx:list, chain_type:str, mode:str = "score",
            mask = None, mask_terminal_dels = False, mask_gaps = False,
            nthreads:int = 2):
        """Scores a batch of sequences -- either heavy or light --
        with the provided mixmodel. Operations are in place so nothing
        is returned.

        Args:
            seq_array (list): A list of np.uint8 numpy arrays of shape (1, n_positions).
            output_scores (list): A list of scores or clusters. Will be modified in-place.
            assigned_idx (list): A list of the indices to which each element of
                seq_array corresponds. Must be the same length as seq_array.
            chain_type (str): One of "H" or "L".
            mode (str): One of 'score', 'assign', 'assign_no_weights', 'classifier'.
                If score, returns the human generative model score.
                If 'assign', provides the most likely cluster number
                for each input sequence. If 'assign_no_weights',
                assigns the closest cluster ignoring mixture weights,
                so that the closest cluster is assigned even if that
                cluster is a low-probability one. If 'classifier',
                assigns a score using the Bayes' rule classifier.
            mask (np.ndarray): Either None or an np.bool array of shape
                (n_positions). If not None, the positions marked False in
                the mask are ignored.
            mask_terminal_dels (bool): If True, terminal deletions are masked.
                This is useful when a sequence contains large unusual terminal
                deletions.
            mask_gaps (bool): If True, all non-filled IMGT positions in the sequence
                are ignored when calculating the score. This is useful when your
                sequence has unusual deletions and you would like to ignore these.
            nthreads (int): The number of threads to use.
        """
        input_array = np.vstack(seq_array)

        n_threads = min(nthreads, len(seq_array))

        if mode == "classifier":
            scores = []
            for species in ("human", "mouse", "rhesus"):
                scores.append(self.models[species][chain_type].score(input_array,
                        mask, mask_terminal_dels, mask_gaps, self.normalize_scores,
                        n_threads = n_threads))
            if chain_type == "H":
                scores.append(self.models["rat"]["H"].score(input_array,
                        mask, mask_terminal_dels, mask_gaps, self.normalize_scores,
                        n_threads = n_threads))
            else:
                scores.append(None)

            scores = self._convert_score_to_classifier(scores[0], scores[1],
                        scores[2], scores[3], chain_type)

        elif mode == "score":
            scores = self.models["human"][chain_type].score(input_array, mask,
                    mask_terminal_dels, mask_gaps, self.normalize_scores,
                    n_threads = n_threads)
            if not self.normalize_scores:
                scores -= self.score_adjustments[chain_type]

        elif mode == "assign":
            scores = self.models["human"][chain_type].predict(input_array, mask = mask,
                    mask_terminal_dels = mask_terminal_dels,
                    mask_gaps = mask_gaps, n_threads = n_threads)
        elif mode == "assign_no_weights":
            scores = self.models["human"][chain_type].predict(input_array, mask = mask,
                    mask_terminal_dels = mask_terminal_dels,
                    mask_gaps = mask_gaps, n_threads = n_threads,
                    use_mixweights = False)

        for i, score in enumerate(scores.tolist()):
            output_scores[assigned_idx[i]] = score



    def _convert_score_to_classifier(self, human_score, mouse_score,
            rhesus_score, rat_score = None, chain_type = "H"):
        """Converts the input generative model scores into a Bayes'
        rule classifier score. Should only be used for testing,
        since a classification approach of this kind is not
        useful for sequences of unknown origin."""
        num_positions = self.models["human"][chain_type].sequence_length

        #This is a default probability assuming equal probability of all
        #amino acids at all positions.
        default_prob = np.exp(np.log(1/21) * num_positions)

        human_probs = np.exp(human_score)
        if rat_score is not None:
            norm_constant = 0.2 * np.exp(mouse_score) + 0.2 * np.exp(rhesus_score) + \
                0.2 * np.exp(rat_score) + 0.2 * human_probs + 0.2 * default_prob
            return 0.2 * human_probs / norm_constant

        norm_constant = 0.25 * np.exp(mouse_score) + 0.25 * np.exp(rhesus_score) + \
                0.25 * human_probs + 0.25 * default_prob
        return 0.25 * human_probs / norm_constant
