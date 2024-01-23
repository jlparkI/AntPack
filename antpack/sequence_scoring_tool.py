"""Provides a tool for scoring individual sequences."""
import os
import copy
import numpy as np
from antpack import SingleChainAnnotator
from .constants.scoring_constants import scoring_constants as constants
from .constants.scoring_constants import allowed_imgt_pos as ahip
from .utilities.model_loader_utils import load_model
from .utilities.data_compile_utils import get_position_dict, get_reverse_position_dict


class SequenceScoringTool():
    """Provides a tool for scoring sequences to determine the likelihood
    they are human, for humanizing them, for assigning them to clusters,
    and for retrieving specific model clusters, together with other
    convenience functions.

    Attributes:
        position_dict (dict): Maps allowed IMGT positions to position numbers.
        rev_position_dict (dict): Maps position numbers back to IMGT positions.
        chain_map (dict): Maps from chain names to a condensed list (both K and L
            map to L).
        models (dict): Contains species- and chain-specific mixture models. This
            will ordinarily contain only human, unless you are considering using
            the model in classifier mode, i.e. by setting offer_classifier_option
            to True (this is more accurate for species present in the training set
            but less accurate if unexpected species are present, so we don't
            recommend it as general practice).
        aa_list (list): The list of allowed / expected AAs. AAs not in here will
            generate an exception.
        aa_dict (dict): Maps AAs to numbers.
        aligner (SingleChainAnnotator): A SingleChainAnnotator for assigning numbering
            to input sequences.
        cdr_mask (dict): Maps the CDR positions from the constants file to positions.
        score_adjustments (dict): Constants to subtract from assigned scores, by chain
            type. These are the median values from training data. Subtracting them just
            places heavy and light chains on the same scale. This behavior can be
            disabled if desired by setting adjusted_scores to False.
        """

    def __init__(self, adjusted_scores = True, offer_classifier_option = False):
        """Class constructor.

        Args:
            adjusted_scores (bool): If True, the median heavy and light chain
                scores from the training data are subtracted from model
                scores. This ensures the heavy and light chain scores are
                more directly comparable, but also means that the scores
                can no longer be directly converted to probabilities.
            offer_classifier_option (bool): If True, the object loads
                additional species-specific models (mouse, rhesus monkey,
                rat) and can use these to function as a classifier. If this
                is False, setting mode='classifier' when calling a scoring
                function will result in an error.
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
        self.aligner = SingleChainAnnotator(chains=["H", "K", "L"], scheme = "imgt",
                compress_init_gaps = False)

        self.cdr_mask = {key:np.ones((len(self.position_dict[key])), bool) for
                key in self.position_dict}

        for chain_type, excluded_pos in ahip.CDR_EXCLUDED.items():
            for pos in excluded_pos:
                self.cdr_mask[chain_type][self.position_dict[chain_type][pos]] = False

        if adjusted_scores:
            self.score_adjustments = {"H":constants.HEAVY_MEDIAN_SCORE,
                    "L":constants.LIGHT_MEDIAN_SCORE}
        else:
            self.score_adjustments = {"H":0.0, "L":0.0}



    def _prep_sequence(self, seq:str):
        """Determines the chain type of a sequence and aligns it
        using the IMGT numbering.

        Args:
            seq (str): An input sequence

        Returns:
            seq_extract (str): The extracted and aligned sequence.
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
        numbering, _, chain_name, _ = self.aligner.analyze_seq(seq)
        chain_name = self.chain_map[chain_name]
        seq_extract, bad_gaps, bad_positions, backmap = self._extract_sequence(seq,
                    numbering, chain_name)
        seq_arr = np.array([[self.aa_dict[letter] for letter
                in seq_extract]], dtype=np.uint8)
        return seq_extract, chain_name, seq_arr, bad_gaps, bad_positions, backmap



    def score_seq(self, seq, mask_term_dels = False,
            mask_all_gaps = False, return_diagnostics = False):
        """Scores a single sequence, returning either just the
        score or some additional information if requested. This
        function is considerably slower than batch_score_seqs
        for large numbers of sequences and does not multithread
        but offers more fine-tuned control over the details of
        scoring and can return additional information about a
        sequence.

        Args:
            seq (str): The sequence. May be either heavy chain or
                light.
            mask_term_dels (bool): If True, ignore N or C
                terminal deletions when calculating probabilities.
            mask_gaps (bool): If True, ignore all gaps when calculating
                probabilities. This can occasionally be useful if you
                want to see what the score would look like ignoring
                a large deletion. If this is True it overrides
                mask_term_dels.
            return_diagnostics (bool): If True, return a list of
                unexpected gaps and a list of unexpected insertions
                (if any) found when aligning the sequence and the
                assigned chain type ("H" or "L").

        Returns:
            score (float): log( p(x) ).
            bad_gaps (list): A list of unexpected gaps at sites where
                a deletion is unusual in the IMGT numbering system.
                Only returned if return_diagnostics is True.
            bad_positions (list): A list of unexpected insertions at
                sites where an insertion is not expected in IMGT.
                Only returned if return_diagnostics is True.
        """
        _, chain_name, arr, bad_gaps, bad_positions, _ = self._prep_sequence(seq)
        if chain_name not in ["H", "L"]:
            if return_diagnostics:
                return np.nan, bad_gaps, bad_positions, chain_name
            return np.nan

        if mask_all_gaps:
            human_score = float(self.models["human"][chain_name].gapped_score(
                    arr, n_threads = 1)[0])

        elif mask_term_dels:
            start_col, end_col = self._get_term_cols(bad_gaps,
                    self.models["human"][chain_name].sequence_length)
            human_score = float(self.models["human"][chain_name].terminal_masked_score(
                    arr, n_threads = 1, start_col = start_col, end_col = end_col)[0])

        else:
            human_score = float(self.models["human"][chain_name].score(arr, n_threads = 1)[0])

        human_score -= self.score_adjustments[chain_name]

        if return_diagnostics:
            return human_score, bad_gaps, bad_positions, chain_name
        return human_score



    def batch_score_seqs(self, seq_list, mode = "score"):
        """Scores a list of sequences in batches. Substantially faster than
        single seq scoring but does not offer the option to retrieve
        diagnostic info. Can also be used to assign a large number of
        sequences to clusters as well.

        Args:
            seq_list (str): The list of input sequences. May contain both
                heavy and light.
            mode (str): One of 'score', 'assign', 'classifier'.
                If score, returns the human generative model score.
                If 'assign', provides the most likely cluster number
                for each input sequence. If 'classifier',
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

        for i, seq in enumerate(seq_list):
            _, chain_name, arr, _, _, _ = self._prep_sequence(seq)

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
                    heavy_idx, "H", mode)

        if len(light_arr) > 0:
            self._batch_score(light_arr, output_scores,
                    light_idx, "L", mode)

        return output_scores



    def suggest_mutations(self, seq:str, s_thresh:float = 1.25):
        """Takes an input sequence, scores it per position,
        uses the nclusters closest clusters to determine which
        modification would be most likely to have an impact,
        suggest mutations and report both the mutations and the
        new score. CDRs are excluded.

        Args:
            seq (str): The sequence to update.
            s_thresh (float): The maximum percentage by which
                the score can shift before backmutation stops.
                Smaller values (closer to 1) will prioritize
                increasing the score over preserving the original
                sequence. Larger values will prioritize preserving
                the original sequence.

        Returns:
            initial_score (float): The score of the sequence pre-modification.
            final_scores (float): The scores of the sequence after each mutation
                is adopted (in sequential order).
            mutations (list): The suggested mutations in AA_position_newAA format,
                where position is the IMGT number for the mutation position. These
                are in sequential order (the same as final_scores).
            updated_seq (list): The updated sequences after each mutation with all
                gaps removed. (This may be a different length from the
                input sequence if the suggested mutation is a deletion or
                an insertion).
        """
        original_seq, chain_name, original_arr, _, _, backmap = self._prep_sequence(seq)
        if chain_name not in ["H", "L"]:
            raise ValueError("The sequence provided does not recognizably "
                    "belong as a heavy or light chain.")

        mixmodel = self.models["human"][chain_name]

        updated_arr = original_arr.copy()
        _, _, best_cluster, _ = self.get_closest_clusters(seq, 1)
        best_cluster = np.log(best_cluster[0,...].clip(min=1e-16))
        best_aas = best_cluster.argmax(axis=1)
        mask = self.cdr_mask[chain_name].copy()
        updated_arr[0,mask] = best_aas[mask]
        starting_score = mixmodel.score(updated_arr, n_threads = 1)[0]
        updated_score = copy.copy(starting_score)

        while updated_score > s_thresh * starting_score:
            score_shifts = best_cluster[np.arange(best_cluster.shape[0]), updated_arr.flatten()] - \
                    best_cluster[np.arange(best_cluster.shape[0]), original_arr.flatten()]
            score_shifts[~mask] = np.inf
            score_shifts[updated_arr[0,:]==original_arr[0,:]] = np.inf
            weakest_position = score_shifts.argmin()
            proposal_arr = updated_arr.copy()
            if proposal_arr[0,weakest_position] == original_arr[0,weakest_position]:
                break
            proposal_arr[0,weakest_position] = original_arr[0,weakest_position]
            updated_score = mixmodel.score(proposal_arr, n_threads = 1)[0]
            if updated_score > s_thresh * starting_score:
                updated_arr = proposal_arr
            else:
                break

        updated_seq = [self.aa_list[aa] for aa in updated_arr[0,...].tolist()]

        all_mutations = []

        for i, (original_aa, new_aa) in enumerate(zip(original_seq, updated_seq)):
            if new_aa != original_aa:
                if i in backmap:
                    mutation_description = (f"{original_aa}_"
                        f"{backmap[i]+1}_{new_aa}")
                else:
                    counter = i
                    while counter not in backmap and counter > 0:
                        counter -= 1
                    if counter == 0:
                        mutation_description = (f"{original_aa}_"
                            f"{0}_{new_aa}")
                    else:
                        mutation_description = (f"{original_aa}_"
                            f"{backmap[counter]}insert_{new_aa}")
                all_mutations.append(mutation_description)

        score = mixmodel.score(updated_arr, n_threads = 1) - \
                        self.score_adjustments[chain_name]
        return score[0], all_mutations, "".join(updated_seq).replace("-", "")



    def get_closest_clusters(self, seq:str, nclusters:int = 1, return_ids_only = False):
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
        _, chain_name, arr, _, _, _ = self._prep_sequence(seq)
        if chain_name not in ["H", "L"]:
            raise ValueError("The sequence provided does not recognizably "
                    "belong as a heavy or light chain.")

        # We can flatten here, because only one sequence is used as input.
        cluster_probs = self.models["human"][chain_name].predict_proba(arr).flatten()
        best_clusters = np.argsort(cluster_probs)[-nclusters:]

        if return_ids_only:
            return best_clusters, chain_name

        mu_mix = self.models["human"][chain_name].mu_mix[best_clusters,...].copy()
        mixweights = self.models["human"][chain_name].mix_weights[best_clusters].copy()

        return best_clusters, chain_name, mu_mix, mixweights



    def retrieve_cluster(self, cluster_id, chain_type = "H"):
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
        mixweights = self.models["human"][chain_type].mix_weights[cluster_id].copy()[0]
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
        _, chain_name, arr, _, _, _ = self._prep_sequence(seq)
        return chain_name, arr



    def _batch_score(self, seq_array:list, output_scores:list,
            assigned_idx:list, chain_type:str, mode = "score",
            nthreads = 2):
        """Scores a batch of sequences -- either heavy or light --
        with the provided mixmodel. Operations are in place so nothing
        is returned.

        Args:
            seq_array (list): A list of np.uint8 numpy arrays of shape (1, n_positions).
            output_scores (list): A list of scores or clusters. Will be modified in-place.
            assigned_idx (list): A list of the indices to which each element of
                seq_array corresponds. Must be the same length as seq_array.
            chain_type (str): One of "H" or "L".
            mode (str): One of 'score', 'classifier. Score calculates
                the log likelihood for each datapoint; classifier builds a
                simple Bayes' rule classifier using generative model scores.
            nthreads (int): The number of threads to use.
        """
        input_array = np.vstack(seq_array)
        n_threads = min(nthreads, len(seq_array))

        if mode == "classifier":
            scores = self.models["human"][chain_type].score(input_array, n_threads)
            mouse_scores = self.models["mouse"][chain_type].score(input_array, n_threads)
            rhesus_scores = self.models["rhesus"][chain_type].score(input_array, n_threads)
            if chain_type == "H":
                rat_scores = self.models["rat"]["H"].score(input_array, n_threads)
            else:
                rat_scores = None

            scores = self._convert_score_to_classifier(scores, mouse_scores,
                        rhesus_scores, rat_scores, chain_type)

        elif mode == "score":
            scores = self.models["human"][chain_type].score(input_array, n_threads)
            scores -= self.score_adjustments[chain_type]

        elif mode == "assign":
            scores = self.models["human"][chain_type].predict(input_array, n_threads)

        for i, score in enumerate(scores.tolist()):
            output_scores[assigned_idx[i]] = score




    def _extract_sequence(self, input_seq, input_numbering, chain_type):
        """Converts a numbered sequence into a gapped sequence that can
        be encoded.

        Args:
            input_seq (str): An input sequence that has already been
                run through the aligner.
            input_numbering (list): The IMGT numbering assigned by the
                aligner.
            chain_type (str): One of 'H', 'K', 'L'. If not one of these,
                an all-gap sequence is returned.

        Returns:
            output_seq (list): A gapped output sequence with insertions
                at beginning and end removed and with unexpected insertions
                left out.
            unexpected_positions (list): A list of insertions which are
                unusual in the IMGT scheme.
            unexpected_gaps (list): A list of positions at which a gap was
                encountered which are not normally a gap in the IMGT scheme.
            backmap (dict): A dictionary mapping back from the numbering
                positions to the original position numbers.
        """
        if chain_type == "unknown":
            output_seq = ["-" for i in range(len(input_seq))]
            return output_seq, [], [], {}

        position_dict = self.position_dict[chain_type]
        rev_dict = self.rev_position_dict[chain_type]
        output_seq = ["-" for i in range(len(position_dict))]
        unexpected_positions = []
        backmap = {}
        input_counter = 0

        for position, aa in zip(input_numbering, input_seq):
            if position == "-":
                input_counter += 1
                continue
            if position not in position_dict:
                unexpected_positions.append(position)
                input_counter += 1
                continue
            output_seq[position_dict[position]] = aa
            if aa != "-":
                backmap[position_dict[position]] = input_counter
                input_counter += 1

        unexpected_gaps = [rev_dict[position] for position, aa in
                enumerate(output_seq) if aa == "-"]
        return output_seq, unexpected_gaps, unexpected_positions, backmap


    def _get_term_cols(self, gaps, seqlen):
        """If there are n-terminal and c-terminal deletions,
        this function finds them and assigns start and end
        columns for masking."""
        start_col = 0
        for i in range(1, len(gaps)):
            if str(i) != gaps[i-1]:
                break
            start_col = i

        end_col, counter = seqlen, len(gaps) - 1
        for i in range(seqlen - 1, -1, -1):
            if str(i) != gaps[counter]:
                break
            end_col = i
            counter -= 1
        return start_col, end_col


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
