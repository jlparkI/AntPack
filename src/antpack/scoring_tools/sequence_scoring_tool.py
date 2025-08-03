"""Provides a tool for scoring individual sequences."""
import os
import numpy as np
from ..numbering_tools import SingleChainAnnotator
from .scoring_constants import scoring_constants as constants
from .scoring_constants import allowed_imgt_pos as ahip
from ..utilities.model_loader_utils import load_model
from ..scoring_tools.scoring_constants import allowed_imgt_pos as ahip
from antpack.antpack_cpp_ext import SequenceTemplateAligner



class SequenceScoringTool():
    """Tool for scoring sequences."""

    def __init__(self, offer_classifier_option = False,
            normalization = "none", max_threads = 2):
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
            max_threads (int): The maximum number of threads to use when
                making predictions.
        """
        project_dir = os.path.abspath(os.path.dirname(__file__))

        self.template_aligners = {"H":SequenceTemplateAligner(
            ahip.heavy_allowed_positions, "H", "imgt"),
                "L":SequenceTemplateAligner(
                    ahip.light_allowed_positions, "L", "imgt")
                }
        self.chain_map = {"H":"H", "K":"L", "L":"L", "":"unknown"}

        if offer_classifier_option:
            self.models = {
                    "human":{"H":load_model(project_dir, "heavy",
                                                max_threads=max_threads),
                        "L":load_model(project_dir, "light",
                                                max_threads=max_threads)},
                    "mouse":{"H":load_model(project_dir, "heavy", species="mouse",
                                                max_threads=max_threads),
                        "L":load_model(project_dir, "light", species="mouse",
                                                max_threads=max_threads)},
                    "rhesus":{"H":load_model(project_dir, "heavy", species="rhesus",
                                                max_threads=max_threads),
                        "L":load_model(project_dir, "light", species="rhesus",
                                                max_threads=max_threads)},
                    "rat":{"H":load_model(project_dir, "heavy", species="rat",
                                                max_threads=max_threads) }
                }
        else:
            self.models = {"human":{"H":load_model(project_dir, "heavy",
                                    max_threads=max_threads),
                        "L":load_model(project_dir, "light",
                                    max_threads=max_threads)} }

        self.aa_list = constants.aa_list
        self.aligner = SingleChainAnnotator(chains=["H", "K", "L"],
                scheme = "imgt")

        self.normalize_scores = False
        self.score_adjustments = {"H":0.0, "L":0.0}
        if normalization == "training_set_adjust":
            self.score_adjustments = {"H":constants.HEAVY_MEDIAN_SCORE,
                    "L":constants.LIGHT_MEDIAN_SCORE}
        elif normalization == "normalize":
            self.normalize_scores = True



    def _prep_sequence(self, seq:str):
        """Determines the chain type of a sequence and aligns it
        using the IMGT numbering.

        Args:
            seq (str): An input sequence

        Returns:
            chain_name (str): One of "L", "H". Note that "K" and "L"
                are combined since the same mixture model is used for
                "K" and "L".
            aligned_seq (str): The sequence aligned to the MSA used to
                train the mixture model, with any positions not present
                in the mixture model disregarded and gaps added where
                appropriate.
        """
        numbering, _, chain_type, _ = self.aligner.analyze_seq(seq)
        chain_type = self.chain_map[chain_type]
        if chain_type == "unknown":
            return chain_type, len(seq) * '-'

        output_seq = self.template_aligners[chain_type].align_sequence(
                seq, numbering, False)
        return chain_type, output_seq




    def get_standard_positions(self, chain_type:str):
        """Returns a list of the standard positions used by SAM when scoring a
        sequence. IMGT numbered positions outside this set are ignored when
        assigning a score."""
        if chain_type == "H":
            return self.template_aligners["H"].get_template_numbering()
        if chain_type in ("K", "L"):
            return self.template_aligners["L"].get_template_numbering()
        raise ValueError("Unrecognized chain type supplied.")



    def get_standard_mask(self, chain_type:str, region:str = "fmwk1",
            cdr_labeling_scheme="imgt"):
        """Returns a mask for all regions EXCEPT the one you specify.
        You can then use this mask as input to score_seqs to see what
        the score would be if only that region were included.

        Args:
            chain_type (str): One of "H", "L".
            region (str): One of 'fmwk1', 'fmwk2', 'fmwk3',
                'fmwk4', 'cdr1', 'cdr2', 'cdr3', 'fmwk', 'cdr'.
                This function will construct a mask that excludes
                all other regions.
            cdr_labeling_scheme (str): The numbering scheme used for
                humanness calculations is IMGT, but for generating a
                mask, you can use a different scheme to assign CDRs
                if desired. This value can be one of 'aho', 'imgt',
                'kabat', 'martin'.

        Returns:
            mask (list): A numpy array of the same length as the
                list returned by "get_standard_positions()". Only
                positions marked True will be used when scoring a
                sequence. You can further modify this mask if needed
                and can pass it to "score_seqs" to use when scoring
                sequences.

        Raises:
            ValueError: A ValueError is raised if unexpected inputs are supplied.
        """
        if chain_type not in ["H", "L"]:
            raise ValueError("chain_type must be one of 'H', 'L'.")
        if region not in ("fmwk1", "fmwk2", "fmwk3", "fmwk4",
                "cdr1", "cdr2", "cdr3", "fmwk", "cdr"):
            raise ValueError("Unexpected region supplied.")

        numbering = self.template_aligners[chain_type].get_template_numbering()
        region_labels = self.aligner.assign_cdr_labels(numbering, chain_type,
                cdr_labeling_scheme)
        mask = [r.startswith(region) for r in region_labels]
        return np.array(mask, dtype=np.bool)




    def get_diagnostic_info(self, seq:str):
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
                chain type to which the sequence was aligned. "Unknown" if there
                was an error in numbering.
        """
        numbering, _, chain_type, _ = self.aligner.analyze_seq(seq)
        chain_type = self.chain_map[chain_type]
        if chain_type == "unknown":
            return [], [], chain_type

        template_numbering = self.template_aligners[chain_type].get_template_numbering()

        aligned_seq = self.template_aligners[chain_type].align_sequence(
                seq, numbering, False)
        bad_gaps = [template_numbering[i] for i,l in enumerate(aligned_seq) if l=='-']
        forward_idx = self.template_aligners[chain_type].retrieve_alignment_forward_numbering(
                seq, numbering, False)

        bad_positions = [numbering[i] for i,l in enumerate(forward_idx) if l==-1]
        return bad_gaps, bad_positions, chain_type



    def score_seqs(self, seq_list, custom_light_mask:list = None,
            custom_heavy_mask:list = None, mask_terminal_dels:bool = False,
            mask_gaps:bool = False, mode:str = "score"):
        """Scores a list of sequences in batches or assigns them
        to clusters. Can be used in conjunction with a user-supplied
        mask (for positions to ignore) and in conjunction with Substantially faster than
        single seq scoring but does not offer the option to retrieve
        diagnostic info. Can also be used to assign a large number of
        sequences to clusters as well.

        Args:
            seq_list (str): The list of input sequences. May contain both
                heavy and light.
            custom_light_mask (list): Either None or an array generated by
                self.get_standard_mask indicating positions to ignore. This can be useful
                if you just want to score a specific region, or if there is a large
                deletion that should be ignored.
            custom_heavy_mask (list): Either None or an array generated by
                self.get_standard_mask indicating positions to ignore. This can be useful
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
        chain_data = {"H":{"seqs":[], "idx":[], "mask":custom_heavy_mask},
                "L":{"seqs":[], "idx":[], "mask":custom_light_mask}}
        output_scores = np.zeros((len(seq_list)))

        for i, seq in enumerate(seq_list):
            chain_name, aligned_seq = self._prep_sequence(seq)
            if chain_name not in ("L", "H"):
                output_scores[i] = np.nan
            else:
                chain_data[chain_name]["seqs"].append(aligned_seq)
                chain_data[chain_name]["idx"].append(i)

        for chain_name, relevant_data in chain_data.items():
            if len(relevant_data["seqs"]) == 0:
                continue

            assigned_idx = relevant_data["idx"]
            mask = relevant_data["mask"]
            sequences = relevant_data["seqs"]

            if mode == "classifier":
                scores = []
                for species in ("human", "mouse", "rhesus"):
                    scores.append(self.models[species][chain_name].score(sequences,
                            mask=mask, mask_terminal_dels=mask_terminal_dels,
                            mask_gaps=mask_gaps, normalize_scores=self.normalize_scores))
                if chain_name == "H":
                    scores.append(self.models["rat"]["H"].score(sequences,
                            mask=mask, mask_terminal_dels=mask_terminal_dels,
                            mask_gaps=mask_gaps, normalize_scores=self.normalize_scores))
                else:
                    scores.append(None)

                scores = self._convert_score_to_classifier(scores[0], scores[1],
                            scores[2], scores[3], chain_name)

            elif mode == "score":
                scores = self.models["human"][chain_name].score(sequences,
                        mask=mask, mask_terminal_dels=mask_terminal_dels,
                        mask_gaps=mask_gaps,
                        normalize_scores=self.normalize_scores)
                if not self.normalize_scores:
                    scores -= self.score_adjustments[chain_name]

            elif mode == "assign":
                scores = self.models["human"][chain_name].predict(sequences,
                        mask=mask, mask_terminal_dels=mask_terminal_dels,
                        mask_gaps=mask_gaps)
            elif mode == "assign_no_weights":
                scores = self.models["human"][chain_name].predict(sequences,
                        mask=mask, mask_terminal_dels=mask_terminal_dels,
                        mask_gaps=mask_gaps, use_mixweights=False)

            for i, score in enumerate(scores.tolist()):
                output_scores[assigned_idx[i]] = score

        return output_scores



    def get_closest_clusters(self, seq:str, nclusters:int = 1):
        """Gets the closest cluster(s) for a given sequence. These can be used
        for humanization, to determine which amino acids in the input sequence
        are most problematic, or to generate new sequences containing motifs of interest.
        To get the cluster parameters associated with clusters identified by this
        function, call "retrieve_cluster".

        Args:
            seq (str): The input sequence.
            nclusters (int): The number of clusters to retrieve.

        Returns:
            cluster_idx (np.ndarray): The index number for each cluster in the
                original model that is returned. Sorted from highest probability
                cluster to lowest probability.
            chain_name (str): One of "H" or "L", indicating whether chain is
                heavy or light (K and L are both mapped to L).
        """
        chain_name, aligned_seq = self._prep_sequence(seq)
        if chain_name not in ["H", "L"]:
            raise ValueError("The sequence provided does not recognizably "
                    "belong as a heavy or light chain.")

        # We can flatten here, because only one sequence is used as input.
        cluster_probs = self.models["human"][chain_name].predict_proba(
                [aligned_seq]).flatten()
        best_clusters = np.argsort(cluster_probs)[-nclusters:]
        return best_clusters, chain_name



    def retrieve_cluster(self, cluster_id:int, chain_type:str):
        """A convenience function to get the per-position probabilities
        associated with a particular cluster.

        Args:
            cluster_id (int): The id number of the cluster to retrieve. Can
                be generated by calling self.get_closest_clusters or
                self.batch_score_seqs with mode = "assign".
            chain_type (str): One of "H", "L".

        Returns:
            mu_mix (np.ndarray): An array of shape (1, sequence_length,
                21), where 21 is the number of possible AAs.
            mixweights (float): The probability of this cluster in the mixture.
            aas (list): A list of amino acids in standard order. The last
                dimension of mu_mix corresponds to these aas in the order given.
        """
        params = self.models["human"][chain_type].get_model_parameters()
        mu_mix = params[0][cluster_id:cluster_id+1,...].copy()
        mixweights = params[1][cluster_id].copy()
        return mu_mix, mixweights, self.aa_list



    def convert_sequence_to_array(self, seq:str):
        """Converts an input sequence to a type uint8_t array where
        each integer indicates the amino acid at that position.

        Args:
            seq (str): The sequence of interest.

        Returns:
            chain_name (str): The chain type; one of "H", "L"
                or "unknown" if there is an error.
            arr (np.ndarray): A numpy array of shape (1,M) where M is
                the length after converting to a fixed length array.
                nan is returned if there is an error numbering the
                sequence.
        """
        chain_type, aligned_seq = self._prep_sequence(seq)
        if chain_type == "unknown":
            return chain_type, np.full(len(aligned_seq), np.nan)
        output_array = np.zeros((1, len(aligned_seq)), dtype=np.uint8)
        self.models["human"][chain_type].encode_input_seqs([aligned_seq],
                output_array)
        return chain_type, output_array



    def calc_per_aa_probs(self, seq:str, cluster_id:int):
        """Calculate the log probability of each amino acid in
        the input sequence given a specified cluster number. To
        get the cluster number of the cluster closest to your
        input sequence, call get_closest_clusters.

        Args:
            seq (str): The sequence of interest.
            cluster_id (int): The index of the cluster you
                would like to use to generate these probabilities.
                Call get_closest_clusters to find those closest
                to your input sequence.

        Returns:
            chain_type (str): One of "H", "L", or "unknown" if there
                is an error numbering your sequence.
            logprobs (np.ndarray): An array of shape (M) where M
                is the length of your input sequence. nan is returned
                if there is an error numbering the sequence.
        """
        chain_type, seq_arr = self.convert_sequence_to_array(seq)
        if chain_type == "unknown":
            return chain_type, np.full(seq_arr.shape[1], np.nan)
        mu_mix, _, _ = self.retrieve_cluster(cluster_id,
                chain_type)
        # First, remove all gapped positions...
        idx = np.where(seq_arr<=20)[0]
        mu_mix = mu_mix[0,:,idx]
        seq_arr = seq_arr[idx]

        # Next, extract the probabilities at filled positions and
        # take the log.
        mu_mix = mu_mix[seq_arr, np.arange(seq_arr.shape[0])]
        mu_mix = np.log(mu_mix.clip(min=1e-16))
        return chain_type, mu_mix


    def _convert_score_to_classifier(self, human_score, mouse_score,
            rhesus_score, rat_score = None, chain_type = "H"):
        """Converts the input generative model scores into a Bayes'
        rule classifier score. Should only be used for testing,
        since a classification approach of this kind is not
        useful for sequences of unknown origin."""
        num_positions = self.models["human"][chain_type].get_specs()[1]

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
