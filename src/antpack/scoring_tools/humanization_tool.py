"""Provides a tool for humanizing input sequences according to
a preset formula."""
import copy
import numpy as np
from .sequence_scoring_tool import SequenceScoringTool
from .scoring_constants import allowed_imgt_pos as ahip


class HumanizationTool():
    """Tool for humanizing sequences."""

    def __init__(self):
        """Class constructor."""

        self.score_tool = SequenceScoringTool(offer_classifier_option = False)

        position_dict = self.score_tool.position_dict

        self.cdr_mask = {k:np.ones((len(v)), bool) for
                (k,v) in position_dict.items()}

        for chain_type, excluded_pos in ahip.CDR_EXCLUDED.items():
            for pos in excluded_pos:
                self.cdr_mask[chain_type][position_dict[chain_type][pos]] = False



    def suggest_mutations(self, seq:str, excluded_positions:list = [],
            s_thresh:float = 1.25):
        """Takes an input sequence, scores it per position,
        uses the nclusters closest clusters to determine which
        modification would be most likely to have an impact,
        suggest mutations and report both the mutations and the
        new score. CDRs are excluded, together with user-specified
        excluded positions.

        Args:
            seq (str): The sequence to update.
            s_thresh (float): The maximum percentage by which
                the score can shift before backmutation stops.
                Smaller values (closer to 1) will prioritize
                increasing the score over preserving the original
                sequence. Larger values will prioritize preserving
                the original sequence.
            excluded_positions (list): A list of strings (IMGT position numbers)
                indicating positions which should not be changed. This enables
                the user to mask key residues, Vernier zones etc if so
                desired.

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
        original_seq, chain_name, original_arr, _, _, backmap = \
                self.score_tool._full_prep_sequence(seq)
        if chain_name not in ["H", "L"]:
            raise ValueError("The sequence provided does not recognizably "
                    "belong as a heavy or light chain.")

        mixmodel = self.score_tool.models["human"][chain_name]

        updated_arr = original_arr.copy()
        _, _, best_cluster, _ = self.score_tool.get_closest_clusters(seq, 1, False)

        best_cluster = np.log(best_cluster[0,...].clip(min=1e-16))
        best_aas = best_cluster.argmax(axis=1)
        mask = self.cdr_mask[chain_name].copy()
        for position in excluded_positions:
            mask[self.score_tool.position_dict[chain_name][position]] = False

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

        updated_seq = [self.score_tool.aa_list[aa] for aa
                in updated_arr[0,...].tolist()]

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
                        self.score_tool.score_adjustments[chain_name]
        return score[0], all_mutations, "".join(updated_seq).replace("-", "")
