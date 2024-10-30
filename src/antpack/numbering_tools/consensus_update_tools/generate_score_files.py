"""Contains the tools needed to generate .npy arrays for each
chain type. The .npy arrays are used to score new sequences
so that they are aligned and numbered."""
import os
import numpy as np
from Bio.Align import substitution_matrices
from ..constants import aho_default_params as aho_dp
from ..constants import imgt_default_params as imgt_dp
from ..constants import kabat_default_params as kabat_dp
from ..constants import martin_default_params as martin_dp
from ..constants import all_scheme_default_params as shared_dp


def build_consensus_files():
    fdir = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(fdir, "..", "consensus_data")
    cdir = os.getcwd()

    for scheme_name in ["kabat", "martin", "imgt", "ctermfinder", "aho"]:
        for chain_name in ["H", "K", "L"]:
            confile = f"{scheme_name.upper()}_CONSENSUS_{chain_name}.txt"
            build_scoring_files(data_path, cdir, confile,
                        chain_type = chain_name,
                        scheme = scheme_name)


def build_scoring_files(target_dir, current_dir, consensus_file,
        chain_type = "H", scheme = "kabat"):
    """Builds a scoring matrix for the amino acids at each position for each
    chain type using predefined consensus files.

    Args:
        target_dir (str): The filepath of the output directory.
        current_dir (str): The filepath of the current directory.
    """
    if chain_type not in ["H", "K", "L"]:
        raise ValueError("Currently only H, K, L chains are supported.")
    os.chdir(target_dir)
    consensus_list = load_consensus_file(consensus_file)
    if scheme == "kabat":
        save_consensus_array(consensus_list, chain_type, kabat_dp, "kabat")
    elif scheme == "martin":
        save_consensus_array(consensus_list, chain_type, martin_dp, "martin")
    elif scheme == "imgt":
        save_consensus_array(consensus_list, chain_type, imgt_dp, "imgt")
    elif scheme == "aho":
        save_consensus_array(consensus_list, chain_type, aho_dp, "aho")
    elif scheme == "ctermfinder":
        save_cterm_finder_array(consensus_list, chain_type)

    os.chdir(current_dir)



def load_consensus_file(consensus_file):
    """Loads a specified consensus file and stores it as a list
    of lists which can be converted to a scoring matrix."""
    consensus_list = []
    with open(consensus_file, "r", encoding="utf-8") as fhandle:
        position_number = 0
        for line in fhandle:
            if line.startswith("#") or line.startswith("/"):
                continue
            if len(line) <= 1:
                continue
            new_position_number = int(line.split(",")[0])
            if new_position_number != position_number + 1:
                raise ValueError(f"Consensus file {consensus_file} has incorrect formatting!")
            position_number = int(line.split(",")[0])
            observed_aas = line.strip().split(",")[1:]
            if len(observed_aas) < 1:
                raise ValueError(f"Consensus file {consensus_file} has incorrect formatting!")
            consensus_list.append(observed_aas)
    return consensus_list



def save_consensus_array(consensus_list, chain_type, constants = imgt_dp, scheme = "imgt"):
    """Converts a list of sequences for a specific chain type
    to an array with the score for each possible amino acid substitution at
    each position, including gap penalties. For IMGT (as for other numbering
    schemes), we prefer to place insertions at specific places, so we tailor
    the gap penalties to encourage this. Meanwhile, other positions are
    HIGHLY conserved, so we tailor the penalties to encourage this as well.
    IMGT numbers from 1 so we have to adjust for this."""
    allowed_amino_acids = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N",
            "P", "Q", "R", "S", "T", "V", "W", "Y"}
    #The list of allowed amino acids. Must be in alphabetical order, very important.
    allowed_aa_list = sorted(list(allowed_amino_acids))

    blosum = substitution_matrices.load("BLOSUM62")
    blosum_key = {letter:i for i, letter in enumerate(blosum.alphabet)}


    if chain_type.endswith("K") or chain_type.endswith("L"):
        conserved_positions = constants.light_conserved_positions
        if scheme == "aho" and chain_type.endswith("K"):
            special_positions = constants.kappa_special_positions
        else:
            special_positions = constants.light_special_positions
        cdrs = constants.light_cdrs
        npositions = constants.NUM_LIGHT

    elif chain_type.endswith("H"):
        conserved_positions = constants.heavy_conserved_positions
        special_positions = constants.heavy_special_positions
        cdrs = constants.heavy_cdrs
        npositions = constants.NUM_HEAVY

    else:
        return


    key_array = np.zeros((npositions, 23))

    for i, observed_aas in enumerate(consensus_list[:npositions]):
        position = i + 1

        #Choose the gap penalty
        #for template (column 20) and for query (column 21) of key array.
        if position in conserved_positions:
            key_array[i,21:] = shared_dp.HIGHLY_CONSERVED_GAP_PENALTY
        elif position in special_positions:
            key_array[i,21] = special_positions[position][0]
            key_array[i,22:] = special_positions[position][1]
        elif position in cdrs:
            key_array[i,21:] = cdrs[position]
        elif position in shared_dp.n_terminal_gap_positions:
            key_array[i,21] = shared_dp.n_terminal_gap_positions[position][0]
            key_array[i,22:] = shared_dp.n_terminal_gap_positions[position][1]
        else:
            key_array[i,21] = shared_dp.DEFAULT_QUERY_GAP_PENALTY
            key_array[i,22] = shared_dp.DEFAULT_TEMPLATE_GAP_PENALTY

        #Next, fill in the scores for other amino acid substitutions. If a conserved
        #residue, use the ones we specify here. Otherwise, use the best possible
        #score given the amino acids observed in the alignments. If the only
        #thing observed in the alignments is gaps, no penalty is applied.
        for j, letter in enumerate(allowed_aa_list):
            letter_blosum_idx = blosum_key[letter]
            if position in conserved_positions:
                if letter == conserved_positions[position]:
                    key_array[i,j] = shared_dp.HIGHLY_CONSERVED_BONUS
            else:
                key_array[i,j] = max([blosum[letter_blosum_idx, blosum_key[k]] if
                    k != "-" else 0 for k in observed_aas])

        #Finally, fill in the score for X. X should be assumed to be the best
        #possible amino acid that could be present at that position.
        key_array[i,20] = key_array[i,:20].max()
        if position in conserved_positions:
            key_array[i,20] = shared_dp.HIGHLY_CONSERVED_BONUS

    np.save(f"{scheme.upper()}_CONSENSUS_{chain_type}.npy", key_array)



def save_cterm_finder_array(consensus_list, chain_type):
    """Converts a list of sequences for a specific chain type
    to an array with the score for each possible amino acid substitution at
    each position SPECIFIC to the specialized c-terminal finder. Each
    covers an array of shape """
    allowed_amino_acids = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N",
            "P", "Q", "R", "S", "T", "V", "W", "Y"}
    #The list of allowed amino acids. Must be in alphabetical order, very important.
    allowed_aa_list = sorted(list(allowed_amino_acids))

    blosum = substitution_matrices.load("BLOSUM62")
    blosum_key = {letter:i for i, letter in enumerate(blosum.alphabet)}

    key_array = np.zeros((len(consensus_list), 21))

    for i, observed_aas in enumerate(consensus_list):

        # Fill in the scores for amino acid substitutions. For this (unlike
        # the aligner), we do not provide a large bonus for matching at
        # conserved positions, and we do not need to consider gaps.If the only
        # thing observed in the alignments is gaps, no penalty is applied.
        for j, letter in enumerate(allowed_aa_list):
            letter_blosum_idx = blosum_key[letter]
            key_array[i,j] = max([blosum[letter_blosum_idx, blosum_key[k]] if
                    k != "-" else 0 for k in observed_aas])

        #Finally, fill in the score for X. X should be assumed to be the best
        #possible amino acid that could be present at that position.
        key_array[i,20] = key_array[i,:20].max()

    np.save(f"CTERMFINDER_CONSENSUS_{chain_type}.npy", key_array)
