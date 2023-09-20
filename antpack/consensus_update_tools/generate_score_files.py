"""Contains the tools needed to convert a stockholm alignment of
all the chain types into a set of .npy arrays and a CONSENSUS.txt
file with consensus sequences for each chain type. The .npy arrays
are used to score new sequences so that they are correctly
aligned and numbered."""
import os
import numpy as np
from Bio.Align import substitution_matrices
from ..constants.allowed_inputs import allowed_aa_list
from ..constants import imgt_default_params as imgt_dp
from ..constants import kabat_default_params as kabat_dp
from ..constants import martin_default_params as martin_dp


def build_alternative_scoring(target_dir, current_dir, consensus_file,
        chain_type = "H", scheme = "kabat"):
    """Builds a scoring matrix for the amino acids at each position for each
    chain type for non-IMGT schemes. The consensus file for non-IMGT schemes
    is predefined, so we merely need to load it then convert it to a scoring
    matrix.

    Args:
        target_dir (str): The filepath of the output directory.
        current_dir (str): The filepath of the current directory.
    """
    if chain_type not in ["H", "K", "L"]:
        raise ValueError("For schemes other than IMGT, only H, K, L chains are supported.")
    os.chdir(target_dir)
    consensus_list = load_consensus_file(consensus_file)
    if scheme == "kabat":
        save_consensus_array(consensus_list, chain_type, kabat_dp, "kabat")
    elif scheme == "martin":
        save_consensus_array(consensus_list, chain_type, martin_dp, "martin")
    os.chdir(current_dir)



def build_imgt_consensus_files(target_dir, current_dir, alignment_fname):
    """Builds a consensus for the amino acids at each position for each
    chain type for IMGT. Currently combines all species (it may sometimes be
    desirable to separate species -- will consider this later).
    In target_dir, a file called 'CONSENSUS.txt' is created containing
    the consensus for each chain, while a separate .npy file for each
    chain is saved to the same directory.

    Args:
        target_dir (str): The filepath of the output directory.
        current_dir (str): The filepath of the current directory.
        alignment_fname (str): The name of the alignment file. It should
            already live in target_dir.
    """
    os.chdir(target_dir)
    combined_dict = {}
    separate_species_dict = {}

    with open(alignment_fname, "r", encoding="utf-8") as fhandle:
        read_now = False
        for line in fhandle:
            if line.startswith("#=GF"):
                read_now = True
                species, chain = line.strip().split()[-1].split("_")
                if species not in separate_species_dict:
                    separate_species_dict[species] = {}
                if chain not in combined_dict:
                    combined_dict[chain] = []
                if chain not in separate_species_dict[species]:
                    separate_species_dict[species][chain] = []
            elif line.startswith("#=GC RF"):
                read_now = False
            elif read_now:
                combined_dict[chain].append(line.strip().split()[-1])
                separate_species_dict[species][chain].append(line.strip().split()[-1])

    for chain_type, seq_list in combined_dict.items():
        consensus_list = convert_observed_to_consensus(seq_list, imgt_dp.heavy_cdrs)
        write_imgt_consensus_file(consensus_list, chain_type)
        save_consensus_array(consensus_list, chain_type)

    for species in separate_species_dict:
        for chain_type, seq_list in separate_species_dict[species].items():
            consensus_list = convert_observed_to_consensus(seq_list, imgt_dp.heavy_cdrs)
            chain_name = "_".join([species, chain_type])
            write_imgt_consensus_file(consensus_list, chain_name)
            save_consensus_array(consensus_list, chain_name)

    os.chdir(current_dir)


def load_consensus_file(consensus_file):
    """Loads a specified consensus file and stores it as a list
    of lists which can be converted to a scoring matrix."""
    consensus_list = []
    with open(consensus_file, "r") as fhandle:
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



def convert_observed_to_consensus(sequences, cdrs):
    """Converts sequences from an alignment to a consensus that
    can be written to file or converted to a scoring matrix."""
    len_distro = [len(s) for s in sequences]

    npositions = max(len_distro)
    if npositions != min(len_distro):
        raise ValueError("Sequences of different lengths encountered in the MSA")

    aa_consensus = []
    for i in range(npositions):
        position = i + 1
        if position in cdrs:
            observed_aas = ['-']
        else:
            observed_aas = set()
            for seq in sequences:
                observed_aas.add(seq[i])
            observed_aas = list(observed_aas)
        aa_consensus.append(observed_aas)
    return aa_consensus


def write_imgt_consensus_file(consensus_list, chain_type, scheme = "imgt"):
    """Writes a consensus file with a list of the amino acids observed
    at each position."""
    with open(f"{scheme.upper()}_CONSENSUS_{chain_type}.txt", "w+", encoding="utf-8") as fhandle:
        fhandle.write(f"# CHAIN {chain_type}\n")
        for i, observed_aas in enumerate(consensus_list):
            fhandle.write(f"{i+1},")
            fhandle.write(",".join(sorted(observed_aas)))
            fhandle.write("\n")
        fhandle.write("//\n\n")



def save_consensus_array(consensus_list, chain_type, constants = imgt_dp, scheme = "imgt"):
    """Converts a list of sequences for a specific chain type
    to an array with the score for each possible amino acid substitution at
    each position, including gap penalties. For IMGT (as for other numbering
    schemes), we prefer to place insertions at specific places, so we tailor
    the gap penalties to encourage this. Meanwhile, other positions are
    HIGHLY conserved, so we tailor the penalties to encourage this as well.
    IMGT numbers from 1 so we have to adjust for this."""
    blosum = substitution_matrices.load("BLOSUM62")
    blosum_key = {letter:i for i, letter in enumerate(blosum.alphabet)}


    if chain_type.endswith("K") or chain_type.endswith("L"):
        conserved_positions = constants.light_conserved_positions
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

    key_array = np.zeros((npositions, 22))

    for i, observed_aas in enumerate(consensus_list[:npositions]):
        position = i + 1

        #Choose the gap penalty
        #for template (column 20) and for query (column 21) of key array.
        if position in conserved_positions:
            key_array[i,20:] = constants.HIGHLY_CONSERVED_GAP_PENALTY
        elif position in special_positions:
            key_array[i,20] = special_positions[position][0]
            key_array[i,21:] = special_positions[position][1]
        elif position in cdrs:
            key_array[i,20:] = cdrs[position]
        else:
            key_array[i,20] = constants.DEFAULT_QUERY_GAP_PENALTY
            key_array[i,21] = constants.DEFAULT_TEMPLATE_GAP_PENALTY

        #Next, fill in the scores for other amino acid substitutions. If a conserved
        #residue, use the ones we specify here. Otherwise, use the best possible
        #score given the amino acids observed in the alignments. If the only
        #thing observed in the alignments is gaps, no penalty is applied.
        for j, letter in enumerate(allowed_aa_list):
            letter_blosum_idx = blosum_key[letter]
            if position in conserved_positions:
                if letter == conserved_positions[position]:
                    key_array[i,j] = constants.HIGHLY_CONSERVED_BONUS
            else:
                key_array[i,j] = max([blosum[letter_blosum_idx, blosum_key[k]] if
                    k != "-" else 0 for k in observed_aas])

    np.save(f"{scheme.upper()}_CONSENSUS_{chain_type}.npy", key_array)
