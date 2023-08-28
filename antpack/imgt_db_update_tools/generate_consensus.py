"""Contains the tools needed to convert a stockholm alignment of
all the chain types into a set of .npy arrays and a CONSENSUS.txt
file with consensus sequences for each chain type. The .npy arrays
are used to score new sequences so that they are correctly
aligned and numbered."""
import os
import numpy as np
from Bio.Align import substitution_matrices
from ..constants.allowed_inputs import allowed_aa_list
from ..constants.hmmbuild_constants import insertion_positions, conserved_positions, cdrs


def build_consensus_files(target_dir, current_dir, alignment_fname):
    """Builds a consensus for the amino acids at each position for each
    chain type. Currently combines all species (it may sometimes be
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
    separate_species_dict = {}
    combined_dict = {}

    with open(alignment_fname, "r", encoding="utf-8") as fhandle:
        read_now = False
        for line in fhandle:
            if line.startswith("#=GF"):
                read_now = True
                species, chain = line.strip().split()[-1].split("_")
                if chain not in combined_dict:
                    separate_species_dict[chain] = {}
                    combined_dict[chain] = []
                if species not in separate_species_dict[chain]:
                    separate_species_dict[chain][species] = []
            elif line.startswith("#=GC RF"):
                read_now = False
            elif read_now:
                separate_species_dict[chain][species].append(line.strip().split()[-1])
                combined_dict[chain].append(line.strip().split()[-1])

    for chain_type, seq_list in combined_dict.items():
        write_consensus_file(seq_list, chain_type)
        save_consensus_array(seq_list, chain_type)

    for chain_type in separate_species_dict.keys():
        for species, seq_list in separate_species_dict[chain_type].items():
            output_name = "_".join([species, chain_type])
            save_consensus_array(seq_list, output_name)
            write_consensus_file(seq_list, output_name)

    os.chdir(current_dir)


def write_consensus_file(sequences, chain_type):
    """Writes a consensus file with a list of the amino acids observed
    at each position."""
    with open(f"CONSENSUS_{chain_type}.txt", "w+", encoding="utf-8") as fhandle:
        fhandle.write(f"# CHAIN {chain_type}\n")
        position_key = {i:set() for i in range(len(sequences[0]))}
        for sequence in sequences:
            for i, letter in enumerate(sequence):
                position_key[i].add(letter)
        for i in range(len(sequences[0])):
            observed_aas = sorted(list(position_key[i]))
            fhandle.write(f"{i+1},")
            fhandle.write(",".join(observed_aas))
            fhandle.write("\n")
        fhandle.write("//\n\n")



def save_consensus_array(sequences, chain_type):
    """Converts a list of sequences for a specific chain type
    to an array with the score for each possible amino acid substitution at
    each position, including gap penalties. For IMGT (as for other numbering
    schemes), we prefer to place insertions at specific places, so we tailor
    the gap penalties to encourage this. Meanwhile, other positions are
    HIGHLY conserved, so we tailor the penalties to encourage this as well.
    IMGT numbers from 1 so we have to adjust for this."""
    blosum = substitution_matrices.load("BLOSUM62")
    blosum_key = {letter:i for i, letter in enumerate(blosum.alphabet)}

    key_array = np.zeros((128, 21))
    len_distro = [len(s) for s in sequences]

    if max(len_distro) != min(len_distro):
        raise ValueError("Sequences of different lengths encountered in the MSA")

    for i in range(len_distro[0]):
        position = i + 1
        observed_aas = set()
        for seq in sequences:
            observed_aas.add(seq[i])
        observed_aas = list(observed_aas)

        #First, choose the w2 weight which increases the cost of deviating
        #from expected at conserved positions. Then, choose the gap penalty
        #(column 20 of key array)
        score_weight = 1.0
        if position in conserved_positions:
            score_weight = 5.0
            key_array[i,20] = -55.0
        elif position in cdrs and position not in insertion_positions:
            key_array[i,20] = cdrs[position]
        elif position in insertion_positions:
            key_array[i,20] = -1.0
        elif position in [1,128]:
            key_array[i,20] = 0.0
        else:
            key_array[i,20] = -26.0

        #Next, fill in the scores for other amino acid substitutions. If a conserved
        #residue, use the ones we specify here. Otherwise, use the best possible
        #score given the amino acids observed in the alignments. If the only
        #thing observed in the alignments is gaps, no penalty is applied.
        for j, letter in enumerate(allowed_aa_list):
            letter_blosum_idx = blosum_key[letter]
            if position in conserved_positions:
                score = max([blosum[letter_blosum_idx, blosum_key[k]] for k in
                    conserved_positions[position]])
            else:
                score = max([blosum[letter_blosum_idx, blosum_key[k]] if k != "-" else 0 for k in
                    observed_aas])

            key_array[i,j] = score * score_weight

    np.save(f"CONSENSUS_{chain_type}.npy", key_array)
