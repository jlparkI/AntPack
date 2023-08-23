import os
import sys
import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices

#Note that AAs in this list should be in alphabetical order --
#very important.
standard_aas = ["A", "C", "D", "E", "F", "G", "H", "I", "K",
        "L", "M", "N", "P", "Q", "R", "S", "T",
        "V", "W", "Y"]


def main():
    array_key_dict = dict()
    with open(sys.argv[1], "r", encoding="utf-8") as fhandle:
        read_now = False
        for line in fhandle:
            if line.startswith("#=GF"):
                read_now = True
                species_chain = line.strip().split()[-1]
                sequences = []
            elif line.startswith("#=GC RF"):
                read_now = False
                if species_chain.endswith("H"):
                    process_save_heavy(species_chain, sequences)
                sequences = []
            elif read_now:
                sequences.append(line.strip().split()[-1])


def process_save_heavy(species_chain, sequences):
    """Converts a list of sequences for a specific species-chain combination
    to an array with the score for each possible amino acid substitution at
    each position, including gap penalties. For IMGT (as for other numbering
    schemes), we prefer to place insertions at specific places, so we tailor
    the gap penalties to encourage this. Meanwhile, other positions are
    HIGHLY conserved, so we tailor the penalties to encourage this as well.
    IMGT numbers from 1 so we have to adjust for this."""
    blosum = substitution_matrices.load("BLOSUM62")
    blosum_key = {letter:i for i, letter in enumerate(blosum.alphabet)}

    INSERTION_POSITIONS = {33, 61, 111}
    CONSERVED_POSITIONS = {23:["C"], 41:["W"], 104:["C"], 118:["F", "W"]}
    cdrs = {27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38,
            56, 57, 58, 59, 60, 61, 62, 63, 64, 65,
            105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117}
    key_array = np.zeros((128, 21))
    len_distro = [len(s) for s in sequences]

    if max(len_distro) != min(len_distro):
        raise ValueError("Sequences of different lengths encountered in the MSA")

    for i in range(max(len_distro)):
        position = i + 1
        observed_aas = set()
        for s in sequences:
            observed_aas.add(s[i])
        observed_aas = list(observed_aas)

        #First, choose the w2 weight which increases the cost of deviating
        #from expected at conserved positions. Then, choose the gap penalty
        #(column 20 of key array)
        w2 = 1.0
        if position in CONSERVED_POSITIONS:
            w2 = 5.0
            key_array[i,20] = -55.0
        elif position in cdrs and position not in INSERTION_POSITIONS:
            key_array[i,20] = -26.0
        elif position in INSERTION_POSITIONS:
            key_array[i,20] = -1.0
        elif position in [1,128]:
            key_array[i,20] = 0.0
        else:
            key_array[i,20] = -26.0

        #Next, fill in the scores for other amino acid substitutions. If a conserved
        #residue, use the ones we specify here. Otherwise, use the best possible
        #score given the amino acids observed in the alignments. If the only
        #thing observed in the alignments is gaps, no penalty is applied.
        for j, letter in enumerate(standard_aas):
            letter_blosum_idx = blosum_key[letter]
            if position in CONSERVED_POSITIONS:
                score = max([blosum[letter_blosum_idx, blosum_key[k]] for k in
                    CONSERVED_POSITIONS[position]])
            else:
                score = max([blosum[letter_blosum_idx, blosum_key[k]] if k != "-" else 0 for k in
                    observed_aas])

            key_array[i,j] = score * w2



if __name__ == "__main__":
    main()
