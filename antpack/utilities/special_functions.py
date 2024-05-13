"""Contains special functions implemented here to avoid having to
introduce additional dependencies (logsumexp from scipy and
SeqIO / translate for fasta files from Biopython)."""
import numpy as np


# This is the standard codon table, supplied explicitly here
# so that we do not need Biopython as a dependency. It has
# two modifications: it maps '...' to a gap and it includes
# the stop codons.
codon_table = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TGT': 'C', 'TGC': 'C',
        'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L',
        'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
        'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
        'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R',
        'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
        'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T',
        'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
        'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R',
        'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
        'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A',
        'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E',
        'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
        'GGG': 'G',
        'TAA': '*', 'TAG': '*', 'TGA': '*',
        '...': '.'}


def logsumexp(a, axis):
    """Compute the log of the sum of exponentials of input elements
    in a numerically stable way. This is a simplified version of
    scipy's logsumexp function.

    Args:
        a (np.ndarray): The input array.
        axis (int): The axis over which to apply this operation.

    Returns:
        out (np.ndarray): The result of the logsumexp operation.
    """
    initial_value = -np.inf if np.size(a) == 0 else None
    a_max = np.amax(a, axis=axis, keepdims=True, initial=initial_value)

    if a_max.ndim > 0:
        a_max[~np.isfinite(a_max)] = 0
    elif not np.isfinite(a_max):
        a_max = 0

    tmp = np.exp(a - a_max)

    # suppress warnings about log of zero
    with np.errstate(divide='ignore'):
        s = np.sum(tmp, axis=axis, keepdims=False)
        out = np.log(s)

    a_max = np.squeeze(a_max, axis=axis)
    out += a_max

    return out



def translate_imgt_nt(sequence, description, codon_start):
    """Translates a nucleotide sequence which MAY contain gaps
    into an amino acid sequence. This is not general-purpose and
    should be used only with the IMGT input. In particular,
    it assumes that whenever there is a trailing 1-2 extra
    nucleotides at the end of the sequence these should be
    ignored (true for IMGT input, NOT true in general).

    Args:
        sequence (str): An amino acid string.
        description (str): The IMGT description.
        codon_start (int): The position at which to start
            reading.

    Returns:
        translated_seq (str): The translated sequence. 
    """
    # Subtract 1 here, because IMGT numbers from 1.
    useq = sequence.upper()[codon_start - 1:]
    parse_len = 3 * (len(useq) // 3)
    translated_seq = []

    for i in range(0, parse_len, 3):
        codon = useq[i:i+3]
        
        # The next two statements mirror what IMGT has done when
        # translating their reference sequences.
        if 'N' in codon:
            translated_seq.append("X")
            continue
        if codon not in codon_table:
            if codon.startswith(".") or codon.endswith("."):
                translated_seq.append(".")
                continue
            # Y is unambiguous only IF either C or T at that position would
            # give the same AA. Yes, some IMGT refseqs contain Y...
            if "Y" in codon:
                if codon_table[codon.replace("Y", "T")] == \
                        codon_table[codon.replace("Y", "C")]:
                    codon = codon.replace("Y", "T")
                else:
                    raise RuntimeError("Unexpected input in the IMGT refseq db!")

            else:
                raise RuntimeError("Unexpected input in the IMGT refseq db!")
        translated_seq.append(codon_table[codon])

    return "".join(translated_seq)
