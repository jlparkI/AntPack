"""Contains tools for translating DNA sequence to amino acid
sequence, especially in cases where the reading frame is
unknown, preparatory to AA numbering."""
import os
import gzip
from antpack.antpack_cpp_ext import DNASeqTranslatorCpp



class DNASeqTranslator(DNASeqTranslatorCpp):
    """Contains functions that determine the correct reading
    frame and complement to use for an input DNA sequence
    and translates the input DNA sequence to protein / AA."""

    def __init__(self):
        common_kmers = _load_common_kmers()
        super().__init__(common_kmers)



def _load_common_kmers():
    """A convenience function for loading a set of kmers frequently
    observed in heavy and light chains."""
    filepath = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                "consensus_data", "mabs", "ALL_REGION_KMER_FILES.txt.gz")
    with gzip.open(filepath, "rt") as fhandle:
        common_kmers = {line.strip().split(',')[1]
                for line in fhandle}
    return common_kmers
