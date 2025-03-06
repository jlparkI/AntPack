"""Contains the CTermFinder class. Used for testing only."""
import os
import gzip
from antpack.antpack_cpp_ext import PrefilteringTool



class PyCTermFinder(PrefilteringTool):

    def __init__(self):
        """Class constructor.
        """
        consensus_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                "consensus_data")
        kmer_dict = _load_nterm_kmers()
        super().__init__(consensus_path, kmer_dict)



def _load_nterm_kmers():
    """A convenience function for loading a dictionary of n-terminal
    kmers."""
    filepath = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                "consensus_data", "mabs", "START_REGION_KMER_FILES.txt.gz")

    kmer_dict = {}
    kmers_to_ints = {"H":0, "K":1, "L":2}

    with gzip.open(filepath, "rt") as fhandle:
        for line in fhandle:
            chain, kmer = line.strip().split(",")
            kmer_dict[kmer] = kmers_to_ints[chain]
    return kmer_dict
