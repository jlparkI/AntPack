"""Contains the PairedChainAnnotator class, which extracts variable region
heavy and light chains from an input sequence that is presumed to contain
a single variable heavy and a single variable light chain (but may
contain only a single chain instead)."""
import os
from .cterm_finder import _load_nterm_kmers
from antpack.antpack_cpp_ext import PairedChainAnnotatorCpp



class PairedChainAnnotator(PairedChainAnnotatorCpp):

    def __init__(self, scheme = "imgt"):
        """Class constructor.

        Args:
            scheme (str): The numbering scheme. Must be one of "imgt",
                "martin", "kabat", "aho".

        Raises:
            RuntimeError: A RuntimeError is raised if unacceptable inputs are
                supplied.
        """
        consensus_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                "consensus_data")
        kmer_dict = _load_nterm_kmers()
        super().__init__(scheme, consensus_path, kmer_dict)
