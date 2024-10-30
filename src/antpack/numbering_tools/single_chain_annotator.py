"""Contains the SingleChainAnnotator class, which provides the tools needed
to parse single chains (either heavy or light) or a sequence which likely contains
a heavy or light chain."""
import os
from .cterm_finder import _load_nterm_kmers
from antpack.antpack_cpp_ext import SingleChainAnnotatorCpp




class SingleChainAnnotator(SingleChainAnnotatorCpp):

    def __init__(self, chains=["H", "K", "L"], scheme="imgt"):
        """Class constructor.

        Args:
            chains (list): A list of chains. Each must be one of "H", "K", "L".
                If ["H", "K", "L"] (default), the annotator will automatically
                determine the most appropriate chain type for each input
                sequence.
            scheme (str): The numbering scheme. Must be one of "imgt",
                "martin", "kabat", "aho".

        Raises:
            ValueError: A ValueError is raised if unacceptable inputs are
                supplied.
        """
        consensus_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                "consensus_data")
        kmer_dict = _load_nterm_kmers()
        super().__init__(chains, scheme, consensus_path, kmer_dict)
