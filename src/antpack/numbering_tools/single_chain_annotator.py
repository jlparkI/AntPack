"""Contains the SingleChainAnnotator class, which provides the tools needed
to parse single chains (either heavy or light) or a sequence which likely contains
a heavy or light chain."""
import os
from .cterm_finder import _load_nterm_kmers
from ..antpack_license import get_license_key_info
from antpack.antpack_cpp_ext import SingleChainAnnotatorCpp




class SingleChainAnnotator(SingleChainAnnotatorCpp):

    def __init__(self, chains=["H", "K", "L"], scheme="imgt"):
        """Class constructor.

        Args:
            chains (list): A list of chains. Each must EITHER be one of
                "H", "K", "L" for antibodies or one of "A", "B", "D", "G"
                for TCRs. If ["H", "K", "L"] (default) or ["A", "B", "D", "G"]
                the annotator will automatically determine the most appropriate
                chain type for each input sequence. You cannot supply a mixture
                of TCR and antibody chains (e.g. ["H", "A"]) -- the list you
                supply must contain either TCR or antibody chains but not
                both.
            scheme (str): The numbering scheme. Must be one of "imgt",
                "martin", "kabat", "aho". If TCR chains are supplied,
                only "imgt" is accepted.

        Raises:
            ValueError: A ValueError is raised if unacceptable inputs are
                supplied.
        """
        license_key, user_email = get_license_key_info()
        consensus_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                "consensus_data")
        kmer_dict = _load_nterm_kmers()
        super().__init__(chains, scheme, consensus_path, kmer_dict,
                license_key, user_email)
