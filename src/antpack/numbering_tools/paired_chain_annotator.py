"""Contains the PairedChainAnnotator class, which extracts variable region
heavy and light chains from an input sequence that is presumed to contain
a single variable heavy and a single variable light chain (but may
contain only a single chain instead)."""
import os
from .cterm_finder import _load_nterm_kmers
from ..antpack_license import get_license_key_info
from antpack.antpack_cpp_ext import PairedChainAnnotatorCpp



class PairedChainAnnotator(PairedChainAnnotatorCpp):

    def __init__(self, scheme = "imgt", receptor_type="mab"):
        """Class constructor.

        Args:
            scheme (str): The numbering scheme. Must be one of "imgt",
                "martin", "kabat", "aho". If receptor_type is 'tcr'
                only "imgt" is accepted.
            receptor_type (str): One of "mab", "tcr". Default is "mab"
                (antibody).

        Raises:
            RuntimeError: A RuntimeError is raised if unacceptable inputs are
                supplied.
        """
        license_key, user_email = get_license_key_info()
        consensus_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                "consensus_data")
        kmer_dict = _load_nterm_kmers()
        super().__init__(scheme, consensus_path, kmer_dict, receptor_type,
                license_key, user_email)
