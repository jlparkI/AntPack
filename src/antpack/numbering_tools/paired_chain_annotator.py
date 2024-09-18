"""Contains the PairedChainAnnotator class, which extracts variable region
heavy and light chains from an input sequence that is presumed to contain
a single variable heavy and a single variable light chain."""
import os
from antpack.antpack_cpp_ext import PairedChainAnnotatorCpp



class PairedChainAnnotator(PairedChainAnnotatorCpp):

    def __init__(self, scheme = "imgt"):
        """Class constructor.

        Args:
            scheme (str): The numbering scheme. Must be one of "imgt",
                "martin", "kabat".

        Raises:
            RuntimeError: A RuntimeError is raised if unacceptable inputs are
                supplied.
        """
        consensus_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                "consensus_data")
        super().__init__(scheme, consensus_path)