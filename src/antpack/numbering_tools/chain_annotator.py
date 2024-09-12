"""Contains the ChainAnnotator class, which extracts variable region
heavy and light chains from an input sequence that contains an
unknown number of heavy and light chains."""
import os
from antpack_cpp_ext import ChainAnnotatorCpp



class ChainAnnotator(ChainAnnotatorCpp):

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
