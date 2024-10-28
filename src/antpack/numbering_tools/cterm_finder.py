"""Contains the CTermFinder class."""
import os
from antpack.antpack_cpp_ext import CTermFinder



class PyCTermFinder(CTermFinder):

    def __init__(self):
        """Class constructor.
        """
        consensus_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                "consensus_data")
        super().__init__(consensus_path)
