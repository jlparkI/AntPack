"""Contains the SingleChainAnnotator class, which provides the tools needed
to parse single chains (either heavy or light) or a sequence which may contain
a heavy or light chain. If you want to extract the heavy and light chain variable
regions from a sequence that contains both, use MultiChainAnnotator."""
import os
from antpack_cpp_ext import SingleChainAnnotatorCpp




class SingleChainAnnotator(SingleChainAnnotatorCpp):

    def __init__(self, chains=["H", "K", "L"], scheme="imgt",
            compress_init_gaps=False):
        """Class constructor.

        Args:
            chains (list): A list of chains. Each must be one of "H", "K", "L".
                If ["H", "K", "L"] (default), the annotator will automatically
                determine the most appropriate chain type for each input
                sequence.
            scheme (str): The numbering scheme. Must be one of "imgt",
                "martin", "kabat".
            compress_init_gaps (bool): If True, rearrange gaps in the first 5
                positions post-alignment so that gaps are at the beginning of
                the sequence wherever possible. This is more consistent with
                results from some other tools although it is debatable
                if this is more correct. Defaults to False.

        Raises:
            ValueError: A ValueError is raised if unacceptable inputs are
                supplied.
        """
        project_path = os.path.abspath(os.path.dirname(__file__))
        super().__init__(chains, scheme, compress_init_gaps,
                False, project_path)
