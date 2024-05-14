"""This is the base class for Single and MultiChainAnnotator. Contains routines
shared by both."""
import os
import numpy as np
from .constants import imgt_default_params, martin_default_params, kabat_default_params
from antpack_cpp_ext import IGAligner




class AnnotatorBaseClass:
    """This class contains the tools needed to parse and number
    sequences and is the parent for both Multi and SingleChain
    Annotators.
    """

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
        if len(chains) == 0:
            raise ValueError("Must supply at least one chain.")
        for chain in chains:
            if chain not in ["H", "K", "L"]:
                raise ValueError(f"Unrecognized chain {chain} supplied.")

        if scheme not in ["imgt", "martin", "kabat"]:
            raise ValueError("Unsupported scheme supplied.")

        if scheme == "martin":
            defaults = martin_default_params
        elif scheme == "kabat":
            defaults = kabat_default_params
        else:
            defaults = imgt_default_params

        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        self.scoring_tools = []

        try:
            os.chdir(os.path.join(project_path, "consensus_data"))
            for chain in chains:
                npy_filename = f"{scheme.upper()}_CONSENSUS_{chain}.npy"
                text_file = f"{scheme.upper()}_CONSENSUS_{chain}.txt"
                chain_name = chain

                score_matrix = np.load(npy_filename)
                con_map = self._load_consensus_map(text_file)
                # Note that IGAligner class constructor checks the input score
                # matrix to ensure the right dimensions; if it does not like these,
                # it will throw an exception that the PyBind wrapper will hand
                # off to Python.
                self.scoring_tools.append(IGAligner(score_matrix, con_map,
                                chain_name, scheme, defaults.DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                                defaults.DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                                compress_init_gaps))
        except Exception as exc:
            os.chdir(current_dir)
            raise RuntimeError("The consensus data for the package either has been deleted or "
                    "moved or was never properly installed.") from exc

        os.chdir(current_dir)



    def _load_consensus_map(self, filename):
        """A convenience function that loads the consensus map under
        filename and converts it to a list of lists. Each entry in the
        list is a list of amino acids allowed at that position. The
        resulting structure is passed to the C++ extension. Note
        that an empty list at a given position indicates any AA is
        tolerated at that position. Therefore, if a '-' is found
        at a given position, assume any AA is tolerated there.

        Args:
            filename (str): Path to a valid consensus sequence file.
        """
        with open(filename, "r", encoding="utf-8") as fhandle:
            read_now = False
            consensus = []
            last_position = 0
            for line in fhandle:
                if line.startswith("#"):
                    read_now = True
                    continue
                if line.startswith("//"):
                    break
                if not read_now:
                    continue

                position = line.split(",")[0]
                if int(position) - 1 != last_position:
                    raise RuntimeError("Problem with consensus file.")
                last_position += 1
                aas = line.strip().split(",")[1:]
                if "-" in aas:
                    consensus.append([])
                    continue
                consensus.append(aas)
        return consensus
