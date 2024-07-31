"""Contains the PairedChainAnnotator class, which extracts variable region
heavy and light chains from an input sequence that is presumed to contain
a single variable heavy and a single variable light chain."""
import os
import copy
import numpy as np
from .annotator_base_class import AnnotatorBaseClass
from .single_chain_annotator import SingleChainAnnotator
from antpack_cpp_ext import CTermFinder



class PairedChainAnnotator(AnnotatorBaseClass):
    """This class contains the tools needed to parse and number
    sequences that contain two variable region chains, one of type
    heavy and one of type light (i.e. typical antibodies).
    """

    def __init__(self, scheme = "imgt"):
        """Class constructor.

        Args:
            scheme (str): The numbering scheme. Must be one of "imgt",
                "martin", "kabat".

        Raises:
            RuntimeError: A RuntimeError is raised if unacceptable inputs are
                supplied.
        """
        super().__init__(scheme)
        self.light_chain_analyzer = SingleChainAnnotator(chains = ["L", "K"],
                scheme = scheme)
        self.heavy_chain_analyzer = SingleChainAnnotator(chains = ["H"],
                scheme = scheme)
        self.analyzer = SingleChainAnnotator(chains = ["H", "K", "L"],
                scheme = scheme)

        self.chain_list = ['H', 'K', 'L']

        cterm_score_array = np.zeros((9, 20, 3))
        current_dir = os.getcwd()
        try:
            project_path = os.path.abspath(os.path.dirname(__file__))
            os.chdir(os.path.join(project_path, "consensus_data"))
            for i, chain in enumerate(self.chain_list):
                npy_filename = f"CTERMFINDER_CONSENSUS_{chain}.npy"
                cterm_score_array[:,:,i] = np.load(npy_filename)

        except Exception as exc:
            os.chdir(current_dir)
            raise RuntimeError("The consensus data for the package either has been deleted or "
                    "moved or was never properly installed.") from exc

        os.chdir(current_dir)

        self.cterm_finder = CTermFinder(cterm_score_array)



    def analyze_seq(self, sequence):
        """Extracts and numbers the variable chain regions from a sequence,
        essentially by using one SingleChainAnnotator with chains set to ["K", "L"]
        and another with ["H"].

        Note that this function DOES NOT check that each sequence really does contain
        two variable chains. If you input a sequence that contains only one variable
        chain, you may get an error.

        The extracted light or heavy chains that are returned can be passed to build_msa.

        Args:
            sequence (str): A string which is a sequence
                containing the usual 20 amino acids.

        Returns:
            heavy_chain_result (tuple): A tuple of (numbering, percent_identity,
                chain_name, error_message, start, end), where start and end mark the
                points in the input sequence at which the heavy chain numbering starts
                and ends. Notice the numbering is consequently not the same length as
                the input sequnce (unlike SingleChainAnnotator). A low percent identity
                or an error message may indicate a problem with the input sequence.
            light_chain_result (tuple): A tuple of (numbering, percent_identity,
                chain_name, error_message), where start and end mark the
                points in the input sequence at which the light chain numbering starts
                and ends. A low percent identity or an error message may indicate a problem
                with the input sequence.

        Raises:
            RuntimeError: A RuntimeError is raised if invalid inputs are supplied.
        """
        if not isinstance(sequence, str):
            raise RuntimeError("This function expects a string as input.")

        mscores, mpositions = np.zeros((3)), np.zeros((3), dtype=np.int32)
        err = self.cterm_finder.find_c_terminals(sequence,
                mscores, mpositions)

        if err != "" or mscores.max() == 0:
            invalid_heavy_result = ([], 0.0, "", "Invalid sequence supplied -- "
                    "nonstandard AAs", 0, 0)
            invalid_light_result = copy.copy(invalid_heavy_result)
            return invalid_heavy_result, invalid_light_result

        # We COULD use the type of c-terminal found to try to infer
        # the chain type that occurs first, but there is always
        # the (fairly unlikely) chance that there is some very unusual
        # CDR3 sequence that is a perfect match for a C-terminal and
        # will thereby confound this algorithm. So, instead, we align
        # the region of the sequence after the first c-terminal and
        # determine the type of this, and then align the region left
        # out by this first successful alignment to whatever is
        # remaining. This is more expensive (costs 1-2 additional
        # alignments) but rules out one (unlikely) contingency.

        start_position = mpositions.min() + 9
        # Kill sequences where the number of AAs before or after start
        # is clearly insufficient to be a light or heavy chain without
        # serious deletions.
        if start_position < 85 or len(sequence) - start_position < 85:
            invalid_heavy_result = ([], 0.0, "", "Could not find c-terminal "
                "of first variable chain -- sequence may not be a paired "
                "chain sequence", 0, 0)
            invalid_light_result = copy.copy(invalid_heavy_result)
            return invalid_heavy_result, invalid_light_result

        init_results = self.analyzer.analyze_seq(sequence[start_position:])

        if init_results[2]  == "H":
            heavy_chain_result = self._pad_nterm(sequence, init_results)
            exstart = next((i for i in range(len(heavy_chain_result[0])) if
                heavy_chain_result[0][i] != '-'), 0)
            sub_results = self.light_chain_analyzer.analyze_seq(sequence[:exstart])
            light_chain_result = self._pad_cterm(sequence, sub_results)

        else:
            light_chain_result = self._pad_nterm(sequence, init_results)
            exstart = next((i for i in range(len(light_chain_result[0])) if
                light_chain_result[0][i] != '-'), 0)
            sub_results = self.heavy_chain_analyzer.analyze_seq(sequence[:exstart])
            heavy_chain_result = self._pad_cterm(sequence, sub_results)

        return heavy_chain_result, light_chain_result




    def _pad_nterm(self, sequence, aligner_results):
        """Internal-only function which preps the results of a SingleChainAnnotator
        call for return to the end user."""
        padded_num = ['-' for _ in range(len(sequence) -
            len(aligner_results[0]) )] + aligner_results[0]
        return (padded_num, aligner_results[1], aligner_results[2], aligner_results[3])

    def _pad_cterm(self, sequence, aligner_results):
        """Internal-only function which preps the results of a SingleChainAnnotator
        call for return to the end user."""
        padded_num = aligner_results[0] + ['-' for _ in range(len(sequence) -
            len(aligner_results[0])) ]
        return (padded_num, aligner_results[1], aligner_results[2], aligner_results[3])
