"""Contains the SingleChainAnnotator class, which provides the tools needed
to parse single chains (either heavy or light) or a sequence which may contain
a heavy or light chain. If you want to extract the heavy and light chain variable
regions from a sequence that contains both, use MultiChainAnnotator."""
import os
import numpy as np
from .constants import imgt_default_params, martin_default_params, kabat_default_params
from .annotator_base_class import AnnotatorBaseClass
from antpack_cpp_ext import validate_sequence, IGAligner




class SingleChainAnnotator(AnnotatorBaseClass):

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
        super().__init__(scheme)

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
        self.score_mat_sizes = []

        try:
            os.chdir(os.path.join(project_path, "consensus_data"))
            for chain in chains:
                npy_filename = f"{scheme.upper()}_CONSENSUS_{chain}.npy"
                text_file = f"{scheme.upper()}_CONSENSUS_{chain}.txt"
                chain_name = chain

                score_matrix = np.load(npy_filename)
                con_map = self._load_consensus_map(text_file)
                self.scoring_tools.append(IGAligner(score_matrix, con_map,
                                chain_name, scheme, defaults.DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                                defaults.DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                                compress_init_gaps))
                self.score_mat_sizes.append(score_matrix.shape[0])

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



    def analyze_seqs(self, sequences, get_region_labels = False):
        """Numbers and scores a list of input sequences. The outputs
        can be passed to build_msa if desired.

        Args:
            sequences (list): A list of strings, each of which is a sequence
                containing the usual 20 amino acids.
            get_region_labels (bool): If True, get a list of labels: "-", "fmwk1",
                "cdr1", "fmwk2", "cdr2", "fmwk3", "cdr3", "fmwk4" to indicate
                to which region each numbered amino acid belongs. The cdr definitions
                that are used are the same as those for the numbering scheme (i.e.
                if using IMGT numbering IMGT CDR definitions are used). If False
                this list is not generated.

        Returns:
            sequence_results (list): A list of tuples of (sequence_numbering, percent_identity,
                chain_name, error_message) if get_region_labels is False or
                (sequence_numbering, percent_identity, chain_name, error_message, region_labels)
                if get_region_labels is True. If no error was encountered, the error
                message is "". An alignment with low percent identity (e.g. < 0.85)
                may indicate a sequence that is not really an antibody, that contains
                a large deletion, or is not of the selected chain type.
        """
        if not isinstance(sequences, list):
            raise ValueError("sequences should be a list of strings.")


        sequence_results = []
        for sequence in sequences:
            if not validate_sequence(sequence):
                sequence_results.append(([], 0.0, "",
                    "Invalid sequence supplied -- nonstandard AAs"))
                continue
            results = [scoring_tool.align(sequence, get_region_labels) for
                    scoring_tool in self.scoring_tools]
            results = sorted(results, key=lambda x: x[1])
            if get_region_labels:
                sequence_results.append(results[-1])
            else:
                sequence_results.append(results[-1][:-1])

        return sequence_results


    def analyze_seq(self, sequence, get_region_labels = False):
        """Numbers and scores a single input sequence. A list of
        outputs from this function can be passed to build_msa
        if desired.

        Args:
            sequence (str): A string which is a sequence containing the usual 20 amino acids.
            get_region_labels (bool): If True, get a list of labels: "-", "fmwk1",
                "cdr1", "fmwk2", "cdr2", "fmwk3", "cdr3", "fmwk4" to indicate
                to which region each numbered amino acid belongs. The cdr definitions
                that are used are the same as those for the numbering scheme (i.e.
                if using IMGT numbering IMGT CDR definitions are used). If False
                this list is not generated.

        Returns:
            sequence_results (tuple): A tuple of (sequence_numbering, percent_identity,
                chain_name, error_message) if get_region_labels is False or
                (sequence_numbering, percent_identity, chain_name, error_message, region_labels)
                if get_region_labels is True. If no error was encountered, the error
                message is "". An alignment with low percent identity (e.g. < 0.85)
                may indicate a sequence that is not really an antibody, that contains
                a large deletion, or is not of the selected chain type.
        """
        if not validate_sequence(sequence):
            return ([], 0.0, "", "Invalid sequence supplied -- nonstandard AAs")

        results = [scoring_tool.align(sequence, get_region_labels) for
                scoring_tool in self.scoring_tools]
        results = sorted(results, key=lambda x: x[1])
        sequence_results = results[-1]
        if get_region_labels:
            return sequence_results
        return sequence_results[:-1]




    def _alignment_test(self, sequence):
        """Numbers and scores a single input sequence. Also returns the raw scoring
        matrix. This should therefore be used for testing only since this
        is not useful to typical end users and is slower than the analyze_seq(s)
        functions.

        Args:
            sequence (str): A string which is a sequence containing the usual 20 amino acids.

        Returns:
            sequence_results (tuple): A tuple of (sequence numbering, percent_identity,
                chain_name, error_message, score_matrix, path_trace).
                If no error was encountered, the error message is "".
        """
        results = []

        for score_mat_size, scoring_tool in zip(self.score_mat_sizes, self.scoring_tools):
            score_mat = np.zeros((score_mat_size + 1, len(sequence) + 1))
            path_trace = np.zeros((score_mat_size + 1, len(sequence) + 1), dtype=np.uint8)
            output = scoring_tool.align_test_only(sequence, False, score_mat, path_trace)
            output_tuple = (output[0], output[1], output[2], output[3], score_mat, path_trace)
            results.append(output_tuple)

        results = sorted(results, key=lambda x: x[1])[-1]
        return results
