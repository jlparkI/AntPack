"""Contains the SingleChainAnnotator class, which provides the tools needed
to parse antibody sequences of known chain type. A separate class is maintained
to handle sequences of unknown chain type or that may contain multiple
chains. SingleChainAnnotator is however substantially faster since it does
not need to determine the type of chain or whether multiple chains are
present, so it is useful when you have a large number of sequences of
known chain type to process."""
import os
import numpy as np
from Bio import SeqIO
from .constants import allowed_inputs
from ant_ext import IMGTAligner, validate_sequence







class SingleChainAnnotator:
    """This class contains the tools needed to parse and number
    either a single sequence, a list of sequences, or a fasta file containing
    many sequences that belong to one of several possible chain types. If your
    sequences may contain multiple chains e.g. both heavy and
    light chains (e.g. an scFv with linker), use MultiChainAnnotator instead.
    """

    def __init__(self, species = "all", chains = ["H"]):
        """Class constructor.

        Args:
            species (str): Must be one of "alpaca", "cow", "human",
                "mouse", "pig", "rabbit", "rat", "rhesus", "all".
            chain (list): A list of chains. Each must be one of "H", "K", "L",
                "KL", "A", "B", "D", "G". Note that not all chains are
                supported for all species. "A", "B", "D", "G" for example
                are supported for human, mouse and all but not for rabbit,
                rat, rhesus or pig.

        Raises:
            ValueError: A ValueError is raised if unacceptable inputs are
                supplied.
        """
        if species not in ["alpaca", "cow", "human", "mouse", "pig", "rabbit",
                "rat", "rhesus", "all"]:
            raise ValueError(f"Unrecognized species {species} supplied.")
        if len(chains) == 0:
            raise ValueError("Must supply at least one chain.")
        for chain in chains:
            if chain not in ["H", "K", "L", "A", "B", "D", "G"]:
                raise ValueError(f"Unrecognized chain {chain} supplied.")
            if chain in ["A", "B", "D", "G"] and species not in ["human", "mouse"]:
                raise ValueError(f"Unsupported chain-species combo {chain} {species} supplied.")
            if chain in ["H"] and species in ["rat"]:
                raise ValueError(f"Unsupported chain-species combo {chain} {species} supplied.")

        self.chains = chains

        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        score_matrices, consensus_maps = [], []

        try:
            os.chdir(os.path.join(project_path, "consensus_data"))
            for chain in chains:
                if species == "all":
                    npy_filename = f"CONSENSUS_{chain}.npy"
                    text_file = f"CONSENSUS_{chain}.txt"
                else:
                    npy_filename = f"CONSENSUS_{species}_{chain}.npy"
                    text_file = f"CONSENSUS_{species}_{chain}.txt"
                score_matrices.append(np.load(npy_filename))
                consensus_maps.append(self._load_consensus_map(text_file))
        except Exception as exc:
            os.chdir(current_dir)
            raise ValueError("The consensus data for the package either has been deleted or "
                    "moved or was never properly installed.") from exc

        os.chdir(current_dir)
        self.scoring_tools = []

        for (score_matrix, con_map) in zip(score_matrices, consensus_maps):
            if len(score_matrix.shape) != 2:
                raise ValueError("A score matrix was located but has an unexpected shape. "
                        "Please report this error to the package maintainer.")
            if score_matrix.shape[0] != 128 or score_matrix.shape[1] != 21:
                raise ValueError("A score matrix was located but has an unexpected shape. "
                        "Please report this error to the package maintainer.")

            try:
                self.scoring_tools.append(IMGTAligner(score_matrix, con_map))
            except:
                import pdb
                pdb.set_trace()



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
        # Check to make sure all expected IMGT positions are found.
        expected_imgt_positions = {str(i):0 for i in range(1,129)}

        with open(filename, "r", encoding="utf-8") as fhandle:
            read_now = False
            consensus = []
            last_position = 0
            for line in fhandle:
                if line.startswith("# CHAIN"):
                    read_now = True
                    continue
                if line.startswith("//"):
                    break
                if not read_now:
                    continue

                position = line.split(",")[0]
                if position not in expected_imgt_positions:
                    raise RuntimeError("Problem with consensus file.")
                if int(position) - 1 != last_position:
                    raise RuntimeError("Problem with consensus file.")
                last_position += 1
                expected_imgt_positions[position] += 1
                aas = line.strip().split(",")[1:]
                if "-" in aas:
                    consensus.append([])
                    continue
                consensus.append(aas)
        return consensus



    def analyze_online_seqs(self, sequences, scheme="imgt"):
        """Numbers and scores a list of input sequences.

        Args:
            sequences (list): A list of strings, each of which is a sequence
                containing the usual 20 amino acids.
            scheme (str): The numbering scheme. Should be one of 'imgt', 'martin',
                'chothia', 'kabat', 'aho', 'wolfguy'.

        Returns:
            sequence_results (list): A list of tuples of (sequence numbering, score,
                None). If a sequence is invalid, the tuple is (None, None, error message),
                where error message is a string. The results are in the same order
                as the inputs.
        """
        if scheme not in allowed_inputs.allowed_schemes:
            raise ValueError(f"Input scheme {scheme} is not allowed; "
                "should be one of {allowed_inputs.allowed_schemes}")
        if not isinstance(sequences, list):
            raise ValueError("sequences should be a list of strings.")


        sequence_results = []
        for sequence in sequences:
            if not validate_sequence(sequence):
                sequence_results.append((None, None, "invalid_sequence"))
                continue
            results = []
            for (chain, scoring_tool) in self.scoring_tools:
                numbering, percent_ident, err_code = scoring_tool.align(sequence)
                results.append((numbering, percent_ident, chain, err_code))
            results = sorted(results, key=lambda x: x[1])
            numbering, percent_ident, chain, err_code = results[-1]
            if err_code == 0:
                sequence_results.append((None, None, None, "invalid_sequence"))
            elif err_code == 2 or len(sequence) != len(numbering):
                sequence_results.append((None, None, None, "fatal_alignment_error"))
            else:
                sequence_results.append((numbering, percent_ident, chain, None))

        return sequence_results


    def analyze_fasta(self, fasta_file, scheme="imgt"):
        """A generator that numbers and scores sequences from a fasta
        file. Since it is a generator it will only load one sequence at
        a time. You should in your code when retrieving results from
        this generator ensure that the error code for each sequence is None
        before doing anything else with the results.

        Args:
            fasta_file (str): A filepath to a valid fasta file.
            scheme (str): The numbering scheme. Should be one of 'imgt', 'martin',
                'chothia', 'kabat', 'aho', 'wolfguy'.

        Returns:
            sequence_results (list): A list of tuples of (sequence numbering, score,
                None). If a sequence is invalid, the tuple is (None, None, error message),
                where error message is a string. The results are in the same order
                as the inputs.

        Raises:
            ValueError: A ValueError is raised if unacceptable inputs are supplied.
        """
        if scheme not in allowed_inputs.allowed_schemes:
            raise ValueError(f"Input scheme {scheme} is not allowed; "
                "should be one of {allowed_inputs.allowed_schemes}")
        if not os.path.isfile(fasta_file):
            raise ValueError("Nonexistent file path supplied.")

        with open(fasta_file, "r", encoding="utf-8") as fhandle:
            for seqrecord in SeqIO.parse(fhandle, "fasta"):
                sequence = str(seqrecord.seq)
                if not validate_sequence(sequence):
                    best_result = (None, None, "invalid_sequence")
                else:
                    results = []
                    for (chain, scoring_tool) in zip(self.chains, self.scoring_tools):
                        numbering, percent_ident, err_code = scoring_tool.align(sequence)
                        results.append((numbering, percent_ident, chain, err_code))
                    results = sorted(results, key=lambda x: x[1])
                    numbering, percent_ident, chain, err_code = results[-1]
                    if err_code != 1:
                        best_result = (None, None, None, "alignment_error")
                    else:
                        best_result = (numbering, percent_ident, chain, None)
                yield seqrecord, best_result
