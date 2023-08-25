"""Contains the SingleChainAnnotator class, which provides the tools needed
to parse antibody sequences of known chain type. A separate class is maintained
to handle sequences of unknown chain type or that may contain multiple
chains. SingleChainAnnotator is however substantially faster since it does
not need to determine the type of chain or whether multiple chains are
present, so it is useful when you have a large number of sequences of
known chain type to process."""
import os
import numpy as np
from .constants import allowed_inputs
from ant_ext import IMGTAligner, validate_sequence







class SingleChainAnnotator:
    """This class contains the tools needed to parse and number
    either a single sequence, a list of sequences, or a fasta file containing
    many sequences that belong to a single known chain type. If you do not
    know the chain type, or if your sequences contain e.g. both heavy and
    light chains (e.g. an scFv with linker), use MultiChainAnnotator instead.
    """

    def __init__(self, species = "all", chain = "H"):
        """Class constructor.

        Args:
            species (str): Must be one of "alpaca", "cow", "human",
                "mouse", "pig", "rabbit", "rat", "rhesus", "all".
            chain (str): Must be one of "H", "K", "L", "KL", "A", "B", "D", "G".
                Note that not all chains are supported for all species. "A",
                "B", "D", "G" for example are supported for human, mouse and all
                but not for rabbit, rat, rhesus or pig.
        """
        if species not in ["alpaca", "cow", "human", "mouse", "pig", "rabbit",
                "rat", "rhesus", "all"]:
            raise ValueError(f"Unrecognized species {species} supplied.")
        if chain not in ["H", "K", "L", "A", "B", "D", "G"]:
            raise ValueError(f"Unrecognized chain {chain} supplied.")
        if chain in ["A", "B", "D", "G"] and species not in ["human", "mouse"]:
            raise ValueError(f"Unsupported chain-species combo {chain} {species} supplied.")
        if chain in ["H"] and species in ["rat"]:
            raise ValueError(f"Unsupported chain-species combo {chain} {species} supplied.")

        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        try:
            os.chdir(os.path.join(project_path, "consensus_data"))
            if species == "all":
                filename = f"CONSENSUS_{chain}.npy"
            else:
                filename = f"CONSENSUS_{species}_{chain}.npy"
            self.score_matrix = np.load(filename)
        except Exception as exc:
            os.chdir(current_dir)
            raise ValueError("The consensus data for the package either has been deleted or "
                    "moved or was never properly installed.") from exc

        os.chdir(current_dir)
        if len(self.score_matrix) != 2:
            raise ValueError("The score matrix was located but has an unexpected shape. "
                    "Please report this error to the package maintainer.")
        if self.score_matrix.shape[0] != 128 or self.score_matrix.shape[1] != 21:
            raise ValueError("The score matrix was located but has an unexpected shape. "
                    "Please report this error to the package maintainer.")



    def _in_memory_checks(self, sequences, scheme):
        """Checks user inputs for in-memory sequence analysis.

        Args:
            sequences (list): A list of sequence strings.
            sequence_names (list): A list of names for each sequence.
            scheme (str): The numbering scheme. Should be one of 'imgt', 'martin',
                'chothia', 'kabat', 'aho', 'wolfguy'.

        Raises:
            ValueError: A ValueError is raised if inappropriate inputs are detected.
        """
        if scheme not in allowed_inputs.allowed_schemes:
            raise ValueError(f"Input scheme {scheme} is not allowed; "
                "should be one of {allowed_inputs.allowed_schemes}")
        if not isinstance(sequences, list):
            raise ValueError("sequences should be a list of strings.")
        for sequence in sequences:
            if not validate_sequence(sequence):
                raise ValueError(f"Sequence {sequence} contains nonstandard "
                        "amino acids.")


    def analyze_online_seqs(self, sequences, sequence_names, scheme="imgt",
            get_all_scores=False, assign_germline=False, bit_score_threshold=80,
            ncpu=0):
        """Numbers and scores a list of input sequences.

        Args:
            sequences (list): A list of strings, each of which is a sequence
                containing the usual 20 amino acids.
            sequence_names (list): A list of strings, each of which is the
                name of the corresponding sequence in sequences.
            scheme (str): The numbering scheme. Should be one of 'imgt', 'martin',
                'chothia', 'kabat', 'aho', 'wolfguy'.
            get_all_scores (bool): If True, in addition to numbering the top hit for
                each sequence, bitscores are returned for all hits rather than just the
                top one. This is useful in instances where we would like to know for
                example what the bitscore is against mouse even though the sequence
                is in fact human.
            assign_germline (bool): If True, indicate the germline genes that most
                closely correspond to the search sequences.
            bit_score_threshold (float): The bit score threshold for a hit to
                be considered a hit.
            ncpu (int): The number of threads hmmer is allowed to use. If 0, it will
                be autoselected.
        """


        self._in_memory_checks(sequences, sequence_names, scheme)

        loaded_seqs = [TextSequence(name=name.encode(), sequence=seq).digitize(self.alphabet)
                for (name, seq) in zip(sequence_names, sequences)]

        all_hits = list(hmmscan(self.hmm_db, loaded_seqs, cpus=ncpu))
        alignments, score_list = self._parse_hmmer_hits(all_hits, scheme, bit_score_threshold,
                                get_all_scores, assign_germline)

        if get_all_scores:
            return alignments, score_list
        return alignments
