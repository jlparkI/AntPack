"""Provides tools for scoring an antibody sequence against germlines of
different species to determine which is the most likely species of origin."""
import pyhmmer
from pyhmmer.hmmer import hmmscan
from pyhmmer.easel import TextSequence

from .constants import allowed_inputs









class AntibodyScoring:
    """This class contains the tools needed to score an antibody against
    profiles representing different species (and thereby determine the
    most likely species of origin).

    Attributes:
        hmm_db: A pyhmmer optimized profile.
        alphabet: A pyhmmer alphabet (generally standard AAs).
    """

    def __init__(self, hmm_path):
        """Class constructor.

        Args:
            hmm_path (str): The filepath for an HMM profile (if you don't have
                a profile, build one using IMGTBuilder).
        """
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            self.hmm_db = hmm_file.optimized_profiles()
        self.alphabet = pyhmmer.easel.Alphabet.amino()


    def _validate_sequence(self, sequence):
        """
        Check whether a sequence is a protein sequence or is exceptionally long / contains
        nonstandard characters.
        """
        if len(sequence) > 1000:
            raise ValueError(f"The submitted sequence {sequence} is too long "
                "to be an antibody sequence.")
        if len([s for s in sequence if s not in allowed_inputs.allowed_amino_acids]) > 0:
            raise ValueError(f"Submitted sequence {sequence} contains non-amino-acid characters.")


    def _in_memory_checks(self, sequences, sequence_names):
        """Checks user inputs for in-memory sequence analysis.

        Args:
            sequences (list): A list of sequence strings.
            sequence_names (list): A list of names for each sequence.

        Raises:
            ValueError: A ValueError is raised if inappropriate inputs are detected.
        """
        if self.hmm_db is None:
            raise ValueError("An hmm profile has not been added! Consider using "
                    "IMGT Builder to generate a profile if you do not already "
                    "have one.")
        if not isinstance(sequences, list) or not isinstance(sequence_names, list):
            raise ValueError("Both sequences and sequence_names should be lists of strings.")
        if len(sequences) != len(sequence_names):
            raise ValueError("Sequences and sequence_names must have the same length.")

        for sequence in sequences:
            self._validate_sequence(sequence)


    def score_online_seqs(self, sequences, sequence_names, ncpu=0, skip_checks=False):
        """Numbers and scores a list of input sequences.

        Args:
            sequences (list): A list of strings, each of which is a sequence
                containing the usual 20 amino acids.
            sequence_names (list): A list of strings, each of which is the
                name of the corresponding sequence in sequences.
            ncpu (int): The number of threads hmmer is allowed to use. If 0, it will
                be autoselected.
            skip_checks (bool): If True, skip performing validity checks on the input.
                DANGER -- setting this to true can slightly reduce overhead but may
                cause unexpected, hard-to-diagnose errors. Think twice before setting
                to True!
        """

        if not skip_checks:
            self._in_memory_checks(sequences, sequence_names)

        loaded_seqs = [TextSequence(name=name.encode(), sequence=seq).digitize(self.alphabet)
                for (name, seq) in zip(sequence_names, sequences)]

        return [(h.name.decode(), h.score) for h in
                hmmscan(self.hmm_db, loaded_seqs, cpus=ncpu)]
