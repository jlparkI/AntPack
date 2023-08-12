"""Provides tools for scoring an antibody sequence against germlines of
different species to determine which is the most likely species of origin."""
import os
from Bio import SeqIO
import pyhmmer
from pyhmmer.hmmer import hmmscan
from pyhmmer.easel import TextSequence

from .constants import allowed_inputs





class AntibodyScoring:
    """This class contains the tools needed to score an antibody against
    profiles representing different species (and thereby determine the
    most likely species of origin).

    Attributes:
        all_chains_db: A pyhmmer optimized profile for all chains.
        heavy_ig_db: A pyhmmer optimized profile for heavy chains only.
        light_ig_db: A pyhmmer optimized profile for light chains only.
        alphabet: A pyhmmer alphabet (generally standard AAs).
    """

    def __init__(self, all_chains_path = None, ig_heavy_only_path = None,
            ig_light_only_path = None):
        """Class constructor.

        Args:
            all_chains_path (str): The filepath for an HMM db containing all chains,
                including T cell receptors, for all species. If None, the default
                db built into the package is used; the default db was constructed on 8/11/23.
            ig_heavy_only_path (str): The filepath for an HMM db containing IG heavy chains
                for all species. If None, the default db built into the package is used;
                the default db was constructed on 8/11/23.
            ig_light_only_path (str): The filepath for an HMM db containing IG light chains
                species. If None, the default db built into the package is used;
                the default db was constructed on 8/11/23.
        """
        f_dir = os.path.abspath(os.path.dirname(__file__))
        if all_chains_path is None:
            all_chains_path = os.path.join(f_dir, "hmm_dbs", "all", "ALL.hmm")
        if ig_heavy_only_path is None:
            ig_heavy_only_path = os.path.join(f_dir, "hmm_dbs", "heavy", "HEAVY.hmm")
        if ig_light_only_path is None:
            ig_light_only_path = os.path.join(f_dir, "hmm_dbs", "light", "LIGHT.hmm")

        with pyhmmer.plan7.HMMFile(all_chains_path) as hmm_file:
            self.all_chains_db = hmm_file.optimized_profiles()
        with pyhmmer.plan7.HMMFile(ig_heavy_only_path) as hmm_file:
            self.heavy_ig_db = hmm_file.optimized_profiles()
        with pyhmmer.plan7.HMMFile(ig_light_only_path) as hmm_file:
            self.light_ig_db = hmm_file.optimized_profiles()
        self.alphabet = pyhmmer.easel.Alphabet.amino()


    def _check_chain_type(self, chain_types):
        """Checks the input chain type and returns a reference to the appropriate
        database."""
        if chain_types == "ALL":
            db_file = self.all_chains_db
        elif chain_types == "HEAVY":
            db_file = self.heavy_ig_db
        elif chain_types == "LIGHT":
            db_file = self.light_ig_db
        else:
            raise ValueError("Unrecognized chain_types supplied.")
        return db_file


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
        if self.all_chains_db is None or self.heavy_ig_db is None or self.light_ig_db is None:
            raise ValueError("An hmm profile has not been added! Consider using "
                    "IMGT Builder to generate a profile if you do not already "
                    "have one.")
        if not isinstance(sequences, list) or not isinstance(sequence_names, list):
            raise ValueError("Both sequences and sequence_names should be lists of strings.")
        if len(sequences) != len(sequence_names):
            raise ValueError("Sequences and sequence_names must have the same length.")

        for sequence in sequences:
            self._validate_sequence(sequence)


    def score_online_seqs(self, sequences, sequence_names, chain_types = "ALL",
            ncpu=0, skip_checks=False):
        """Scores a list of input sequences.

        Args:
            sequences (list): A list of strings, each of which is a sequence
                containing the usual 20 amino acids.
            sequence_names (list): A list of strings, each of which is the
                name of the corresponding sequence in sequences.
            chain_types (str): One of 'ALL', 'HEAVY', 'LIGHT'. Determines which db
                is used. ALL can be used to determine chain type as well if in doubt,
                while HEAVY or LIGHT will be faster if you know for example that this
                is a heavy chain.
            ncpu (int): The number of threads hmmer is allowed to use. If 0, it will
                be autoselected.
            skip_checks (bool): If True, skip performing validity checks on the input.
                DANGER -- setting this to true can slightly reduce overhead but may
                cause unexpected, hard-to-diagnose errors. Think twice before setting
                to True!

        Returns:
            score_results (list): A list of lists of tuples. Each sublist is the
                list of (species_chain, bitscore) for the corresponding input sequence.

        Raises:
            ValueError: A ValueError is raised if unexpected arguments are passed.
        """
        if not skip_checks:
            self._in_memory_checks(sequences, sequence_names)

        db_file = self._check_chain_type(chain_types)

        loaded_seqs = [TextSequence(name=name.encode(), sequence=seq).digitize(self.alphabet)
                for (name, seq) in zip(sequence_names, sequences)]

        return [[(h.name.decode(), h.score) for h in tophit] for tophit in
                hmmscan(loaded_seqs, db_file, cpus=ncpu)]


    def score_offline_seqs(self, fasta_file, chain_types = "ALL", ncpu=0):
        """Scores all the sequences in a fasta file.

        Args:
            fasta_file (str): A valid filepath to a fasta file.
            chain_types (str): One of 'ALL', 'HEAVY', 'LIGHT'. Determines which db
                is used. ALL can be used to determine chain type as well if in doubt,
                while HEAVY or LIGHT will be faster if you know for example that this
                is a heavy chain.
            ncpu (int): The number of threads hmmer is allowed to use. If 0, it will
                be autoselected.

        Returns:
            score_results (list): A list of lists of tuples. Each sublist is the
                list of (species_chain, bitscore) for the corresponding input sequence.

        Raises:
            ValueError: A ValueError is raised if unexpected arguments are passed.
        """
        try:
            with open(fasta_file, "r", encoding="utf-8") as fhandle:
                pass
        except Exception as exc:
            raise ValueError("Invalid filepath supplied.") from exc

        db_file = self._check_chain_type(chain_types)

        scores = []
        with open(fasta_file, "r", encoding="utf-8") as fhandle:
            for record in SeqIO.parse(fhandle, 'fasta'):
                loaded_seq = TextSequence(name=record.name.encode(),
                        sequence=str(record.seq)).digitize(self.alphabet)
                hits = list(hmmscan([loaded_seq], db_file, cpus=ncpu))
                scores.append( [(hit.name.decode(), hit.score) for hit in hits[0]] )

        return scores
