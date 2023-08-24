"""Contains the AntibodyAnnotator class, which provides the tools needed
to parse antibody sequences (and makes use of the functions supplied in
numbering_toolkit)."""
import pyhmmer
from pyhmmer.hmmer import hmmscan
from pyhmmer.easel import TextSequence

from .constants import allowed_inputs








class AntibodyAnnotator:
    """This class contains the tools needed to parse and number
    either a single sequence, a list of sequences, or a fasta file containing
    many sequences. If the latter, the class provides a generator so that
    you can analyze one sequence at a time (without loading all into memory
    simultaneously).
    """

    def __init__(self):
        """Class constructor.

        Args:
        """

        self.hmm_db = None
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


    def _in_memory_checks(self, sequences, sequence_names, scheme):
        """Checks user inputs for in-memory sequence analysis.

        Args:
            sequences (list): A list of sequence strings.
            sequence_names (list): A list of names for each sequence.
            scheme (str): The numbering scheme. Should be one of 'imgt', 'martin',
                'chothia', 'kabat', 'aho', 'wolfguy'.

        Raises:
            ValueError: A ValueError is raised if inappropriate inputs are detected.
        """
        if self.hmm_db is None:
            raise ValueError("An hmm profile has not been added! Consider using "
                    "IMGT Builder to generate a profile if you do not already "
                    "have one.")
        if scheme not in allowed_inputs.allowed_schemes:
            raise ValueError(f"Input scheme {scheme} is not allowed; "
                "should be one of {allowed_inputs.allowed_schemes}")
        if not isinstance(sequences, list) or not isinstance(sequence_names, list):
            raise ValueError("Both sequences and sequence_names should be lists of strings.")
        if len(sequences) != len(sequence_names):
            raise ValueError("Sequences and sequence_names must have the same length.")

        for sequence in sequences:
            self._validate_sequence(sequence)


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




    ## Parsing and recognising domain hits from hmmscan ##
    def _check_for_overlap(self, domain1, domain2):
        """Ensure that two domains do not overlap / are not redundant."""
        domain1, domain2 = sorted( [domain1, domain2], key=lambda x: x.target_from  )
        return domain2.target_from < domain1.target_to


    def _parse_hmmer_hits(self, all_hits, scheme, bit_score_threshold=80, keep_all_scores=False,
                assign_germline=False):
        """Parses the hits returned by pyhmmer.

        Args:
            all_hits (list): A list of TopHits objects. Each TopHits is the set of hits
                for the corresponding sequence from the input list.
            scheme (str): The numbering scheme. Should be one of 'imgt', 'martin',
                'chothia', 'kabat', 'aho', 'wolfguy'.
            bit_score_threshold (int): Bit scores below this threshold are
                not considered hits and are removed.
            keep_all_scores (bool): If True, store scores for all profile hits
            assign_germline (bool): If True, indicate the germline genes that most
                closely correspond to the search sequences.
        """
        all_alignments = []
        all_scores = []

        for hit_list in all_hits:
            clean_hit_list = [hit for hit in hit_list if hit.score >
                    bit_score_threshold]
            if keep_all_scores:
                all_scores.append([(h.name.decode(), h.score) for h in clean_hit_list])
            if len(clean_hit_list) == 0:
                all_alignments.append([])
                continue
            best_hit  = max( (v, i) for i, v in enumerate(clean_hit_list) )[0]

            retained_domains = []
            for domain in best_hit.domains:
                is_overlap = [self._check_for_overlap(previous_domain, domain) for previous_domain
                        in retained_domains]
                if sum(is_overlap) == 0:
                    retained_domains.append(domain)
            retained_domains = sorted(retained_domains, key=lambda x: domain.target_from)

        return all_alignments, all_scores
