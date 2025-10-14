"""Contains tools for searching and updating a local database."""
import os
import numpy as np
from ..numbering_tools.cterm_finder import _load_nterm_kmers
from ..antpack_license import get_license_key_info
from ..utilities.vj_utilities import load_vj_gene_consensus_db
from antpack.antpack_cpp_ext import LocalDatabaseToolCpp



class LocalDBTool:
    """Contains tools for searching a local database,
    exporting or clustering its contents and adding
    or deleting records."""

    def __init__(self, database_path:str):
        """Class constructor.

        Args:
            database_path (str): The filepath for the
                database. All other necessary info about
                the database will be retrieved from the
                database.

        Raises:
            RuntimeError: A runtime error is raised if the
                database is not found, cannot be opened or
                does not contain expected tables and metadata.
        """
        license_key, user_email = get_license_key_info()
        project_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                "..")
        consensus_path = os.path.join(project_path,
                "numbering_tools", "consensus_data")
        nterm_kmer_dict = _load_nterm_kmers()

        vj_db_path = os.path.join(project_path, "vj_tools", "consensus_data")
        vj_names, vj_seqs, _ = load_vj_gene_consensus_db(os.getcwd(),
                vj_db_path, "imgt")

        blosum_matrix = np.load(os.path.join(project_path,
            "numbering_tools", "consensus_data", "mabs",
            "blosum_matrix.npy")).astype(np.float64)
        self.local_db_manager = LocalDatabaseToolCpp(database_path,
                license_key, user_email, consensus_path,
                nterm_kmer_dict, vj_names, vj_seqs, blosum_matrix)



    def search(self, seq:str, annotation:tuple,
            mode="3", cdr_cutoff=0.25,
            max_cdr_length_shift:int=2, max_hits:int=10,
            retrieve_closest_only:bool=True,
            vgene="", species=""):
        """Searches the database and returns a list of nearest
        neighbors that meet the input criteria. The search can
        be conducted using the full sequence or a sub-region
        of it (e.g. cdr3, the cdrs etc).

        Args:
            seq (str): A sequence. Can contain the normal 20 AAs
                and 'X' (although 'X' should be used sparingly).
            annotation (tuple): A tuple containing (numbering,
                percent_identity, chain_name, error_message). This
                is what you will get as output if you pass sequences
                to the analyze_seq method of SingleChainAnnotator or
                PairedChainAnnotator.
                
                Note that MUST be generated using the same numbering
                scheme and set of cdr definitions used to create the
                database, or you may get unexpected / nonsensical
                results. If you are unsure what set of definitions was
                used when the database was created, call
                get_database_metadata().
        """
        return self.local_db_manager.search(seq, annotation,
                mode, cdr_cutoff, max_cdr_length_shift,
                max_hits, retrieve_closest_only, vgene, species)

    def _retrieve_hit_dict(self, seq:str, annotation:tuple,
            cdr_cutoffs = [-1, -1, 0.25], max_cdr_length_shift = 2):
        """Used only for testing.
        """
        return self.local_db_manager.dummy_candidate_search(seq, annotation,
                cdr_cutoffs, max_cdr_length_shift,
                1000, "")


    def search_by_metadata(self, metadata:str, max_hits=100):
        """Searches the database and returns a list of sequences for
        which the metadata contains a string supplied in the input.
        Notice that this search is slow because it requires scanning
        the full database (no shortcuts).

        Args:
            metadata (str): A string that you want to search for.
                Sequences whose metadata contains this string will
                be returned. This is case-sensitive.
            max_hits (int): The maximum number of sequences to retrieve.
                Any hits above and beyond this will not be returned.

        Returns:
            hit_seqs (list): A list of sequences of length max_hits or
                less. An empty list if no hits were found. If more hits
                were found than max_hits, the *closest* hits are returned.
            hit_metadata (list): The metadata associated with each hit
                sequence (this is the sequence description from the
                fasta file used to create the database).
        """
        return self.local_db_manager.metadata_search(metadata, max_hits)




    def add_seq(self, seq:str, metadata:str):
        """Adds a sequence to the database. The sequence will
        be numbered and its cdrs assigned using the same schemes
        used to originally construct the db.

        Args:
            seq (str): A sequence which should be either a cdr3,
                all three concatenated cdrs, or a full heavy or
                light chain.
            metadata (str): Any metadata associated with the sequence
                (e.g. a description).
        """
        self.local_db_manager.add_sequence(seq, metadata)



    def get_database_metadata(self):
        """Returns metadata about the database (numbering scheme,
        receptor type, cdr definition type, and memo entered when
        creating the db if any).

        Returns:
            numbering_scheme (str): The numbering scheme (e.g. aho, imgt etc.)
                used for this database.
            receptor_type (str): The receptor type (mab or tcr).
            sequence_type (str): Either 'single', 'paired' or 'unknown'.
                If 'unknown', it is assumed that each sequence added may have
                either one chain or two and they are analyzed accordingly.
            cdr_scheme (str): The cdr definitions used for this database
                (e.g. aho, imgt etc.)
            user_memo (str): A memo supplied by the user when creating the
                database with some info about what it is for; may be "" if
                no memo was supplied.
        """
        return self.local_db_manager.get_database_metadata()


    def get_database_counts(self):
        """Returns key cdr positions for heavy and light chain
        mapped to the number of occurrences of each amino acid
        at each position. Mainly used for testing, not really
        intended for use by end users (but maybe occasionally
        useful).

        Returns:
            heavy_position_table (dict): A map from cdr positions to
                counts.
            light_position_table (dict): A map from cdr positions to
                counts.
        """
        return self.local_db_manager.get_position_tables()



    def export_database_to_csv(self, csv_filepath, gzip_file=False):
        """Takes all sequences from the database and exports
        them to a csv file with the following columns:
        Metadata:        The metadata for that sequence
        Sequence:        The raw sequence
        Heavy numbering: The numbering tokens concatenated with
                         '_' as a divider
        Light numbering: The numbering tokens concatenated with
                         '_' as a divider
        Heavy chain:     The heavy chain (if any)
        Light chain:     The light chain (if any)
        Vgene family:    The vgene family
        Jgene family:    The jgene family
        Vgenes:          All matching vgenes (by percent identity)
        Jgenes:          All matching jgenes (by percent identity)

        Args:
            csv_filepath (str): The filepath for the csv output.
            gzip_file (bool): If True, the file is gzipped and .gz is
                added to the specified filepath.

        Raises:
            RuntimeError: A RuntimeError is raised if the csv cannot be
                opened or an error is encountered.
        """
        pass
