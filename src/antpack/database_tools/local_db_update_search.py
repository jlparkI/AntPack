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

    def __init__(self, database_path:str,
            thread_mode:str="single"):
        """Class constructor.

        Args:
            database_path (str): The filepath for the
                database. All other necessary info about
                the database will be retrieved from the
                database.
            thread_mode (str): One of 'single', 'multi',
                'prodcon'. 'multi' uses multithreading
                and works well on solid-state drives; it should not
                be used on hard drives where it may be slightly
                slower than single threading due to seek cost.
                'prodcon' uses two threads at all times
                and works better for hard drives. 'single' uses
                only a single thread.

        Raises:
            RuntimeError: A runtime error is raised if the
                database is not found, cannot be opened or
                does not contain expected tables and metadata.
        """
        license_key, user_email = get_license_key_info()
        self.local_db_manager = LocalDatabaseToolCpp(database_path,
                license_key, user_email, thread_mode)



    def search(self, seq:str, annotation:tuple,
            mode="3", cdr_cutoff=0.25,
            max_cdr_length_shift:int=2, max_hits:int=10,
            vgene="", species=""):
        """Searches the database and returns a list of nearest
        neighbors that meet the input criteria. The search can
        be conducted using the full sequence or a sub-region
        of it (e.g. cdr3, the cdrs etc).

        Args:
        """
        return self.local_db_manager.search(seq, annotation,
                mode, cdr_cutoff, max_cdr_length_shift,
                max_hits, vgene, species)


    def search_batch(self, seqs:list, annotations:list,
            mode="3", cdr_cutoff=0.25,
            max_cdr_length_shift:int=2, max_hits:int=10,
            vgenes:list=[], species:list=[]):
        """Searches the database and returns a list of nearest
        neighbors that meet the input criteria for each of the input
        sequences. This is essentially a batch version of search.
        The search can be conducted using the full sequence or a
        sub-region of it (e.g. cdr3, the cdrs etc).

        Args:
        """
        return self.local_db_manager.search_seqs(seqs, annotations,
                mode, cdr_cutoff, max_cdr_length_shift,
                max_hits, vgenes, species)


    def get_sequence(self, seq_id:int):
        """Retrieves the sequence and metadata associated with
        a given sequence id. Sequence ids can be obtained from
        any of the search functions.

        Args:
            seq_id (int): A sequence id associated with one of the
                sequences in the database.

        Returns:
            sequence (str): Either a blank string if the id that was supplied
                does not match any database sequences, or the sequence from
                the database.
            metadata (str): Either a blank string if the id that was supplied
                does not match any database sequences, or the metadata associated
                with the sequence from the database.
        """
        return self.local_db_manager.retrieve_sequence(seq_id)




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
        mapped to the number of occurrences of each dimer
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
