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
                'multi' uses multithreading and works well on
                solid-state drives; it should not be used on hard
                drives where it may be slightly slower than single
                threading due to seek cost. 'single' uses
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
            mode:str="3", cdr_cutoff:float=0.25,
            blosum_cutoff:float=-1, max_cdr_length_shift:int=2,
            max_hits:int=1000, vgene:str="", species:str=""):
        """Searches the database and returns a list of nearest
        neighbors that meet the input criteria. The search can
        be conducted using the full sequence or a sub-region
        of it (e.g. cdr3, the cdrs etc).

        Args:
            seq (str): The input amino acid sequence. May contain
                'X', must not contain '*'.
            annotation (tuple): A tuple of (numbering, percent id,
                chain type, error message) that is returned by
                any of the analyze_seq functions for annotators
                in AntPack.
            mode (str): One of "3" or "123", indicating whether to
                look for sequences that are similar based on cdr3 or
                based on all three cdrs.
            cdr_cutoff (float): The maximum number of differences
                allowed from the query seq as a percentage of the
                query seq cdr length(s) (excluding gaps). This
                is assessed separately for cdr3 and cdrs 1/2 if
                using mode "123". In other words, if searching
                all three cdrs with cutoff 0.25, the maximum
                number of differences allowed in cdr3 is 25%
                of that cdr length rounded to the nearest integer,
                and the maximum number of differences in cdrs 1 & 2
                combined is 25% of their combined lengths rounded
                to the nearest integer. cdr_cutoff is allowed to
                range from 0 to 0.4 (if you want to use a larger
                value you will need to do a slower full-table scan
                instead of using this function).
            blosum_cutoff (float): Either -1 or a value between 0 and
                1 (inclusive). If < 0, this argument is ignored. If
                between 0 and 1, any match that has BLOSUM mismatch score
                less than (cdrlength * cutoff) for the query is
                discarded. Specifying a value between 0 and 1 makes
                the search more restrictive by excluding results that
                are within CDR cutoff but have a highly improbable
                mutation (as determined by BLOSUM scoring).
            max_cdr_length_shift (int): The maximum +/- amount by
                which the length of cdr3 for a match can differ from
                the query. If 0, the results are required to have cdr3
                be the same length as the query.
            max_hits (int): The maximum number of hits to return.
            vgene (str): One of "" or a valid vgene. If a valid vgene,
                hits are required to belong to the same v-gene family.
            species (str): One of "" or a valid species (human, mouse,
                alpaca, rabbit). If "", hits are required to belong to
                the same assigned species.

        Returns:
            hits (list): A list of the sequence id numbers for the hits.
                These can be used to retrieve the sequences and their
                metadata using the get_sequence call.
            distances (list): A list of the hamming distances of each
                hit to the query, unless a blosum cutoff >= 0 was supplied,
                in which case BLOSUM mismatch distance is calculated.
        """
        return self.local_db_manager.search(seq, annotation,
                mode, cdr_cutoff, blosum_cutoff,
                max_cdr_length_shift, max_hits,
                vgene, species)


    def search_batch(self, seqs:list, annotations:list,
            mode:str="3", cdr_cutoff:float=0.25,
            blosum_cutoff:float=-1,
            max_cdr_length_shift:int=2, max_hits:int=10,
            vgenes:list=[], species:list=[]):
        """Searches the database and returns a list of nearest
        neighbors that meet the input criteria for each of the input
        sequences. This is essentially a batch version of search.
        The search can be conducted using the full sequence or a
        sub-region of it (e.g. cdr3, the cdrs etc). This is faster
        than doing individual searches if multithreading is in use
        and you have a large number of queries.

        Args:
            seq (str): A list of input amino acid sequences. May contain
                'X', must not contain '*'.
            annotation (tuple): A list of tuples of (numbering, percent id,
                chain type, error message) that is returned by
                any of the analyze_seq functions for annotators
                in AntPack.
            mode (str): One of "3" or "123", indicating whether to
                look for sequences that are similar based on cdr3 or
                based on all three cdrs.
            cdr_cutoff (float): The maximum number of differences
                allowed from each query seq as a percentage of the
                query seq cdr length(s) (excluding gaps). This
                is assessed separately for cdr3 and cdrs 1/2 if
                using mode "123". In other words, if searching
                all three cdrs with cutoff 0.25, the maximum
                number of differences allowed in cdr3 is 25%
                of that cdr length rounded to the nearest integer,
                and the maximum number of differences in cdrs 1 & 2
                combined is 25% of their combined lengths rounded
                to the nearest integer. cdr_cutoff is allowed to
                range from 0 to 0.4 (if you want to use a larger
                value you will need to do a slower full-table scan
                instead of using this function).
            blosum_cutoff (float): Either -1 or a value between 0 and
                1 (inclusive). If < 0, this argument is ignored. If
                between 0 and 1, any match that has BLOSUM mismatch score
                less than (cdrlength * cutoff) for the query is
                discarded. Specifying a value between 0 and 1 makes
                the search more restrictive by excluding results that
                are within CDR cutoff but have a highly improbable
                mutation (as determined by BLOSUM scoring).
            max_cdr_length_shift (int): The maximum +/- amount by
                which the length of cdr3 for a match can differ from
                a query. If 0, the results are required to have cdr3
                be the same length as the query.
            max_hits (int): The maximum number of hits to return for
                each query.
            vgene (str): Either an empty list or a list of valid vgenes
                of the same length as seqs. If not empty, each hit is
                required to belong to the same vgene family as the vgene
                for that query.
            species (str): Either an empty list or a list of valid species
                of the same length as seqs. If not empty, each hit is
                required to belong to the same species as the species
                for that query.

        Returns:
            hits (list): A list of lists of sequence id numbers for hits
                for each query. List[i] is the hits for query [i].
                These can be used to retrieve the sequences and their
                metadata using the get_sequence call.
            distances (list): A list of lists of hamming distances of each
                hit to the query, unless a blosum cutoff >= 0 was supplied,
                in which case BLOSUM mismatch distance is calculated. List[i]
                is the distance for query[i].
        """
        return self.local_db_manager.search_seqs(seqs, annotations,
                mode, cdr_cutoff, blosum_cutoff,
                max_cdr_length_shift, max_hits,
                vgenes, species)


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



    def build_sparse_distance_matrix(self, chain_type:str="heavy",
            mode:str="3", cdr_cutoff:float=0.25,
            blosum_cutoff:float=-1,
            max_hits_per_query:int=10000,
            max_cdr_length_shift:int=0, distance_type:str="pid",
            filter_by_vgene:bool=True, filter_by_species:bool=True,
            verbose=True):
        """Assembles the lists needed to build a sparse
        distance matrix for the existing database. The matrix
        is sparse because distances for sequences outside the cutoff
        is not calculated. The resulting matrix can be used as input
        to sklearn's DBSCAN.

        It is important to realize the resulting matrix may not be
        symmetric depending on the settings that are selected. This will
        not however be a problem for the DBSCAN algorithm.

        Returns:
            chain_type (str): One of "heavy", "light".
            mode (str): One of "3" or "123", indicating whether to
                assess cutoffs and calculate distances based on cdr 3
                only or all 3 cdrs.
            cdr_cutoff (float): The maximum number of differences
                allowed from each sequence as a percentage of its
                cdr length(s) (excluding gaps). This
                is assessed separately for cdr3 and cdrs 1/2 if
                using mode "123". In other words, if searching
                all three cdrs with cutoff 0.25, the maximum
                number of differences allowed in cdr3 is 25%
                of that cdr length rounded to the nearest integer,
                and the maximum number of differences in cdrs 1 & 2
                combined is 25% of their combined lengths rounded
                to the nearest integer. cdr_cutoff is allowed to
                range from 0 to 0.4 (if you want to use a larger
                value you will need to do a slower full-table scan
                instead of using this function).
            blosum_cutoff (float): Either -1 or a value between 0 and
                1 (inclusive). If < 0, this argument is ignored. If
                between 0 and 1, any match that has BLOSUM mismatch score
                less than (cdrlength * cutoff) for the query is
                discarded. Specifying a value between 0 and 1 makes
                the search more restrictive by excluding results that
                are within CDR cutoff but have a highly improbable
                mutation (as determined by BLOSUM scoring).
            max_hits_per_query (int): Only up to this many neighbors
                will be stored per datapoint.
            max_cdr_length_shift (int): The maximum +/- amount by
                which the length of cdr3 for a match can differ from
                a query. If 0, the results are required to have cdr3
                be the same length as the query. Notice that the
                matrix may not be symmetric if this value is not 0,
                since it is possible for A to be inside B's cutoff
                but not vice versa.
            distance_type (str): One of 'pid', 'hamming', 'blosum'. If
                'pid', the number of differences from each sequence to each
                partner is calculated as a percentage of that sequence's
                cdr length(s) (depending on whether using cdr3 only or
                all three). If 'hamming' the hamming distance is
                calculated. If 'blosum' a blosum cutoff of 0 or greater
                must be supplied, and the blosum mismatch distance is
                calculated.
            filter_by_vgene (bool): If True, each sequence only considers
                possible partners that are in the same vgene family.
            filter_by_species (bool): If True, each sequence only considers
                possible partners that are assigned to the same species.
            verbose (bool): If True, print out updates every 10,000 sequences.

        Returns:
            distances (list): A list of length S for S distances in the dataset
                that are inside the cdr cutoff.
            row_idx (list): The row id (to use in the resulting sparse
                distance matrix) for each distance in distances.
            col_idx (list): The col idx (to use in the resulting sparse
                distance matrix) for each distance in distances. To build
                a sparse distance matrix using scipy's csr_matrix,
                call `csr_matrix((distances, (row_idx, col_idx)),
                    [shape=(num_seqs, num_seqs)])`.
        """
        return self.local_db_manager.build_sparse_dmat_lists(chain_type,
                mode, cdr_cutoff, blosum_cutoff, max_hits_per_query,
                max_cdr_length_shift, distance_type, filter_by_vgene,
                filter_by_species, verbose)



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


    def get_num_seqs(self, chain_type="all"):
        """Gets either the total number of sequences OR
        the total number of entries in the heavy chain table
        or the light chain table, as requested.

        Args:
            chain_type (str): One of "all", "heavy", "light".
                Determines whether to count all chains or just
                heavy or light.

        Returns:
            count (int): The number of sequences in the table
                matching the specifications.
        """
        return self.local_db_manager.get_sequence_count(chain_type)


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
