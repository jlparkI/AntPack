"""Contains tools for searching a local database."""
from ..antpack_license import get_license_key_info
from antpack.antpack_cpp_ext import LocalDatabaseToolCpp



class LocalDBSearchTool:
    """Contains tools for searching a local database and
    clustering its contents."""

    def __init__(self, database_path:str,
                 multithread:bool=False):
        """Class constructor.

        Args:
            database_path (str): The filepath for the
                database. All other necessary info about
                the database will be retrieved from the
                database.
            multithread (bool): If True, open database in such
                a way that multithreading / batch search is
                available, otherwise do not. We strongly recommend
                setting this to False on hard drives; on solid-state
                drives, by contrast, batch search is faster than
                single search if querying a large number of
                sequences so this may be advantageous.

        Raises:
            RuntimeError: A runtime error is raised if the
                database is not found, cannot be opened or
                does not contain expected tables and metadata.
        """
        license_key, user_email = get_license_key_info()
        self.local_db_manager = LocalDatabaseToolCpp(database_path,
                license_key, user_email, multithread)



    def search(self, seq:str, annotation:tuple,
            mode:str="3", cdr_cutoff:float=0.25,
            blosum_cutoff:float=-1, max_cdr_length_shift:int=2,
            use_family_only=True, symmetric_search=False,
            vgene:str="", species:str="", jgene:str=""):
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
                range from 0 to 0.33 (if you want to use a larger
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
            use_family_only (bool): If True, when filtering by
                vgene, hits are required to have the same vgene *family*,
                but can have a different vgene within that family. If
                False, hits are required to have the same vgene AND
                family. This argument is ignored if vgene or species
                is "".
            symmetric_search (bool): If True, the cdr cutoff is applied
                 both forward and reverse -- in other words, the Hamming
                 distance from hit to query must be less than either
                 cdr_cutoff * query_cdr_length OR cdr_cutoff *
                 hit_cdr_length. This ensures that distance matrices
                 constructed using search are symmetric. If not planning
                 to build a distance matrix, this is unnecessary.
            vgene (str): One of "" or a valid vgene. If a valid vgene,
                hits are required to match either the vgene family or
                the specific vgene (depending on the
                use_vjgene_family_only argument). Ignored if no species
                is supplied.
            species (str): One of "" or a valid species (human, mouse,
                alpaca, rabbit). If "", hits are required to belong to
                the same assigned species.
            jgene (str): One of "" or a valid jgene. If a valid jgene,
                hits are required to match either the jgene family or
                the specific jgene (depending on the
                use_vjgene_family_only argument). Ignored if no
                species or vgene is supplied.

        Returns:
            hits (list): A list of tuples containing (sequence_id,
                distance, num_non_canonical_positions). The sequence
                ids can be used to retrieve the sequences and their
                metadata using the get_sequence call. Non-canonical
                positions are rare unusual insertions that are not
                included in the distance calculation; this value will
                normally be zero. If it is not, it indicates the hit
                contains unusual non-canonical positions.
            error_message (str): An error message which is "" if
                all went well and is informative otherwise.
        """
        return self.local_db_manager.search(seq, annotation,
                mode, cdr_cutoff, blosum_cutoff,
                max_cdr_length_shift, use_family_only,
                symmetric_search, vgene, species, jgene)



    def search_from_preprocessed_data(self, prepped_inputs:tuple,
            mode:str="3", cdr_cutoff:float=0.25,
            blosum_cutoff:float=-1, max_cdr_length_shift:int=2,
            use_family_only=True, symmetric_search=False):
        """Searches the database using preprocessed query info
        generated by ```retrieve_preprocessed_search_data```
        and returns a list of nearest neighbors that meet the
        input criteria. This function is a good alternative
        to ```search``` if you are retrieving specific sequences
        from one database and searching them against another,
        in which case it is much faster than numbering / prepping
        the retrieved sequences yourself.

        Args:
            prepped_inputs (tuple): A tuple of information returned
                by ```retrieve_preprocessed_search_data```.
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
                range from 0 to 0.33 (if you want to use a larger
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
            use_family_only (bool): If True, when filtering by
                vgene, hits are required to have the same vgene *family*,
                but can have a different vgene within that family. If
                False, hits are required to have the same vgene AND
                family. This argument is ignored if vgene or species
                is "".
            symmetric_search (bool): If True, the cdr cutoff is applied
                 both forward and reverse -- in other words, the Hamming
                 distance from hit to query must be less than either
                 cdr_cutoff * query_cdr_length OR cdr_cutoff *
                 hit_cdr_length. This ensures that distance matrices
                 constructed using search are symmetric. If not planning
                 to build a distance matrix, this is unnecessary.

        Returns:
            hits (list): A list of tuples containing (sequence_id,
                distance, num_non_canonical_positions). The sequence
                ids can be used to retrieve the sequences and their
                metadata using the get_sequence call. Non-canonical
                positions are rare unusual insertions that are not
                included in the distance calculation; this value will
                normally be zero. If it is not, it indicates the hit
                contains unusual non-canonical positions.
            error_message (str): An error message which is "" if
                all went well and is informative otherwise.
        """
        return self.local_db_manager.search_from_preprocessed_data(
                prepped_inputs, mode, cdr_cutoff,
                blosum_cutoff, max_cdr_length_shift,
                symmetric_search)



    def retrieve_preprocessed_search_data(self,
                    seq_id:int, chain_type:str="heavy",
                    filter_by_vgene:bool=True,
                    filter_by_jgene:bool=True,
                    use_vgene_family_only:bool=False):
        """Retrieves preprocessed search data that can be supplied
        directly to ```search_from_preprocessed_data```. This is
        useful if you are retrieving sequences from one database
        to search against another database created with the same
        numbering and cdr scheme.

        Args:
            seq_id (int): A sequence id associated with one of the
                sequences in the database.
            chain_type (str): One of 'heavy', 'light'.
            filter_by_vgene (bool): If True, the returned results
                are set up so that the resulting search will find
                only sequences with the same vgene and vgene family.
                If False, the returned results are configured to
                retrieve sequences regardless of vgene.
            filter_by_jgene (bool): If True, the returned results
                are set up so that the resulting search will find
                only sequences with the same jgene and jgene family.
                If False, the returned results are configured to
                retrieve sequences regardless of jgene.
            use_vgene_family_only (bool): If True, the returned
                results are set up so that the resulting search
                will filter for sequences with the same vgene family
                (but possibly a different vgene). Ignored if
                filter_by_vgene is False.

        Returns:
            preprocessed_inputs (tuple): A tuple of query_seq
                (the three concatenated cdrs), cdr_lengths (the
                length of all three cdrs), chain code, vgene_info
                and jgene_info. This can be supplied directly to
                ```search_from_preprocessed_data```.
        """
        return self.local_db_manager.retrieve_search_data(
                seq_id, chain_type, filter_by_vgene,
                filter_by_jgene, use_vgene_family_only)





    def get_sequence(self, seq_id:int):
        """Retrieves the sequence and metadata associated with
        a given sequence id. Sequence ids can be obtained from
        any of the search functions.

        Args:
            seq_id (int): A sequence id associated with one of the
                sequences in the database.

        Returns:
            sequence (str): Either a blank string if the id that was
                supplied does not match any database sequences, or the
                sequence from the database.
            metadata (str): Either a blank string if the id that
                was supplied does not match any database sequences,
                or the metadata associated with the sequence from the
                database.
        """
        return self.local_db_manager.retrieve_sequence(seq_id)


    def get_vgene_jgene(self, seq_id:int, chain_type:str):
        """Retrieves the vgene and jgene associated with
        a given sequence id and chain type. Sequence ids can be
        obtained from any of the search functions.

        Args:
            seq_id (int): A sequence id associated with one of the
                sequences in the database.
            chain_type (str): One of 'heavy', 'light'. If using TCRs,
                B is heavy.

        Returns:
            vgene (str): Either a blank string if the id that was supplied
                does not match any database sequences or there is no vgene
                for this sequence, or the vgene from the database.
            metadata (str): Either a blank string if the id that was supplied
                does not match any database sequences or there is no jgene
                for this sequence, or the jgene from the database.
        """
        return self.local_db_manager.retrieve_vgene_jgene(seq_id, chain_type)



    def basic_clustering(self, chain_type:str="heavy",
        mode:str="123", cdr_cutoff:float=0.2, blosum_cutoff:float=-1,
        verbose:bool=True):
        """Clusters the linked database using single linkage, by searching
        every sequence once against the database. Only sequences that have the
        same cdr3 length and the same vgenes, jgenes and species can
        be assigned to a cluster.

        Args:
            chain_type (str): One of "heavy", "light". For TCRs, TCRB is
                heavy.
            mode (str): One of "123" or "3", indicating whether to use
                cdr3 only or cdrs 1, 2 and 3 when determining whether a
                hit meets neighborhood criteria. Generally use 123 UNLESS
                you are working with sequences where you have data for cdr3
                only (e.g. TCRs).
            cdr_cutoff (float): A value in the range from 0 to 0.25
                (inclusive). Determines the percent identity for two
                sequences to be in the same cluster.
            blosum_cutoff (float): Either -1 or a value >=0. If -1, this
                argument is ignored. If >=0, sequences must have a BLOSUM
                distance < this cutoff at all positions to be in the same
                cluster.
            verbose (bool): If True, print periodic updates while clustering.

        Returns:
            cluster_assignments (list): A list of ints of the same length as the
                number of sequences in the database, including all chain
                types. Each sequence has a cluster label which may be -2 for
                sequence does not have chain of this type, -1 for sequence is
                not assigned to any cluster / has no neighbors within the
                search radius, or a number >= 0 indicating the cluster to which
                it belongs.
        """
        return self.local_db_manager.basic_clustering(
                chain_type, mode, cdr_cutoff, blosum_cutoff,
                verbose)


    def get_database_metadata(self):
        """Returns metadata about the database (numbering scheme,
        receptor type, cdr definition type, and memo entered when
        creating the db if any).

        Returns:
            numbering_scheme (str): The numbering scheme (e.g. aho,
                imgt etc.) used for this database.
            receptor_type (str): The receptor type (mab or tcr).
            sequence_type (str): Either 'single', 'paired' or 'unknown'.
                If 'unknown', it is assumed that each sequence added
                may have either one chain or two and they are analyzed
                accordingly.
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


    def get_max_seqid(self):
        """Gets the largest seqid in the database. This will be
        the same as the number of stored sequences IF no database
        rows have been deleted.

        Returns:
            seq_id (int): The largest sequence id in the database.
        """
        return self.local_db_manager.get_max_seqid()


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
