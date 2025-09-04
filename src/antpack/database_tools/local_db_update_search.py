"""Contains tools for searching and updating a local database."""
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
        self.local_db_manager = LocalDatabaseToolCpp(database_path)



    def search_seq(self, seq:str, annotation:tuple,
            max_hamming_dist:int=5, max_hits:int=100,
            region_label:str="cdr3", require_same_length_cdrs:bool=False):
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
            max_hamming_dist (int): The maximum number of differences
                between the input sequence and any hits in the target
                region (inclusive). If searching cdr3, you can use
                any number up to 5. If searching all three cdrs or the
                whole sequence, you can use any number up to 10.
            max_hits (int): The maximum number of hits to retrieve.
                The closest number of nearest neighbors up to this
                value will be reported.
            region_label (str): One of 'cdr3', 'cdr3' or 'all'. If
                'cdr3', only cdr3 is searched. If 'cdr', all 3 cdrs
                are used. 'all' searches for hits using the full
                sequence.
            require_same_length_cdrs (bool): If True, only sequences whose
                cdr3 (if querying by cdr3) or all cdrs (for other queries)
                have the same length as query are returned.

        Returns:
            hit_seqs (list): A list of sequences of length max_hits or
                less. An empty list if no hits were found. If more hits
                were found than max_hits, the *closest* hits are returned.
            hit_metadata (list): The metadata associated with each hit
                sequence (this is the sequence description from the
                fasta file used to create the database).
        """
        return self.local_db_manager.search(seq, annotation,
                max_hamming_dist, max_hits, region_label,
                require_same_length_cdrs)



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
            
        """
        return self.local_db_manager.get_database_metadata()



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
