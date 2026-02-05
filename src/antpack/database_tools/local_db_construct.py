"""Wraps tools for constructing a searchable database from
various input types."""
import os
import gc
from ..utilities import read_fasta, read_csv
from ..numbering_tools.cterm_finder import _load_nterm_kmers
from ..antpack_license import get_license_key_info
from ..utilities.vj_utilities import load_vj_gene_consensus_db
from antpack.antpack_cpp_ext import DatabaseConstructionTool



def build_database_from_fasta(fasta_filepaths:list,
        database_filepath:str, numbering_scheme:str="imgt",
        cdr_definition_scheme:str="imgt",
        sequence_type:str="single", receptor_type:str="mab",
        pid_threshold:float=0.7, user_memo:str="",
        reject_file:str = None,
        verbose:bool=True):
    """Builds a database from a list of fasta files which may or may
    not be gzipped. The database is constructed so it can be
    searched quickly and the sequence descriptions for each sequence
    in the fasta file are saved as metadata.

    Args:
        fasta_filepath (str): A list of fasta filepath locations, which
            may or may not be gzipped.
        database_filepath (str): The desired location and filename
            for the database.
        numbering_scheme (str): One of 'imgt', 'kabat', 'martin' or
            'aho'. If receptor_type is 'tcr', 'imgt' is the only
            allowed option.
        cdr_definition_scheme (str): One of 'imgt', 'kabat', 'martin',
            'aho' or 'north'. If receptor_type is 'tcr', 'imgt' only is allowed.
        sequence_type (str): One of 'single', 'paired', 'unknown'. If
            'paired' each sequence is assumed to be paired. If 'unknown'
            it is assumed each sequence MAY be paired and should be analyzed
            as paired just in case.
        receptor_type (str): One of 'mab', 'tcr'.
        pid_threshold (float): A value between 0 and 1 for percent identity
            threshold. If sequence_type is 'single' or 'unknown', sequences
            not meeting this threshold are rejected. If sequence_type is
            'paired' the sequences are still retained as long as one of
            the chains meets this threshold.
        user_memo (str): A string describing the purpose of the database / anything
            important you want your future self or other users to know about the
            contents. Will be saved as part of the database metadata.
        reject_file (str): Either None or a filepath. If None, any sequences that
            are rejected are silently ignored. If a filepath, rejected sequences
            are written to that filepath which is saved as a fasta file.
        verbose (bool): If True, print regular updates while running.

    Raises:
        RuntimeError: A RuntimeError is raised if invalid arguments are
            supplied.
    """
    if not isinstance(fasta_filepaths, list):
        raise RuntimeError("The fasta filepaths should be a list of filepaths. "
                "If you have only one filepath, you can enclose it in brackets "
                "to make a list.")
    if os.path.exists(database_filepath):
        raise RuntimeError("The database already exists.")

    license_key, user_email = get_license_key_info()
    project_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
            "..")
    consensus_path = os.path.join(project_path,
            "numbering_tools", "consensus_data")
    nterm_kmer_dict = _load_nterm_kmers()

    vj_db_path = os.path.join(project_path, "vj_tools", "consensus_data")
    vj_names, vj_seqs, _ = load_vj_gene_consensus_db(os.getcwd(),
            vj_db_path, "imgt")

    print("Estimating number of sequences.")
    nseqs = 0

    for fasta_filepath in fasta_filepaths:
        for (seqinfo, seq) in read_fasta(fasta_filepath):
            nseqs += 1

    if nseqs == 0:
        print("No sequences found.")
        return

    db_construct_tool = DatabaseConstructionTool(database_filepath,
            numbering_scheme, cdr_definition_scheme,
            sequence_type, receptor_type,
            pid_threshold, license_key, user_email,
            consensus_path, nterm_kmer_dict,
            vj_names, vj_seqs, user_memo)

    if nseqs > 4200000000:
        raise RuntimeError("Current cap on number of sequences "
                           "per database file is 4.2 billion. "
                           f"You have {nseqs}. Consider splitting "
                           "this up into multiple database files, "
                           "e.g. splitting by species or vgene family.")

    print(f"Found {nseqs} sequences. Starting db construction.")

    seqcount = 0
    if reject_file is not None:
        reject_handle = open(reject_file, "w+", encoding="utf-8")
    else:
        reject_handle = None

    db_construct_tool.open_transaction()

    for fasta_filepath in fasta_filepaths:
        for (seqinfo, seq) in read_fasta(fasta_filepath):
            rcode = db_construct_tool.add_sequence(seq, seqinfo)
            if len(rcode) > 0:
                if reject_handle is not None:
                    reject_handle.write(f">{seqinfo}\t{rcode}\n{seq}\n")
            seqcount += 1
            if seqcount % 10000 == 0:
                db_construct_tool.close_transaction()
                db_construct_tool.open_transaction()
                if verbose:
                    print(f"{seqcount} complete.")

    if reject_handle is not None:
        reject_handle.close()

    if verbose:
        print("Now constructing database indices...")

    db_construct_tool.close_transaction()
    db_construct_tool.open_transaction()

    db_construct_tool.finalize_db_construction(verbose)
    gc.collect()



def build_database_from_full_chain_csv(csv_filepaths:list,
        database_filepath:str, column_selections:dict,
        header_rows:int=1,
        numbering_scheme:str="imgt",
        cdr_definition_scheme:str="imgt",
        delimiter:str=',',
        receptor_type:str="mab", pid_threshold:float=0.7,
        user_memo:str="", reject_file:str = None,
        verbose:bool=True):
    """Builds a database from a list of csv files which may or may
    not be gzipped. The database is constructed so it can be
    searched quickly. The csv files should already contain heavy
    and light chains split up into appropriate columns. The csv files
    can optionally contain additional columns with numbering or Vgenes
    -- if present these will save time for database construction since
    AntPack will not need to perform numbering or vgene identification.
    It is ok for some rows to contain both a heavy and light chain and
    others to contain only one or the other with "" for the missing chain.
    You can also load files that contain a heavy or light chain only.

    Use this option if you have the full chains. If you only have pre-
    extracted cdrs, use build_database_from_cdr_only_csv.

    Note that the beta chain for TCRs is assumed to be "heavy" and alpha
    is light.

    Args:
        csv_filepaths (list): A list of csv filepaths, which may or
            may not be gzipped.
        database_filepath (str): The desired location and filename
            for the database.
        column_selections (dict): A dictionary which should contain at least
            some of the following keys:

            * ``"heavy_chain"``: The number (from 0) of the column in each
              csv file containing the heavy chains. Do not include this
              if there are no heavy chains. If heavy chains are supplied,
              heavy vgenes and jgenes must be supplied also.
            * ``"light_chain"``: The number (from 0) of the column in each
              csv file containing the light chains. Do not include this
              if there are no light chains. If light chains are supplied,
              light vgenes and jgenes must be supplied also.
            * ``"heavy_numbering"``: The number (from 0) of the column (if any)
              containing numbering for the heavy chains. The numbering
              should be delimited using the '_' character. If you number
              your sequences using AntPack you can obtain this using
              '_'.join(antpack_numbering) for each sequence. If this is
              not supplied each heavy chain will be numbered by this
              function.
            * ``"light_numbering"``: The number (from 0) of the column (if any)
              containing numbering for the light chains. The numbering
              should be delimited using the '_' character. If you number
              your sequences using AntPack you can obtain this using
              '_'.join(antpack_numbering) for each sequence. If this is
              not supplied each light chain will be numbered by this
              function.
            * ``"heavy_vgene"``: The number (from 0) of the column (if any)
              containing vgene assignments for the heavy chains.
            * ``"light_vgene"``: The number (from 0) of the column (if any)
              containing vgene assignments for the light chains.
            * ``"heavy_jgene"``: The number (from 0) of the column (if any)
              containing jgene assignments for the heavy chains.
            * ``"light_jgene"``: The number (from 0) of the column (if any)
              containing jgene assignments for the light chains.
            * ``"species"``: The number (from 0) of the column (if any)
              containing species assignments. If this is not supplied, it
              will be assumed that all are human.
            * ``"metadata"``: The number (from 0) of the column (if any)
              containing metadata for each heavy / light chain or chain pair.

        header_rows (int): The number of header rows. Header rows are
            skipped.
        numbering_scheme (str): One of 'imgt', 'kabat', 'martin' or
            'aho'. If receptor_type is 'tcr', 'imgt' is the only
            allowed option.
        cdr_definition_scheme (str): One of 'imgt', 'kabat', 'martin',
            'aho' or 'north'. If receptor_type is 'tcr', 'imgt' only is allowed.
        delimiter (str): The delimiter for the csv files.
        receptor_type (str): One of 'mab', 'tcr'.
        pid_threshold (float): A value between 0 and 1 for percent identity
            threshold. If you have not supplied pre-written numbering,
            sequences are numbered during construction, and any not meeting
            this threshold are rejected.
        user_memo (str): A string describing the purpose of the database / anything
            important you want your future self or other users to know about the
            contents. Will be saved as part of the database metadata.
        reject_file (str): Either None or a filepath. If None, any sequences that
            are rejected are silently ignored. If a filepath, rejected sequences
            are written to that filepath which is saved as a csv file.
        verbose (bool): If True, print regular updates while running.

    Raises:
        RuntimeError: A RuntimeError is raised if invalid arguments are
            supplied.
    """
    if not isinstance(csv_filepaths, list):
        raise RuntimeError("The csv filepaths should be a list of filepaths. "
                "If you have only one filepath, you can enclose it in brackets "
                "to make a list.")
    if os.path.exists(database_filepath):
        raise RuntimeError("The database already exists.")

    license_key, user_email = get_license_key_info()
    project_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
            "..")
    consensus_path = os.path.join(project_path,
            "numbering_tools", "consensus_data")
    nterm_kmer_dict = _load_nterm_kmers()

    vj_db_path = os.path.join(project_path, "vj_tools", "consensus_data")
    vj_names, vj_seqs, _ = load_vj_gene_consensus_db(os.getcwd(),
            vj_db_path, "imgt")

    print("Estimating number of sequences.")
    nseqs = 0

    for csv_filepath in csv_filepaths:
        for i, _ in enumerate(read_csv(csv_filepath)):
            if i >= header_rows:
                break
        for i, row in enumerate(read_csv(csv_filepath)):
            nseqs += 1

    if nseqs > 4200000000:
        raise RuntimeError("Current cap on number of sequences "
                           "per database file is 4.2 billion. "
                           f"You have {nseqs}. Consider splitting "
                           "this up into multiple database files, "
                           "e.g. splitting by species or vgene family.")

    if nseqs == 0:
        print("No sequences found.")
        return

    db_construct_tool = DatabaseConstructionTool(database_filepath,
            numbering_scheme, cdr_definition_scheme,
            "single", receptor_type,
            pid_threshold, license_key, user_email,
            consensus_path, nterm_kmer_dict,
            vj_names, vj_seqs, user_memo)

    settings_list = []
    for expected_key in ["heavy_chain", "light_chain",
            "heavy_numbering", "light_numbering",
            "heavy_vgene", "light_vgene", "species",
            "metadata", "heavy_jgene", "light_jgene"]:
        if expected_key in column_selections:
            settings_list.append(column_selections[expected_key])
        else:
            settings_list.append(-1)

    if max(settings_list) < 0:
        raise RuntimeError("There must be at least some column "
                "selections supplied. Please refer to the docs "
                "for the build_database_from_csv function.")
    if max(settings_list[:2]) < 0:
        raise RuntimeError("Either a heavy column, a light column "
                "or both should be supplied under column_selections.")
    if settings_list[0] >= 0:
        if settings_list[4] < 0 or settings_list[8] < 0:
            raise RuntimeError("If supplying a heavy chain, "
                               "heavy vgenes and jgenes must be "
                               "supplied also.")
    if settings_list[1] >= 0:
        if settings_list[5] < 0 or settings_list[9] < 0:
            raise RuntimeError("If supplying a light chain, "
                               "light vgenes and jgenes must be "
                               "supplied also.")

    print(f"Found {nseqs} sequences. Starting db construction.")

    seqcount = 0
    if reject_file is not None:
        reject_handle = open(reject_file,
                "w+", encoding="utf-8")
    else:
        reject_handle = None

    db_construct_tool.open_transaction()

    for csv_filepath in csv_filepaths:
        for row in read_csv(csv_filepath, skiprows=header_rows,
                            delimiter=delimiter):
            if len(row) == 0:
                continue
            rcode = db_construct_tool.add_full_chain_csv_sequence(
                    row, settings_list)
            if len(rcode) > 0:
                if reject_handle is not None:
                    reject_handle.write(f"{','.join(row)}\t{rcode}\n")
            seqcount += 1
            if seqcount % 10000 == 0:
                db_construct_tool.close_transaction()
                db_construct_tool.open_transaction()
                if verbose:
                    print(f"{seqcount} complete.")

    if reject_handle is not None:
        reject_handle.close()

    if verbose:
        print("Now constructing database indices...")

    db_construct_tool.close_transaction()
    db_construct_tool.open_transaction()

    db_construct_tool.finalize_db_construction(verbose)
    gc.collect()





def build_database_from_cdr_only_csv(csv_filepaths:list,
        database_filepath:str, column_selections:dict, delimiter=',',
        receptor_type = "mab", header_rows:int=1, user_memo:str="",
        reject_file:str = None, verbose:bool=True):
    """TCR data and some antibody data is often stored with cdr3 and/or other
    cdrs pre-extracted. This function builds a database specifically using
    this type of input format where you have already extracted the cdrs. If
    using this function, you MUST use IMGT numbering and cdr definitions.
    Note that the beta chain for TCRs is assumed to be "heavy" and alpha
    is light.

    Args:
        csv_filepaths (list): A list of csv filepaths, which may or
            may not be gzipped.
        database_filepath (str): The desired location and filename
            for the database.
        column_selections (dict): A dictionary which should contain at least
            some of the following keys. Do not include any key for which
            information is not present. Note that if you supply some light
            or heavy chain cdr columns, you MUST supply the cdr3 and vgene
            for that chain type.

            * ``"light_cdr3"``: The number (from 0) of the column in each
              csv file containing light cdr3 (alpha for TCRs).
            * ``"heavy_cdr3"``: The number (from 0) of the column in each
              csv file containing heavy cdr3 (beta for TCRs).
            * ``"light_cdr2"``: The number (from 0) of the column in each
              csv file containing light cdr2 (alpha for TCRs).
            * ``"heavy_cdr2"``: The number (from 0) of the column in each
              csv file containing heavy cdr2 (beta for TCRs).
            * ``"light_cdr1"``: The number (from 0) of the column in each
              csv file containing light cdr1 (alpha for TCRs).
            * ``"heavy_cdr1"``: The number (from 0) of the column in each
              csv file containing heavy cdr1 (beta for TCRs).
            * ``"light_vgene"``: The number (from 0) of the column
              containing vgene assignments for the light chains.
            * ``"heavy_vgene"``: The number (from 0) of the column
              containing vgene assignments for the heavy chains.
            * ``"light_jgene"``: The number (from 0) of the column (if any)
              containing jgene assignments for the light chains. If this is
              supplied, a light_vgene column must be supplied as well.
            * ``"heavy_jgene"``: The number (from 0) of the column
              containing jgene assignments for the heavy chains.
            * ``"species"``: The number (from 0) of the column
              containing species assignments. If not supplied,
              the reader will assume species is always human.
            * ``"metadata"``: The number (from 0) of the column (if any)
              containing metadata for each sequence.

        delimiter (str): The delimiter for the csv files.
        receptor_type (str): One of 'mab', 'tcr'.
        header_rows (int): The number of header rows. Header rows are
            skipped.
        user_memo (str): A string describing the purpose of the database / anything
            important you want your future self or other users to know about the
            contents. Will be saved as part of the database metadata.
        reject_file (str): Either None or a filepath. If None, any sequences that
            are rejected are silently ignored. If a filepath, rejected sequences
            are written to that filepath which is saved as a csv file.
        verbose (bool): If True, print regular updates while running.

    Raises:
        RuntimeError: A RuntimeError is raised if invalid arguments are
            supplied.
    """
    if not isinstance(csv_filepaths, list):
        raise RuntimeError("The csv filepaths should be a list of filepaths. "
                "If you have only one filepath, you can enclose it in brackets "
                "to make a list.")
    if os.path.exists(database_filepath):
        raise RuntimeError("The database already exists.")

    license_key, user_email = get_license_key_info()
    project_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), "..")
    consensus_path = os.path.join(project_path,
            "numbering_tools", "consensus_data")
    nterm_kmer_dict = _load_nterm_kmers()

    vj_db_path = os.path.join(project_path, "vj_tools", "consensus_data")
    vj_names, vj_seqs, _ = load_vj_gene_consensus_db(os.getcwd(),
            vj_db_path, "imgt")

    print("Estimating number of sequences.")
    nseqs = 0

    for csv_filepath in csv_filepaths:
        for i, _ in enumerate(read_csv(csv_filepath)):
            if i >= header_rows:
                break
        for i, row in enumerate(read_csv(csv_filepath)):
            nseqs += 1

    if nseqs > 4200000000:
        raise RuntimeError("Current cap on number of sequences "
                           "per database file is 4.2 billion. "
                           f"You have {nseqs}. Consider splitting "
                           "this up into multiple database files, "
                           "e.g. splitting by species or vgene family.")

    if nseqs == 0:
        print("No sequences found.")
        return

    pairing_type = "single"

    if column_selections["light_cdr3"] >= 0 and \
            column_selections["heavy_cdr3"] >= 0:
        pairing_type = "paired"

    db_construct_tool = DatabaseConstructionTool(database_filepath,
            "imgt", "imgt", pairing_type, receptor_type, 0.7,
            license_key, user_email, consensus_path,
            nterm_kmer_dict, vj_names, vj_seqs, user_memo)

    settings_list = []
    for expected_key in ["heavy_cdr3", "heavy_cdr2", "heavy_cdr1",
            "light_cdr3", "light_cdr2", "light_cdr1",
            "heavy_vgene", "light_vgene", "heavy_jgene",
            "light_jgene", "species", "metadata"]:
        if expected_key in column_selections:
            settings_list.append(column_selections[expected_key])
        else:
            settings_list.append(-1)

    if max(settings_list) < 0:
        raise RuntimeError("There must be at least some column "
                "selections supplied. Please refer to the docs "
                "for the build_database_from_csv function.")

    print(f"Found {nseqs} sequences. Starting db construction.")

    seqcount = 0

    if reject_file is not None:
        reject_handle = open(reject_file, "w+", encoding="utf-8")
    else:
        reject_handle = None

    db_construct_tool.open_transaction()

    for csv_filepath in csv_filepaths:
        for row in read_csv(csv_filepath, skiprows=header_rows,
                            delimiter=delimiter):
            if len(row) == 0:
                continue
            rcode = db_construct_tool.add_extracted_cdr_csv_sequence(
                    row, settings_list)
            if len(rcode) > 0:
                if reject_handle is not None:
                    reject_handle.write(f"{','.join(row)}\t{rcode}\n")
            seqcount += 1
            if seqcount % 10000 == 0:
                db_construct_tool.close_transaction()
                db_construct_tool.open_transaction()
                if verbose:
                    print(f"{seqcount} complete.")

    if reject_handle is not None:
        reject_handle.close()

    if verbose:
        print("Now constructing database indices...")

    db_construct_tool.close_transaction()
    db_construct_tool.open_transaction()

    db_construct_tool.finalize_db_construction(verbose)
    gc.collect()
