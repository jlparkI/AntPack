"""Test database construction for errors.
This test currently depends on one method for the
LocalDBTool which is tested separately --
this is unfortunate but was adopted
because existing wrappers for RocksDB were hard
to use for what we are doing here, so it's easier
to use LocalDatabaseTool. TODO: Create a separate
lightweight wrapper for RocksDB for testing."""
import os
import gzip
import pytest
import numpy as np
from antpack import (build_database_from_fasta,
        build_database_from_csv,
        build_tcr_database_from_csv,
        SingleChainAnnotator, VJGeneTool,
        LocalDBTool)
from antpack.utilities import read_fasta
from .database_utilities import get_vgene_code, setup_canonical_numbering



def test_local_db_construct(build_local_mab_lmdb):
    """Check that local database construction yields a
    database containing the information we expect."""
    seqs, seqinfos, db_filepath, params = build_local_mab_lmdb
    ldb = LocalDBTool(db_filepath, 'single')

    assert (ldb.local_db_manager.retrieve_raw_kv_test("database_info",
                b"numbering_scheme").decode())==params["nmbr_scheme"]
    assert (ldb.local_db_manager.retrieve_raw_kv_test("database_info",
                b"sequence_type").decode())==params["sequence_type"]
    assert (ldb.local_db_manager.retrieve_raw_kv_test("database_info",
            b"cdr_definition_scheme").decode())==params["cdr_scheme"]
    assert (ldb.local_db_manager.retrieve_raw_kv_test("database_info",
                b"receptor_type").decode())=="mab"
    assert (ldb.local_db_manager.retrieve_raw_kv_test("database_info",
                b"user_memo").decode())==params["memo"]

    for i, (seq, seqinfo) in enumerate(zip(seqs, seqinfos)):
        ck_key = i.to_bytes(4, byteorder='big')
        assert seq==ldb.local_db_manager.retrieve_raw_kv_test(
                "sequences", ck_key).decode()
        assert seqinfo==ldb.local_db_manager.retrieve_raw_kv_test(
                "metadata", ck_key).decode()

    sca = SingleChainAnnotator(scheme=params["nmbr_scheme"])
    annotations = sca.analyze_seqs(seqs)
    vj_tool = VJGeneTool(scheme=params["nmbr_scheme"])

        # Check each set of chain tables for expected info.
    for chain_code, chain_symbol in enumerate([("H",), ("L", "K")]):
        aligned_seqs, cdrs, unusual_positions, cdr3_region_len, vgenes, vspecies, child_ids = \
                prep_seqs_for_comparison(params["nmbr_scheme"],
                    params["cdr_scheme"], chain_symbol, seqs,
                    annotations, sca, vj_tool)

        for i, (child_id, cdr_grp) in enumerate(
                zip(child_ids, aligned_seqs)):
            ck_key = child_id.to_bytes(4, byteorder='big')
            value = ldb.local_db_manager.retrieve_raw_kv_test(
                f"_{chain_code}_cdrs", ck_key)
            vgene_code = get_vgene_code(vgenes[i], vspecies[i])
            assert vgene_code[0]==int(value[4])
            assert vgene_code[1]==int(value[5])
            assert vgene_code[2]==int(value[6])
            assert vgene_code[3]==int(value[7])
            assert unusual_positions[i]==int(value[3])
            assert cdr_grp[0]==value[12:].decode()

        kmer_profile = eval_nmbr_table_row_contents(ldb,
                    cdrs, child_ids, vgenes, vspecies,
                    chain_code)

        for chain_code, code_table in kmer_profile.items():
            for cdrlen, array in code_table.items():
                table = f"_{chain_code}_diversity"
                key = cdrlen.to_bytes(4, byteorder="big")
                raw_buff = ldb.local_db_manager.retrieve_raw_kv_test(
                        table, key)
                raw_buff = [int.from_bytes(raw_buff[k:k+4], byteorder="big")
                            for k in range(0, len(raw_buff), 4)]
                test_arr = np.array(raw_buff, dtype=np.uint32)
                test_arr = test_arr.reshape((array.shape[0],
                                             array.shape[1]))
                assert np.allclose(test_arr, array)



def test_low_quality_seqs(tmp_path):
    """Test what happens if we try to write low-quality
    sequences to the database."""
    input_file = os.path.join(tmp_path, "temp_data_file.fa")
    db_path = os.path.join(tmp_path, "TEMP_DB")
    reject_file = os.path.join(tmp_path, "REJECTS")

    with open(input_file, "w+", encoding="utf-8") as fhandle:
        fhandle.write(">NAME\nAAAAAAAATTTTTTTTT\n")
        fhandle.write(">NAME\ntesting123\n")
        fhandle.write(">NAME\nEVQLEVQLEVQL\n")

    build_database_from_fasta([input_file], db_path,
        numbering_scheme="imgt",
        cdr_definition_scheme="imgt",
        sequence_type="paired", receptor_type="mab",
        pid_threshold=0.7, user_memo="test",
        reject_file=reject_file)

    # There should be no error in writing the db;
    # we expect that the sequences that failed the
    # threshold test will be written to the reject file,
    # and the relevant database tables will be empty.
    with open(reject_file, "r", encoding="utf-8") as reject_handle:
        with open(input_file, "r", encoding="utf-8") as \
                fhandle:
            for line1, line2 in zip(reject_handle, fhandle):
                assert line1.strip().split('\t')[0]==line2.strip() 

    input_file = os.path.join(tmp_path, "temp_data_file.csv")
    db_path = os.path.join(tmp_path, "TEMP_DB2")
    reject_file = os.path.join(tmp_path, "REJECTS2")

    with open(input_file, "w+") as fhandle:
        fhandle.write("AAAAAAAATTTTTTTTT,H,VJ123,humanoid\n")
        fhandle.write(">NAME,H,XYZ123,alien\n")
        fhandle.write("$$ALT,H,,unknown\n")

    build_database_from_csv([input_file], db_path,
            {"heavy_chain":0, "heavy_chaintype":1,
                "heavy_vgene":2, "species":3},
            numbering_scheme='imgt',
            cdr_definition_scheme='imgt',
            receptor_type="mab", pid_threshold=0.7,
            user_memo="", reject_file=reject_file,
            header_rows=0)

    # There should be no error in writing the db;
    # we expect that the sequences that failed the
    # threshold test will be written to the reject file,
    # and the relevant database tables will be empty.
    with open(reject_file, "r", encoding="utf-8") as reject_handle:
        with open(input_file, "r", encoding="utf-8") as \
                fhandle:
            for line1, line2 in zip(reject_handle, fhandle):
                assert line1.strip().split('\t')[0]==line2.strip() 


def test_tcr_fmt_loading(get_test_data_filepath, tmp_path):
    """Check that tcr standard format data is loaded
    as expected."""
    csv_filepath = os.path.join(get_test_data_filepath,
            "non_antibody_test_data", "tcr_simpletest.csv.gz")
    db_filepath = os.path.join(tmp_path, "TEMP_DB")

    expected_cdrs, expected_mainseq = [], []
    build_tcr_database_from_csv(
            [csv_filepath], db_filepath,
            {"beta_cdr3":0, "beta_vgene":1,
             "beta_jgene":2, "species":4},
            user_memo="", reject_file=None)
    with gzip.open(csv_filepath, "rt") as fhandle:
        _ = fhandle.readline()
        for line in fhandle:
            elements = line.split(",")
            expected_mainseq.append(elements[0] +
                "," + elements[1] + "," +
                elements[2])
            expected_cdrs.append(elements[3])

    ldb = LocalDBTool(db_filepath, "single")

    for i, expected_cdr in enumerate(expected_cdrs):
        ck_key = i.to_bytes(4, byteorder='big')
        value = ldb.local_db_manager.retrieve_raw_kv_test(
            "_0_cdrs", ck_key)
        assert expected_cdr==value[12:].decode()

    for i, expected_seq in enumerate(expected_mainseq):
        ck_key = i.to_bytes(4, byteorder='big')
        value = ldb.local_db_manager.retrieve_raw_kv_test(
            "sequences", ck_key)
        assert expected_seq==value.decode()




def eval_nmbr_table_row_contents(local_db_tool, cdrs, child_ids,
        vgenes, vspecies, chain_code, numbering_scheme="imgt"):
    """Tests the contents of the rows from the numbering
    table."""
    AAMAP = {k:i for i,k in enumerate("ACDEFGHIKLMNPQRSTVWY-")}
    letter_position_map = {"A":[48], "C":[32],
                           "D":[16], "E":[16],
                           "F":[0], "G":[48,32],
                           "H":[48,16], "I":[48,0],
                           "K":[32,16], "L":[32,0],
                           "M":[32], "N":[16,0],
                           "P":[32], "Q":[16,0],
                           "R":[32,16], "S":[48,32,16],
                           "T":[48,32,0], "V":[48,16,0],
                           "W":[32,16,0], "Y":[48,32,16,0],
                           "-":[], "X":[32]
                           }
    profile_counts = {}

    for i, (child_id, cdr_group) in enumerate(zip(child_ids, cdrs)):
        cdr = cdr_group[2]
        cdr3len = len(cdr.replace('-', ''))
        vgene_code = get_vgene_code(vgenes[i], vspecies[i])

        if chain_code not in profile_counts:
            profile_counts[chain_code] = {}

        if cdr3len not in profile_counts[chain_code]:
            profile_counts[chain_code][cdr3len] = \
                    np.zeros(( len(cdr) - 1, 21*21))

        # Construct the augmented child id.
        if numbering_scheme in ("imgt", "aho"):
            cdr_extract = cdr[:8] + cdr[-8:]
        else:
            cdr_extract = cdr[:16]

        bytestring = ['0' for j in range(64)]
        for j, letter in enumerate(cdr_extract):
            for offset_pos in letter_position_map[letter]:
                bytestring[j+offset_pos] = '1'

        # The database is designed to work regardless of endianness,
        # but it seems like overkill to take the same precautions
        # in the unit test.
        # These unit tests are always run on little-endian machines.
        unsigned_tag_filter = int(''.join(bytestring)[::-1], 2).to_bytes(
                8, byteorder='big')

        for position in range(0, len(cdr) - 1):
            kmer = cdr[position:position+2]
            codeval = AAMAP[kmer[0]] * 21 + AAMAP[kmer[1]]
            profile_counts[chain_code][cdr3len][position, codeval] += 1
            if kmer == "--":
                continue

            tag = AAMAP[kmer[0]].to_bytes() + AAMAP[kmer[1]].to_bytes() + \
                cdr3len.to_bytes() + vgene_code[0].to_bytes() + \
                vgene_code[1].to_bytes() + vgene_code[2].to_bytes() + \
                vgene_code[3].to_bytes() + child_id.to_bytes(4, byteorder='big')
            table = f"_{chain_code}_{position}_len_dimer_vgene_key"
            alt_filter = local_db_tool.local_db_manager.retrieve_raw_kv_test(
                    table, tag)
            assert alt_filter==unsigned_tag_filter

    return profile_counts






def prep_seqs_for_comparison(numbering_scheme,
        cdr_scheme, chain_type, sequences,
        annotations, annotator,
        vj_tool):
    """Prep sequences which are of the specified chain
    type for analysis by aligning to the template,
    extracting cdrs etc."""
    selected_seqs = [(s, a) for (s, a) in zip(
        sequences, annotations) if a[2] in chain_type]
    child_ids = [i for i,a in enumerate(annotations)
            if a[2] in chain_type]
    vgenes, vjspecies = [], []

    for selected_seq in selected_seqs:
        vgene, _, _, _, species = vj_tool.assign_vj_genes(
                selected_seq[1], selected_seq[0], "unknown",
                "identity")
        vgenes.append(vgene)
        vjspecies.append(species)

    template_aligner, canon_nmbr = setup_canonical_numbering(
            numbering_scheme, cdr_scheme, chain_type[0])
    recognized_positions = set(canon_nmbr)
    recognized_positions.add('-')
    cdr_labels = annotator.assign_cdr_labels(canon_nmbr,
            chain_type[0], cdr_scheme)

    cdr3_codes = [a for (a,c) in zip(canon_nmbr, cdr_labels)
            if c == "cdr3"]
    cdr3_region_len = len(cdr3_codes)

    aligned_seqs, cdrs, unusual_positions = \
            [], [], []

    for selected_seq in selected_seqs:
        aligned_seq = template_aligner.align_sequence(
                selected_seq[0], selected_seq[1][0], False)
        unusual_positions.append( len([a for a in selected_seq[1][0]
                if a not in recognized_positions]) )
        cdr1 = ''.join([a for (a,c) in zip(aligned_seq,
            cdr_labels) if c == "cdr1"])
        cdr2 = ''.join([a for (a,c) in zip(aligned_seq,
            cdr_labels) if c == "cdr2"])
        cdr3 = ''.join([a for (a,c) in zip(aligned_seq,
            cdr_labels) if c == "cdr3"])
        cdrs.append((cdr1, cdr2, cdr3))
        aligned_seqs.append((cdr1 + cdr2 + cdr3, cdr3))

    return aligned_seqs, cdrs, unusual_positions, \
        cdr3_region_len, vgenes, vjspecies, child_ids




@pytest.fixture(params=[
    {"filepath":"addtnl_test_data.fasta.gz",
     "nmbr_scheme":"imgt", "cdr_scheme":"imgt",
     "sequence_type":"single", "memo":"testing123"},
    {"filepath":"addtnl_test_data.fasta.gz",
     "nmbr_scheme":"imgt", "cdr_scheme":"north",
     "sequence_type":"unknown", "memo":"testing123"},
    {"filepath":"db_csv_load_testing.csv.gz",
     "nmbr_scheme":"imgt", "cdr_scheme":"imgt",
     "sequence_type":"single", "memo":"testing123"},
    ])
def build_local_mab_lmdb(tmp_path, get_test_data_filepath,
        request):
    """Builds a local db using a set of arguments
    appropriate for testing LMDB-based local dbs
    and returns all the info the test function will need
    to evaluate this local db."""
    db_filepath = os.path.join(tmp_path, "TEMP_DB")
    seqs, seqinfos = [], []

    if "csv" in request.param["filepath"]:
        csv_filepath = os.path.join(get_test_data_filepath,
                        request.param["filepath"])
        # Using hardcoded column selections here (since the
        # test data is always the same)
        build_database_from_csv([csv_filepath], db_filepath,
                {"heavy_chain":0, "light_chain":3,
                    "heavy_chaintype":2, "light_chaintype":5,
                    "heavy_numbering":6, "light_numbering":7},
                numbering_scheme=request.param["nmbr_scheme"],
                cdr_definition_scheme=request.param["cdr_scheme"],
                receptor_type="mab", pid_threshold=0.7,
                user_memo=request.param["memo"], reject_file=None)
        with gzip.open(csv_filepath, "rt") as fhandle:
            _ = fhandle.readline()
            for line in fhandle:
                seqinfos.append("")
                elements = line.split(",")
                if len(elements[0]) == 0:
                    seqs.append(elements[3])
                else:
                    seqs.append(elements[0])

    else:
        fasta_filepath = os.path.join(get_test_data_filepath,
                    request.param["filepath"])
        build_database_from_fasta([fasta_filepath], db_filepath,
            numbering_scheme=request.param["nmbr_scheme"],
            cdr_definition_scheme=request.param["cdr_scheme"],
            sequence_type=request.param["sequence_type"],
            receptor_type="mab",
            pid_threshold=0.7, user_memo=request.param["memo"],
            reject_file=None)
        for seqinfo, seq in read_fasta(fasta_filepath):
            seqs.append(seq)
            seqinfos.append(seqinfo)

    return seqs, seqinfos, db_filepath, request.param
