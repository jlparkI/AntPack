"""Test database construction for errors.
Use sqlite for testing and prototyping."""
import os
import gzip
import sqlite3
import pytest
import antpack
from antpack import (build_database_from_fasta,
        build_database_from_full_chain_csv,
        build_database_from_cdr_only_csv,
        SingleChainAnnotator, VJGeneTool,
        LocalDBSearchTool)
from antpack.utilities import read_fasta
from .database_utilities import (get_vgene_code,
                setup_canonical_numbering, int_to_bin)
from ..conftest import (get_test_data_filepath, get_test_base_filepath)



def test_local_db_construct(build_local_mab_lmdb):
    """Check that local database construction yields a
    database containing the information we expect."""
    seqs, seqinfos, db_filepath, params, data_dict = build_local_mab_lmdb
    con = sqlite3.connect(db_filepath)
    cur = con.cursor()

    cur.execute("SELECT * FROM database_info;")
    metadata = cur.fetchone()
    assert metadata[0]==params["nmbr_scheme"]
    assert metadata[1]=="mab"
    assert metadata[2]==params["sequence_type"]
    assert metadata[3]==params["cdr_scheme"]
    assert metadata[4]==params["memo"]
    assert metadata[5]==0.7

    # Temporarily create a tool to back-convert compressed
    # sequences. TODO: Fix this so that local db search tool
    # which is tested separately is not required here.
    ldb = LocalDBSearchTool(db_filepath)
    for i, row in enumerate(
            cur.execute("SELECT * FROM sequences;")):
        if not params["mode"].endswith("extract"):
            test_seq = ldb.get_sequence(i+1)[0]
            assert test_seq==seqs[i]
        assert row[1]==seqinfos[i]

    # Check each set of chain tables for expected info.
    for table_code in [0,1]:
        chain_dict = data_dict[table_code]

        for i, (child_id, cdr_grp) in enumerate(
                zip(chain_dict["child_ids"], chain_dict["cdrs"])):
            cur.execute(f"SELECT * FROM _{table_code}_cdrs "
                        f"WHERE child_id = {child_id+1};")
            row_values = cur.fetchone()
            vgene_code = get_vgene_code(chain_dict["vgenes"][i],
                                        chain_dict["vspecies"][i])
            jgene_code = get_vgene_code(chain_dict["jgenes"][i],
                                        chain_dict["vspecies"][i])
            assert len(chain_dict["cdrs"][i][0].
                   replace('-', ''))==int(row_values[0][0])
            assert len(chain_dict["cdrs"][i][1].
                   replace('-', ''))==int(row_values[1][0])
            assert len(chain_dict["cdrs"][i][2].
                    replace('-', ''))==int(row_values[2][0])
            assert row_values[0][1:].decode() + \
                    row_values[1][1:].decode() == \
                    cdr_grp[0] + cdr_grp[1]
            assert row_values[2][1:].decode() == cdr_grp[2]
            assert vgene_code[0][0]==row_values[3]
            assert vgene_code[1]==row_values[4]
            assert vgene_code[2]==row_values[5]
            assert vgene_code[3]==row_values[6]

            assert jgene_code[1]==row_values[7]
            assert jgene_code[3]==row_values[8]
            assert vgene_code[0][3]==row_values[9]
            assert chain_dict["unusual_positions"][i]==\
                    row_values[10]

        tdict = eval_search_table_row_contents(data_dict[table_code]['cdrs'], table_code,
                                       params['nmbr_scheme'], params['cdr_scheme'])
        # 20 here is arbitrary and is determined by the search depth on the db.
        # Update this if that changes. TODO: Make this test independent of arbitrary
        # constants, or fetch the search depth from the db.
        for i in range(20):
            for row in cur.execute(f"SELECT * FROM _trie_{table_code}_{i};"):
                assert row[0] in tdict[i]
                assert tdict[i][row[0]][0] == row[1]
                assert tdict[i][row[0]][1] == row[2]
                assert tdict[i][row[0]][2] == row[3]



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
        fhandle.write("AAAAAAAATTTTTTTTT,VJ123,JJ126,humanoid\n")
        fhandle.write(">NAME,H,XYZ123,alien\n")
        fhandle.write("$$ALT,H,,unknown\n")

    build_database_from_full_chain_csv([input_file], db_path,
            {"heavy_chain":0, "heavy_vgene":1,
             "heavy_jgene":2, "species":3},
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





def eval_search_table_row_contents(cdrs, chain_code,
    numbering_scheme, cdr_scheme):
    """Tests the contents of the rows from the search tables."""
    AAMAP = {k:i for i,k in enumerate("ACDEFGHIKLMNPQRSTVWY-")}
    sca = SingleChainAnnotator(scheme=numbering_scheme)

    if chain_code == 0:
        func_name = numbering_scheme + "_heavy_nmbr_" +  cdr_scheme + "_cdr_arrangement"
    else:
        func_name = numbering_scheme + "_light_nmbr_" +  cdr_scheme + "_cdr_arrangement"

    ordered_nmbr = getattr(antpack.antpack_cpp_ext, func_name)()
    sorted_nmbr = sca.sort_position_codes(ordered_nmbr)
    position_extracts = []
    # 20 here is arbitrary and is determined by the search depth on the db.
    # Update this if that changes. TODO: Make this test independent of arbitrary
    # constants, or fetch the search depth from the db.
    for i in range(len(ordered_nmbr) - 20, len(ordered_nmbr)):
        position_extracts.append(sorted_nmbr.index(ordered_nmbr[i]))

    assert len(position_extracts) == 20

    # Build a trie stored as a dict. This is horribly wasteful space-wise
    # but is good for a test run quickly on a small sequence set, since it
    # is easy to troubleshoot.
    tdict = {k:{} for k in range(20)}
    next_child_code = 1

    

    for cdr_grp in cdrs:
        cdr3_extract = [cdr_grp[-1][e] for e in position_extracts]
        parent_code = 0

        for i, letter in enumerate(cdr3_extract):
            if parent_code not in tdict[i]:
                tdict[i][parent_code] = [next_child_code, 0, 0]
                next_child_code += 21
            if tdict[i][parent_code][2] < 100 or (tdict[i][parent_code][1] &
                                                  (1 << AAMAP[letter])) == 0:
                tdict[i][parent_code][1] = (tdict[i][parent_code][1] | (1 << AAMAP[letter]))
                tdict[i][parent_code][2] += 1
            parent_code = tdict[i][parent_code][0] + AAMAP[letter]

    return tdict






def prep_seqs_for_comparison(numbering_scheme,
        cdr_scheme, chain_type, sequences,
        annotations, annotator,
        vj_tool, cdr_extract_mode):
    """Prep sequences which are of the specified chain
    type for analysis by aligning to the template,
    extracting cdrs etc."""
    selected_seqs = [(s, a) for (s, a) in zip(
        sequences, annotations) if a[2] in chain_type]
    child_ids = [i for i,a in enumerate(annotations)
            if a[2] in chain_type]
    vgenes, jgenes, vjspecies = [], [], []

    for selected_seq in selected_seqs:
        vgene, jgene, _, _, species = vj_tool.assign_vj_genes(
                selected_seq[1], selected_seq[0], "unknown",
                "identity")
        vgenes.append(vgene)
        jgenes.append(jgene)
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
        # If we are extracting cdr3 only, set the other
        # cdrs to all gaps.
        if cdr_extract_mode == "cdr3_extract":
            cdr1 = '-' * len(cdr1)
            cdr2 = '-' * len(cdr2)
        cdrs.append((cdr1, cdr2, cdr3))
        aligned_seqs.append((cdr1 + cdr2 + cdr3, cdr3))

    return aligned_seqs, cdrs, unusual_positions, \
        cdr3_region_len, vgenes, jgenes, vjspecies, child_ids




@pytest.fixture(params=[
    {"filepath":"addtnl_test_data.fasta.gz",
     "nmbr_scheme":"imgt", "cdr_scheme":"imgt",
     "sequence_type":"single", "memo":"testing123",
     "mode":"full_chain"},
    {"filepath":"addtnl_test_data.fasta.gz",
     "nmbr_scheme":"imgt", "cdr_scheme":"north",
     "sequence_type":"unknown", "memo":"testing123",
     "mode":"full_chain"},
    {"filepath":"test_data.csv.gz",
     "nmbr_scheme":"imgt", "cdr_scheme":"imgt",
     "sequence_type":"single", "memo":"testing123",
     "mode":"full_chain"},
    {"filepath":"addtnl_test_data.fasta.gz",
     "nmbr_scheme":"imgt", "cdr_scheme":"imgt",
     "sequence_type":"single", "memo":"testing123",
     "mode":"all_cdr_extract"},
    {"filepath":"addtnl_test_data.fasta.gz",
     "nmbr_scheme":"imgt", "cdr_scheme":"imgt",
     "sequence_type":"single", "memo":"testing123",
     "mode":"cdr3_extract"},
    ])
def build_local_mab_lmdb(tmp_path, get_test_data_filepath,
        request):
    """Builds a local db using a set of arguments
    appropriate for testing LMDB-based local dbs
    and returns all the info the test function will need
    to evaluate this local db."""
    db_filepath = os.path.join(tmp_path, "TEMP_DB")
    temp_csv_filepath = os.path.join(tmp_path, "TEMP_CSV_FILE.csv")
    seqs, seqinfos = [], []

    if request.param["mode"] == "full_chain":
        sca = SingleChainAnnotator(scheme=request.param["nmbr_scheme"])
        vj_tool = VJGeneTool(scheme=request.param["nmbr_scheme"])

        if request.param["filepath"].endswith("test_data.csv.gz"):
            csv_filepath = os.path.join(get_test_data_filepath,
                        request.param["filepath"])
            with gzip.open(csv_filepath, "rt") as fhandle:
                _ = fhandle.readline()
                with open(temp_csv_filepath, 'w') as outhandle:
                    outhandle.write("heavy_sequence,light_sequence,"
                        "heavy_vgene,light_vgene,heavy_jgene,light_jgene,"
                        "heavy_numbering,light_numbering,species\n")
                    for line in fhandle:
                        seq = line.split(',')[0]
                        annotation = sca.analyze_seq(seq)
                        vgene, jgene, _, _, species = vj_tool.assign_vj_genes(
                            annotation, seq, "unknown")
                        if annotation[1] <= 0.7 or vgene == "":
                            continue

                        seqs.append(seq)
                        seqinfos.append("")

                        if annotation[2] == "H":
                            outhandle.write(f"{seq},,{vgene},,"
                                f"{jgene},,{'_'.join(annotation[0])},,"
                                f"{species}\n")
                        else:
                            outhandle.write(f",{seq},,{vgene},,"
                                f"{jgene},,{'_'.join(annotation[0])},"
                                f"{species}\n")

            build_database_from_full_chain_csv(
                [temp_csv_filepath], db_filepath,
                {"heavy_chain":0, "light_chain":1,
                 "heavy_vgene":2, "light_vgene":3,
                 "heavy_jgene":4, "light_jgene":5,
                 "heavy_numbering":6, "light_numbering":7,
                 "species":8},
                numbering_scheme=request.param["nmbr_scheme"],
                cdr_definition_scheme=request.param["cdr_scheme"],
                receptor_type="mab", pid_threshold=0.7,
                user_memo=request.param["memo"],
                reject_file=os.path.join(tmp_path, "REJECTS"))

        else:
            fasta_filepath = os.path.join(get_test_data_filepath,
                    request.param["filepath"])
            build_database_from_fasta([fasta_filepath], db_filepath,
                numbering_scheme=request.param["nmbr_scheme"],
                cdr_definition_scheme=request.param["cdr_scheme"],
                sequence_type=request.param["sequence_type"],
                receptor_type="mab",
                pid_threshold=0.7, user_memo=request.param["memo"],
                reject_file=os.path.join(tmp_path, "fasta_REJECTS"))
            for seqinfo, seq in read_fasta(fasta_filepath):
                seqs.append(seq)
                seqinfos.append(seqinfo)

    else:
        sca = SingleChainAnnotator(scheme="imgt")
        vj_tool = VJGeneTool(scheme="imgt")
        fasta_filepath = os.path.join(get_test_data_filepath,
                    request.param["filepath"])

        with open(temp_csv_filepath, 'w') as outhandle:
            outhandle.write("heavy_cdr3,heavy_cdr2,heavy_cdr1,"
                    "light_cdr3,light_cdr2,light_cdr1,heavy_vgene,"
                    "light_vgene,heavy_jgene,light_jgene,species\n")

            for seqinfo, seq in read_fasta(fasta_filepath):
                annotation = sca.analyze_seq(seq)
                vgene, jgene, _, _, species = vj_tool.assign_vj_genes(
                    annotation, seq, "unknown")
                if annotation[1] <= 0.7 or vgene == "":
                    continue

                seqs.append(seq)
                seqinfos.append("")
                position_labels = sca.assign_cdr_labels(
                        annotation[0], annotation[2])
                cdr3 = ''.join([a for (a,l) in zip(seq,
                            position_labels) if l == "cdr3"])
                cdr2 = ''.join([a for (a,l) in zip(seq,
                            position_labels) if l == "cdr2"])
                cdr1 = ''.join([a for (a,l) in zip(seq,
                            position_labels) if l == "cdr1"])

                if annotation[2] == "H":
                    outhandle.write(f"{cdr3},{cdr2},{cdr1},,,,"
                        f"{vgene},,{jgene},,{species}\n")
                else:
                    outhandle.write(f",,,{cdr3},{cdr2},{cdr1},,"
                        f"{vgene},,{jgene},{species}\n")

        if request.param["mode"] == "all_cdr_extract":
            column_dict = {"heavy_cdr3":0, "heavy_cdr2":1,
                "heavy_cdr1":2, "light_cdr3":3,
                "light_cdr2":4, "light_cdr1":5,
                "heavy_vgene":6, "light_vgene":7,
                "heavy_jgene":8, "light_jgene":9,
                "species":10}
        else:
            column_dict = {"heavy_cdr3":0, "light_cdr3":3,
                "heavy_vgene":6, "light_vgene":7,
                "heavy_jgene":8, "light_jgene":9,
                "species":10}
        build_database_from_cdr_only_csv(
            [temp_csv_filepath], db_filepath,
            column_dict, receptor_type="mab",
            user_memo=request.param["memo"],
            reject_file=os.path.join(tmp_path, "REJECTS"))

    annotations = sca.analyze_seqs(seqs)
    data_dict = {0:{}, 1:{}}
    # Check each set of chain tables for expected info.
    for table_code, chain_coding in enumerate([("H",),
                                    ("L", "K")]):
        seq_prep_results = prep_seqs_for_comparison(
                request.param["nmbr_scheme"],
                request.param["cdr_scheme"],
                chain_coding, seqs,
                annotations, sca, vj_tool,
                request.param["mode"])

        data_dict[table_code]["aligned_seqs"] = seq_prep_results[0]
        data_dict[table_code]["cdrs"] = seq_prep_results[1]
        data_dict[table_code]["unusual_positions"] = seq_prep_results[2]
        data_dict[table_code]["cdr3_region_len"] = seq_prep_results[3]
        data_dict[table_code]["vgenes"] = seq_prep_results[4]
        data_dict[table_code]["jgenes"] = seq_prep_results[5]
        data_dict[table_code]["vspecies"] = seq_prep_results[6]
        data_dict[table_code]["child_ids"] = seq_prep_results[7]
        # Little bit of a hack, but...the cdrs in this test data
        # do not have any unusual inserts, while framework occasionally
        # does. If we are extracting cdrs only, set all unusual
        # position values to 0.
        if request.param["mode"].endswith("extract"):
            data_dict[table_code]["unusual_positions"] = [0]*len(
                    data_dict[table_code]["unusual_positions"])
    return seqs, seqinfos, db_filepath, request.param, data_dict
