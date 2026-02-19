"""Test database management & search."""
import os
import gzip
from copy import deepcopy
import random
from math import floor
import pytest
import numpy as np
from antpack import (build_database_from_fasta,
        build_database_from_full_chain_csv,
        SingleChainAnnotator, LocalDBSearchTool, VJGeneTool)
from antpack.antpack_cpp_ext import (SequenceTemplateAligner,
        return_imgt_canonical_numbering_cpp,
        get_blosum_mismatch_matrix_cpp)
from antpack.utilities import read_fasta
from .database_utilities import (get_vgene_code,
                        setup_canonical_numbering)
from ..conftest import (get_test_data_filepath,
                        standard_aa_list)



@pytest.mark.parametrize("thread_scheme", [False, True])
def test_local_db_search(build_local_mab_lmdb,
        standard_aa_list, std_blosum_matrix, thread_scheme):
    """Check that searches retrieve the sequences that
    we would expect based on a brute-force exact search."""
    seqs, seqinfos, db_filepath, msa, msa_codes, _ = \
            build_local_mab_lmdb
    print(f"Testing {thread_scheme} scheme...")
    local_db = LocalDBSearchTool(db_filepath, thread_scheme)

    # Check that get_num_seqs works correctly.
    assert local_db.get_num_seqs("all") == len(seqs)
    assert (local_db.get_num_seqs("heavy") ==
            len([s for s in msa if s[4] == "H"]))
    assert (local_db.get_num_seqs("light") ==
            len([s for s in msa if s[4] != "H"]))

    aamap = {k:i for i,k in enumerate(standard_aa_list)}

    # Check that vgene reassembly works correctly.
    for i, seqinfo in enumerate(msa):
        if seqinfo[4] == "H":
            chain_type = "heavy"
        else:
            chain_type = "light"
        test_vgene, test_jgene = local_db.get_vgene_jgene(i+1, chain_type)
        assert test_vgene == seqinfo[1].split("_")[0]
        assert test_jgene == seqinfo[2].split("_")[0]

    # Randomly select a sequence, randomly mutate it and
    # find its closest matches using a random setting for
    # distance and other parameters. Check it on the fly
    # against the search tool. Simultaneously store the
    # results grouped by search criteria so that we can
    # later test all of the results against the batch
    # search.
    allowed_search_settings = {
            "cdr_cutoff":[0.2, 0.25, 0.3],
            "blosum_cutoff":[-1,4],
            "search_mode":["123", "3"],
            "cdr_length_shift":[0,1,2],
            "symmetric_search":[True, False],
            "use_vgene_family_only":[True,False],
            "use_vgene":[True,False],
            "use_jgene":[False]
        }

    random.seed(123)

    for ctr in range(1000):
        idx = random.randint(0, len(msa) - 1)
        query_seq = list(msa[idx][0][0])
        if msa[idx][4] == "H":
            codes = msa_codes[0]
        else:
            codes = msa_codes[1]
        mut_percentage = random.uniform(0, 0.2)

        for i, qletter in enumerate(query_seq):
            if qletter == '-':
                continue
            if codes[1][i].startswith("cdr") and \
                random.uniform(0,1) < mut_percentage:
                query_seq[i] = random.choice(standard_aa_list[:-1])

        query_seq = "".join(query_seq)
        search_settings = {k:random.choice(options)
                for k,options in allowed_search_settings.items()}

        if search_settings["use_vgene"]:
            vgene = msa[idx][1]
            species = msa[idx][3]
            vgene_filter = get_vgene_code(vgene, species)
            if search_settings["use_vgene_family_only"]:
                vgene_filter[0][2] = 255
            if search_settings["use_jgene"]:
                jgene = msa[idx][2]
                species = msa[idx][3]
                jgene_filter = get_vgene_code(jgene, species)
            else:
                jgene = ""
                jgene_filter = ()
        else:
            vgene, species = "", ""
            vgene_filter = ()
            jgene = ""
            jgene_filter = ()

        hits = local_db.search(query_seq,
                (codes[0], 1, msa[idx][4], ""),
                search_settings["search_mode"],
                search_settings["cdr_cutoff"],
                search_settings["blosum_cutoff"],
                search_settings["cdr_length_shift"],
                search_settings["use_vgene_family_only"],
                search_settings["symmetric_search"],
                vgene, species, jgene)[0]
        hits = sorted(hits, key=lambda x: (x[1], x[0], x[2]))

        gt_hit_idx = perform_exact_search(query_seq, msa, msa[idx][4],
                codes, search_settings, vgene_filter,
                jgene_filter, aamap, std_blosum_matrix)
        # If using a BLOSUM cutoff, distance will be floating
        # point, check result is close. Otherwise check for
        # exact match.
        if search_settings["blosum_cutoff"] > 0:
            assert ([h[0] for h in hits]==
                            [h[0]for h in gt_hit_idx])
            assert np.allclose([h[1] for h in hits],
                            [h[1] for h in gt_hit_idx])

        else:
            assert hits==gt_hit_idx

        # Make sure the metadata and sequence retrieved
        # from the database for a given id code match those
        # in the input.
        for hit in hits:
            db_seq, db_metadata = local_db.get_sequence(hit[0])
            assert db_seq==seqs[hit[0]-1]
            assert seqinfos[hit[0]-1]==db_metadata



def test_local_db_search_setup(build_local_mab_lmdb,
    standard_aa_list):
    """Check that the local db management tool is set
    up correctly."""
    _, _, db_filepath, msa, msa_codes, params = \
        build_local_mab_lmdb
    local_db = LocalDBSearchTool(db_filepath)

    # Check that the search tool retrieves the correct metadata.
    metadata = local_db.get_database_metadata()
    assert params["nmbr_scheme"] == metadata[0]
    assert 'mab' == metadata[1]
    assert 'single' == metadata[2]
    assert 'imgt' == metadata[3]
    assert '' == metadata[4]

    # Next check that the position counts table matches
    # expected. First check for heavy. The database
    # returns the tables all together but the Python
    # test code will check them separately.
    dtables = local_db.get_database_counts()

    sub_msa = [m[0][0] for m in msa if m[4] == "H"]
    gt_dtable = get_kmer_counts(sub_msa, msa_codes[0][1],
                                standard_aa_list)

    for key, value in gt_dtable.items():
        assert dtables[0][key]==value

    sub_msa = [m[0][0] for m in msa if m[4] != "H"]
    gt_dtable = get_kmer_counts(sub_msa, msa_codes[1][1],
                                standard_aa_list)

    for key, value in gt_dtable.items():
        assert dtables[1][key]==value





@pytest.mark.parametrize("cdr_scheme", ["imgt"])
def test_basic_clustering(tmp_path_factory,
    get_test_data_filepath, cdr_scheme):
    """Builds a local db and clusters it using a set
    of data with known ground truths."""
    temp_folder = tmp_path_factory.mktemp("clustering")
    input_file = os.path.join(get_test_data_filepath,
                "prepped_clustering_test.csv.gz")
    db_path = os.path.join(temp_folder, "cluster_db")
    reject_file = os.path.join(temp_folder, "REJECTS")

    build_database_from_full_chain_csv([input_file], db_path,
            {"heavy_chain":4, "heavy_vgene":0,
             "heavy_jgene":1, "species":7},
            numbering_scheme='imgt',
            cdr_definition_scheme=cdr_scheme,
            receptor_type="mab", pid_threshold=0.7,
            user_memo="", reject_file=reject_file,
            header_rows=1)

    local_db = LocalDBSearchTool(db_path)
    test_assignments = local_db.basic_clustering(
            "heavy", "3", 0.2, -1, True)
    with gzip.open(input_file, "rt") as fh:
        _=fh.readline()
        actual_assignments = [int(l.strip().split(',')[8]) for
                              l in fh]
    assert actual_assignments == test_assignments







def perform_exact_search(query, msa, chain_code, msa_codes,
        search_params, vgene_filter, jgene_filter, aamap,
        blosum_matrix):
    """Performs an exact search on an input msa, finding
    all sequences that have hamming distance <= specified.
    Additionally filter for BLOSUM distance (if requested)."""
    retained_dists, hit_idx, noncanon_pos, cdrlen = [], [], [], []
    max_hamming = 0

    for j, region in enumerate(["cdr1", "cdr2", "cdr3"]):
        cdrlen.append(len([q for (q,l) in zip(query, msa_codes[1])
            if q != '-' and l == region]))

    max_hamming = [floor((cdrlen[0] + cdrlen[1]) *
                    search_params["cdr_cutoff"]),
        floor(cdrlen[2] * search_params["cdr_cutoff"])]

    for i, seq_data in enumerate(msa):
        if seq_data[4] != chain_code:
            continue

        cdr_dists = [0,0]
        region_lengths = [0,0]

        if len(vgene_filter) > 0:
            if vgene_filter[0][0] > 0 and seq_data[5][0][0] != \
                        vgene_filter[0][0]:
                continue
            if seq_data[5][1] != vgene_filter[1]:
                continue
            if vgene_filter[0][2] != 255 and seq_data[5][2] != \
                    vgene_filter[2]:
                continue
            if vgene_filter[0][1] != 255 and seq_data[5][0][3] != \
                    vgene_filter[0][3]:
                continue

            if len(jgene_filter) > 0:
                if seq_data[6][1] != jgene_filter[1]:
                    continue

        for j, region in enumerate(["cdr1", "cdr2", "cdr3"]):
            if str(j+1) not in search_params["search_mode"]:
                continue
            region_dist = len([a for (a,b,l) in zip(seq_data[0][0], query,
                msa_codes[1]) if l==region and a != b])
            cdr_dists[int(j/2)] += region_dist
            region_len = len([s for (s,l) in
                    zip(seq_data[0][0], msa_codes[1])
                    if l==region and s != '-'])
            region_lengths[int(j/2)] += region_len

            # Add an arbitrary large number if the cdr length is
            # unacceptable or if an individual distance is
            # unacceptable so that the sequence is not saved.
            if region == "cdr3":
                if region_len > cdrlen[j] + \
                        search_params["cdr_length_shift"]:
                    cdr_dists[int(j/2)] += 200
                if region_len < cdrlen[j] - \
                        search_params["cdr_length_shift"]:
                    cdr_dists[int(j/2)] += 200

        # If symmetric search was specified, adjust the cutoff.
        hamming_cutoffs = deepcopy(max_hamming)
        if search_params["symmetric_search"]:
            if region_lengths[0] > 0:
                hamming_cutoffs[0] = min(max_hamming[0],
                        floor(search_params["cdr_cutoff"] * region_lengths[0]))
            hamming_cutoffs[1] = min(max_hamming[1],
                        floor(search_params["cdr_cutoff"] * region_lengths[1]))

        # If we meet Hamming distance criteria, store sequence
        # as a hit UNLESS BLOSUM cutoff was also specified, in
        # which case calculate BLOSUM distance and see if we
        # also meet THAT cutoff.
        if cdr_dists[0] <= hamming_cutoffs[0] and \
                cdr_dists[1] <= hamming_cutoffs[1]:
            if search_params["blosum_cutoff"] < 0:
                hit_idx.append(i+1)
                retained_dists.append(cdr_dists[0] + cdr_dists[1])
                noncanon_pos.append(seq_data[0][1])
            else:
                blosum_dist, max_blosum_dist = 0, 0
                allowed_regions = {f"cdr{k}" for k in
                                   search_params["search_mode"]}
                for l1, l2, code in zip(seq_data[0][0], query,
                                            msa_codes[1]):
                    if code not in allowed_regions:
                        continue
                    l1num = aamap[l1] * 22
                    l2num = aamap[l2]
                    max_blosum_dist = max(blosum_matrix[l1num + l2num],
                                          max_blosum_dist)
                    blosum_dist += blosum_matrix[l1num + l2num]
                blosum_dist = float(blosum_dist)
                if max_blosum_dist <= search_params["blosum_cutoff"]:
                    hit_idx.append(i+1)
                    retained_dists.append(blosum_dist)
                    noncanon_pos.append(seq_data[0][1])


    if len(hit_idx) == 0:
        return []

    hit_idx = list(zip(hit_idx, retained_dists, noncanon_pos))
    return sorted(hit_idx, key=lambda x: (x[1], x[0], x[2]))


def get_kmer_counts(msa, cdr_labels, aalist):
    """Extract kmers from the cdrs of an msa and
    count the number present of each."""
    dimer_count_tables = {}
    aamap = {k:i for i,k in enumerate(aalist)}

    cdr_start = cdr_labels.index('cdr3')
    cdr_end = cdr_labels.index('fmwk4')
    ndimers = cdr_end - cdr_start - 1
    for m in msa:
        cdr = m[cdr_start:cdr_end]
        cdrlen = len([a for a in cdr if a != '-'])
        if cdrlen not in dimer_count_tables:
            dimer_count_tables[cdrlen] = [0] * ndimers * 21 * 21
        for i in range(ndimers):
            kmer_code = aamap[cdr[i]] * 21 + aamap[cdr[i+1]]
            dimer_count_tables[cdrlen][kmer_code+i*21*21] += 1

    return dimer_count_tables



@pytest.fixture
def std_blosum_matrix():
    """Returns the standard BLOSUM matrix for mismatches."""
    return get_blosum_mismatch_matrix_cpp()



@pytest.fixture(scope="module", params=[
    {"filepath":"addtnl_test_data.fasta.gz",
     "nmbr_scheme":"imgt", "cdr_scheme":"imgt"},
    ])
def build_local_mab_lmdb(tmp_path_factory,
        get_test_data_filepath, request):
    """Builds a local db using a set of arguments
    appropriate for testing LMDB-based local dbs
    and returns all the info the test function will need
    to evaluate this local db."""
    temp_folder = tmp_path_factory.mktemp("mgmt_db_folder")
    db_filepath = os.path.join(temp_folder,
                f"{request.param['nmbr_scheme']}_"
                f"{request.param['cdr_scheme']}_TEMP_DB")

    fasta_filepath = os.path.join(get_test_data_filepath,
                request.param["filepath"])
    build_database_from_fasta([fasta_filepath],
        db_filepath, numbering_scheme=request.param["nmbr_scheme"],
        cdr_definition_scheme=request.param["cdr_scheme"],
        sequence_type="single", receptor_type="mab",
        pid_threshold=0.7, user_memo="",
        reject_file=None)

    seqs, seqinfos = [], []
    for seqinfo, seq in read_fasta(fasta_filepath):
        seqs.append(seq)
        seqinfos.append(seqinfo)

    heavy_aligner, heavy_codes = setup_canonical_numbering(
        request.param["nmbr_scheme"], request.param["cdr_scheme"], "H")
    light_aligner, light_codes = setup_canonical_numbering(
        request.param["nmbr_scheme"], request.param["cdr_scheme"], "L")
    sca = SingleChainAnnotator(scheme=request.param["nmbr_scheme"])
    vj = VJGeneTool(scheme=request.param["nmbr_scheme"])

    annotations = sca.analyze_seqs(seqs)
    heavy_msa, light_msa = [], []

    # Note that we must also track the number of
    # unusual (i.e. noncanonical) positions. This
    # info is included in the msa which is therefore
    # stored as a tuple for each sequence.
    for annotation, seq in zip(annotations, seqs):
        if annotation[2] == "H":
            aligned_seq, noncanon = \
                heavy_aligner.align_sequence_save_unusual_positions(
                seq, annotation[0])
            heavy_msa.append( (aligned_seq, len(noncanon)) )
        else:
            aligned_seq, noncanon = \
                light_aligner.align_sequence_save_unusual_positions(
                seq, annotation[0])
            light_msa.append( (aligned_seq, len(noncanon)) )

    heavy_codes = (heavy_codes, sca.assign_cdr_labels(heavy_codes,
        "H", request.param["cdr_scheme"]))
    light_codes = (light_codes, sca.assign_cdr_labels(light_codes,
        "L", request.param["cdr_scheme"]))

    msa, hctr, lctr = [], 0, 0
    for seq, annotation in zip(seqs, annotations):
        vgene, jgene, _, _, species = vj.assign_vj_genes(annotation, seq,
                "unknown", "identity")
        if annotation[2] == "H":
            msa.append( (heavy_msa[hctr], vgene, jgene,
                         species, annotation[2],
                         get_vgene_code(vgene, species),
                         get_vgene_code(jgene, species)) )
            hctr += 1
        else:
            msa.append( (light_msa[lctr], vgene, jgene,
                         species, annotation[2],
                         get_vgene_code(vgene, species),
                         get_vgene_code(jgene, species)) )
            lctr += 1

    return seqs, seqinfos, db_filepath, msa, \
            (heavy_codes, light_codes), request.param
