"""Test database management & search."""
import os
import shutil
import random
import unittest
import numpy as np
from scipy.sparse import csr_matrix
from antpack import (build_database_from_fasta,
        SingleChainAnnotator, LocalDBTool, VJGeneTool)
from antpack.antpack_cpp_ext import (SequenceTemplateAligner,
        return_imgt_canonical_numbering_cpp,
        get_blosum_mismatch_matrix_cpp)
from antpack.utilities import read_fasta


AAMAP = {k:i for i,k in enumerate("ACDEFGHIKLMNPQRSTVWY-")}
BLOSUM = get_blosum_mismatch_matrix_cpp()



class TestLocalDBManagement(unittest.TestCase):


    def test_local_db_search(self):
        """Check that searches retrieve the sequences that
        we would expect based on a brute-force exact search."""
        current_dir = os.path.abspath(os.path.dirname(
            __file__))
        data_filepath = os.path.join(current_dir,
                "test_data", "addtnl_test_data.fasta.gz")
        nmbr_scheme = "imgt"
        cdr_scheme = "imgt"

        if not os.path.exists('TEMP_DB.db'):
            build_database_from_fasta([data_filepath],
                "TEMP_DB.db", "TEMP_FILE",
                numbering_scheme=nmbr_scheme,
                cdr_definition_scheme=cdr_scheme,
                sequence_type="single", receptor_type="mab",
                pid_threshold=0.7, user_memo="")

        seqs, seqinfos = [], []

        for seqinfo, seq in read_fasta(data_filepath):
            seqs.append(seq)
            seqinfos.append(seqinfo)

        for thread_scheme in ["single", "multi"]:
            print(f"Testing {thread_scheme} scheme...")
            local_db = LocalDBTool("TEMP_DB.db", thread_scheme)

            # First annotate the chains and build a heavy and light
            # MSA, then interleave them.
            # Use this msa to do brute force searches and compare
            # the results with those we get back from the database.
            sca = SingleChainAnnotator(scheme=nmbr_scheme)
            vj = VJGeneTool(scheme=nmbr_scheme)

            annotations = sca.analyze_seqs(seqs)
            heavy_msa = [(a,s) for (a,s) in zip(annotations, seqs)
                    if a[2]=="H"]
            light_msa = [(a,s) for (a,s) in zip(annotations, seqs)
                    if a[2]!="H"]
            heavy_codes, heavy_msa = sca.build_msa([t[1] for t in
                heavy_msa], [t[0] for t in heavy_msa])
            light_codes, light_msa = sca.build_msa([t[1] for t in
                light_msa], [t[0] for t in light_msa])

            # Check that get_num_seqs works correctly.
            self.assertTrue(local_db.get_num_seqs("all") == len(seqs))
            self.assertTrue(local_db.get_num_seqs("heavy") ==
                    len(heavy_msa))
            self.assertTrue(local_db.get_num_seqs("light") ==
                    len(light_msa))

            heavy_codes = (heavy_codes, sca.assign_cdr_labels(heavy_codes,
                "H", cdr_scheme))
            light_codes = (light_codes, sca.assign_cdr_labels(light_codes,
                "L", cdr_scheme))

            msa, hctr, lctr = [], 0, 0
            for seq, annotation in zip(seqs, annotations):
                vgene, _, _, _, species = vj.assign_vj_genes(annotation, seq,
                        "unknown", "identity")
                if annotation[2] == "H":
                    msa.append( (heavy_msa[hctr], vgene, species, annotation[2]) )
                    hctr += 1
                else:
                    msa.append( (light_msa[lctr], vgene, species, annotation[2]) )
                    lctr += 1

            self.random_search_function_test(msa, heavy_codes,
                        light_codes, local_db, seqinfos,
                        seqs)
            self.distmat_construction_test(msa, heavy_codes,
                        light_codes, local_db)

            del local_db
        shutil.rmtree("TEMP_DB.db")




    def test_local_db_search_setup(self):
        """Check that the local db management tool is set
        up correctly."""
        current_dir = os.path.abspath(os.path.dirname(
            __file__))
        data_filepath = os.path.join(current_dir,
                "test_data", "addtnl_test_data.fasta.gz")

        if not os.path.exists('TEMP_DB.db'):
            build_database_from_fasta([data_filepath],
                "TEMP_DB.db", "TEMP_FILE",
                numbering_scheme="imgt",
                cdr_definition_scheme="imgt",
                sequence_type="single", receptor_type="mab",
                pid_threshold=0.7, user_memo="")

        seqs, seqinfos = [], []

        for seqinfo, seq in read_fasta(data_filepath):
            seqs.append(seq)
            seqinfos.append(seqinfo)

        local_db = LocalDBTool("TEMP_DB.db")

        # Check that the search tool retrieves the correct metadata.
        metadata = local_db.get_database_metadata()
        self.assertTrue('imgt' == metadata[0])
        self.assertTrue('mab' == metadata[1])
        self.assertTrue('single' == metadata[2])
        self.assertTrue('imgt' == metadata[3])
        self.assertTrue('' == metadata[4])

        # Next check that the position counts table matches
        # expected. First check for heavy. The database
        # returns the tables all together but the Python
        # test code will check them separately.
        dtables = local_db.get_database_counts()

        sca = SingleChainAnnotator(scheme="imgt")
        annotations = sca.analyze_seqs(seqs)

        numbering = return_imgt_canonical_numbering_cpp()
        sta = SequenceTemplateAligner(numbering, "H", "imgt",
                "imgt")

        heavy_chains = [(a,s) for (a,s) in zip(annotations,
            seqs) if a[2]=="H"]
        msa = [sta.align_sequence(h[1], h[0][0], False) for
                h in heavy_chains]

        labels = sca.assign_cdr_labels(numbering,
                "H", "imgt")
        gt_dtable = get_kmer_counts(msa, labels)

        for key, value in gt_dtable.items():
            self.assertTrue(dtables[0][key]==value)

        del numbering, msa, labels, sta



        # ...then light.
        numbering = return_imgt_canonical_numbering_cpp()
        sta = SequenceTemplateAligner(numbering, "L", "imgt",
                "imgt")

        light_chains = [(a,s) for (a,s) in zip(annotations,
            seqs) if a[2]!="H"]
        msa = [sta.align_sequence(h[1], h[0][0], False) for
                h in light_chains]

        labels = sca.assign_cdr_labels(numbering,
                "L", "imgt")
        gt_dtable = get_kmer_counts(msa, labels)

        for key, value in gt_dtable.items():
            self.assertTrue(dtables[1][key]==value)

        del numbering, msa, labels, sta

        del local_db
        shutil.rmtree("TEMP_DB.db")




    def random_search_function_test(self, msa, heavy_codes,
        light_codes, db_tool, seqinfo, raw_seqs):
        """Runs exact tests to serve as ground truth for comparison with
        fast-scheme searches."""
        AALIST = list("ACDEFGHIKLMNPQRSTVWY")

        # Randomly select a sequence, randomly mutate it and
        # find its closest matches using a random setting for
        # distance and other parameters. Check it on the fly
        # against the search tool. Simultaneously store the
        # results grouped by search criteria so that we can
        # later test all of the results against the batch
        # search.
        random.seed(123)
        output_search_res = {}

        for ctr in range(400):
            idx = random.randint(0, len(msa) - 1)
            query_seq = list(msa[idx][0])
            if msa[idx][3] == "H":
                codes = heavy_codes
            else:
                codes = light_codes
            mut_percentage = random.uniform(0, 0.2)

            for i, qletter in enumerate(query_seq):
                if qletter == '-':
                    continue
                if codes[1][i].startswith("cdr") and \
                        random.uniform(0,1) < mut_percentage:
                    query_seq[i] = random.choice(AALIST)

            query_seq = "".join(query_seq)
            cutoff_selector = random.randint(0,2)
            if cutoff_selector == 0:
                cutoff = 0.4
            elif cutoff_selector == 1:
                cutoff = 0.33
            else:
                cutoff = 0.25

            if random.randint(0,1) == 1:
                blosum_cutoff = -1
            else:
                blosum_cutoff = 0.5

            if random.randint(0,1) == 1:
                mode = "123"
            else:
                mode = "3"

            length_selector = random.randint(0,2)
            if length_selector == 0:
                max_cdr_length_shift = 2
            elif length_selector == 1:
                max_cdr_length_shift = 0
            else:
                max_cdr_length_shift = 3

            if random.randint(0,4) == 4:
                vgene = msa[idx][1]
                species = msa[idx][2]
            else:
                species, vgene = "", ""

            check_seq_retrieval = random.randint(0,1) == 0

            hits = db_tool.search(query_seq,
                    (codes[0], 1, msa[idx][3], ""),
                    mode, cutoff, blosum_cutoff,
                    max_cdr_length_shift, 1000,
                    vgene, species)
            hits = sorted(hits, key=lambda x: (x[1], x[0]))

            gt_hit_idx = perform_exact_search(query_seq, msa, msa[idx][3],
                    codes, mode, cutoff, blosum_cutoff,
                    max_cdr_length_shift, 1000,
                    vgene, species)

            criteria = (mode, cutoff, max_cdr_length_shift, blosum_cutoff)
            if criteria not in output_search_res:
                output_search_res[criteria] = []
            output_search_res[ criteria ].append(
                    (query_seq, vgene, species, gt_hit_idx,
                        msa[idx][3], codes[0])
                )

            # If using a BLOSUM cutoff, distance will be floating
            # point, check result is close. Otherwise check for
            # exact match.
            if blosum_cutoff > 0:
                self.assertTrue([h[0] for h in hits]==
                                [h[0]for h in gt_hit_idx])
                self.assertTrue(np.allclose([h[1] for h in hits],
                                [h[1] for h in gt_hit_idx]))

            else:
                if hits != gt_hit_idx:
                    print(f"Error! {ctr}, {mut_percentage} mut percentage, cutoff {cutoff}, "
                            f"max hits 100, vgene_filter {vgene}, "
                            f"num gt hits was {len(gt_hit_idx[0])}", flush=True)
                self.assertTrue(hits==gt_hit_idx)

            # Randomly check some of the sequences to make sure the
            # metadata and sequence retrieved from the database for a given
            # id code match those in the input.
            if check_seq_retrieval and len(hits) > 0:
                for hit in hits:
                    db_seq, db_metadata = db_tool.get_sequence(hit[0])
                    self.assertTrue(db_seq==raw_seqs[hit[0]])
                    self.assertTrue(seqinfo[hit[0]]==db_metadata)


        # Now evaluate batch search tools, using batches with specific settings.
        for criteria, seq_data in output_search_res.items():
            query_seqs = [q[0] for q in seq_data]
            annotations = [ (q[5], 1, q[4], "") for q in seq_data]
            vgene_list = [q[1] for q in seq_data]
            species_list = [q[2] for q in seq_data]

            for i, hit_idx in enumerate(db_tool.search_batch(
                    query_seqs, annotations, criteria[0],
                    criteria[1], criteria[3], criteria[2],
                    1000, vgene_list, species_list)):
                hit_idx = sorted(hit_idx, key=lambda x: (x[1], x[0]))
                # If using BLOSUM which is floating point, check
                # allclose. Hamming distance is integer so ground truth
                # should match test exactly.
                if criteria[3] < 0:
                    if hit_idx != seq_data[i][3]:
                        print(f"Batch error! {criteria}, {hit_idx}, {seq_data[i]}",
                            flush=True)
                    self.assertTrue(hit_idx==seq_data[i][3])
                else:
                    hit_idx = [h[1] for h in hit_idx]
                    gt_hits = [s[1] for s in seq_data[i][3]]
                    self.assertTrue(np.allclose(hit_idx, gt_hits))



    def distmat_construction_test(self, msa, heavy_codes,
            light_codes, local_db):
        """Compare the results of exact distance matrix construction
        using a simple if brutally inefficient procedure
        with those provided by the much more efficient
        procedure used by the API."""
        for chain_type in [("H",), ("L", "K")]:
            if chain_type[0] == "H":
                codes = heavy_codes
                ldb_chain_type = "heavy"
            else:
                codes = light_codes
                ldb_chain_type = "light"

            distances, row_idx, col_idx = [], [], []

            for i, seq_data in enumerate(msa):
                if seq_data[3] not in chain_type:
                    continue
                hits = perform_exact_search(seq_data[0], msa, seq_data[3],
                        codes, "3", 0.25, -1, 1, 1000, seq_data[1], seq_data[2])
                for hit in hits:
                    distances.append(hit[1])
                    row_idx.append(i)
                    col_idx.append(hit[0])

            ldb_distances, ldb_rows, ldb_cols = local_db.build_sparse_distance_matrix(
                    ldb_chain_type, "3", 0.25, -1, 10000, 1, "hamming", True, True)

            csr_test = csr_matrix((ldb_distances, (ldb_rows, ldb_cols)),
                    [(len(msa), len(msa))])
            csr_gt = csr_matrix((distances, (row_idx, col_idx)),
                    [(len(msa), len(msa))])
            self.assertTrue(np.allclose(csr_test.toarray(), csr_gt.toarray()))


def perform_exact_search(query, msa, chain_code, msa_codes,
        mode, cdr_cutoff, blosum_cutoff,
        max_cdr_length_shift, max_hits,
        vgene_filter="", species_filter=""):
    """Performs an exact search on an input msa, finding
    all sequences that have hamming distance <= specified.
    Additionally filter for BLOSUM distance (if requested)."""
    retained_dists, hit_idx, cdrlen = [], [], []
    max_hamming = 0

    for j, region in enumerate(["cdr1", "cdr2", "cdr3"]):
        cdrlen.append(len([q for (q,l) in zip(query, msa_codes[1])
            if q != '-' and l == region]))

    max_hamming = [int(round((cdrlen[0] + cdrlen[1]) * cdr_cutoff + 0.0001)),
            int(round(cdrlen[2] * cdr_cutoff + 0.0001))]

    vgene_cutoff = -1
    for i in range(4, len(vgene_filter)):
        if not vgene_filter[i].isnumeric():
            vgene_cutoff = i
            break

    for i, seq_data in enumerate(msa):
        if seq_data[3] != chain_code:
            continue

        cdr_dists = [0,0]
        if vgene_filter != "":
            if seq_data[1][:vgene_cutoff] != vgene_filter[:vgene_cutoff]:
                continue
        if species_filter != "":
            if seq_data[2] != species_filter:
                continue

        for j, region in enumerate(["cdr1", "cdr2", "cdr3"]):
            if str(j+1) not in mode:
                continue
            region_dist = len([a for (a,b,l) in zip(seq_data[0], query,
                msa_codes[1]) if l==region and a != b])
            cdr_dists[int(j/2)] += region_dist

            # Add an arbitrary large number if the cdr length is
            # unacceptable or if an individual distance is
            # unacceptable so that the sequence is not saved.
            if region == "cdr3":
                region_len = len([s for (s,l) in zip(seq_data[0], msa_codes[1])
                    if l==region and s != '-'])
                if region_len > cdrlen[j] + max_cdr_length_shift or \
                        region_len < cdrlen[j] - max_cdr_length_shift:
                    cdr_dists[int(j/2)] += 200

        # If we meet Hamming distance criteria, store sequence
        # as a hit UNLESS BLOSUM cutoff was also specified, in
        # which case calculate BLOSUM distance and see if we
        # also meet THAT cutoff.
        if cdr_dists[0] <= max_hamming[0] and \
                cdr_dists[1] <= max_hamming[1]:
            if blosum_cutoff < 0:
                hit_idx.append(i)
                retained_dists.append(cdr_dists[0] + cdr_dists[1])
            else:
                blosum_dist, nresidues = 0, 0
                allowed_regions = {f"cdr{k}" for k in mode}
                for l1, l2, code in zip(seq_data[0], query,
                                            msa_codes[1]):
                    if code not in allowed_regions:
                        continue
                    l1num = AAMAP[l1] * 22
                    l2num = AAMAP[l2]
                    blosum_dist += BLOSUM[l1num + l2num]
                    if l2 != '-':
                        nresidues += 1
                blosum_dist = float(blosum_dist) / float(nresidues)
                if blosum_dist <= blosum_cutoff:
                    hit_idx.append(i)
                    retained_dists.append(blosum_dist)


    if len(hit_idx) == 0:
        return []

    hit_idx = list(zip(hit_idx, retained_dists))
    hit_idx = sorted(hit_idx, key=lambda x: (x[1], x[0]))

    return hit_idx[:max_hits]


def get_kmer_counts(msa, cdr_labels):
    """Extract kmers from the cdrs of an msa and
    count the number present of each."""
    dimer_count_tables = {}

    cdr_start = cdr_labels.index('cdr3')
    cdr_end = cdr_labels.index('fmwk4')
    ndimers = cdr_end - cdr_start - 1

    for m in msa:
        cdr = m[cdr_start:cdr_end]
        cdrlen = len([a for a in cdr if a != '-'])
        if cdrlen not in dimer_count_tables:
            dimer_count_tables[cdrlen] = [0] * ndimers * 21 * 21
        for i in range(ndimers):
            kmer_code = AAMAP[cdr[i]] * 21 + AAMAP[cdr[i+1]]
            dimer_count_tables[cdrlen][kmer_code+i*21*21] += 1

    return dimer_count_tables



if __name__ == "__main__":
    unittest.main()
