"""Test database management & search."""
import os
import random
import unittest
from antpack import (build_database_from_fasta,
        SingleChainAnnotator, LocalDBTool)
from antpack.antpack_cpp_ext import (SequenceTemplateAligner,
        return_imgt_canonical_numbering_cpp)
from antpack.utilities import read_fasta


AAMAP = {k:i for i,k in enumerate("ACDEFGHIKLMNPQRSTVWY-")}


class TestLocalDBManagement(unittest.TestCase):


    def test_local_db_search(self):
        """Check that searches retrieve the sequences that
        we would expect based on a brute-force exact search."""
        current_dir = os.path.abspath(os.path.dirname(
            __file__))
        data_filepath = os.path.join(current_dir,
                "test_data", "addtnl_test_data.fasta.gz")

        build_database_from_fasta(data_filepath,
            "TEMP_DB.db", numbering_scheme="imgt",
            cdr_definition_scheme="imgt",
            sequence_type="single", receptor_type="mab",
            pid_threshold=0.7, user_memo="")

        seqs, seqinfos = [], []

        for seqinfo, seq in read_fasta(data_filepath):
            seqs.append(seq)
            seqinfos.append(seqinfo)

        local_db = LocalDBTool("TEMP_DB.db")

        # First annotate the heavy chains and build an MSA.
        # Use this msa to do brute force searches and compare
        # the results with those we get back from the database.
        sca = SingleChainAnnotator(scheme="imgt")
        annotations = sca.analyze_seqs(seqs)

        chains = [(a,s) for (a,s) in zip(annotations,
            seqs) if a[2]=="H"]
        codes, msa = sca.build_msa([t[1] for t in chains],
                [t[0] for t in chains])

        labels = sca.assign_cdr_labels(codes, "H", "imgt")
        possible_regions = ["cdr", "cdr3", "all"]
        # Randomly select a sequence and find its closest
        # matches using a random setting for distance and
        # region of interest.
        random.seed(123)

        for _ in range(250):
            region = random.choice(possible_regions)
            max_hamming_dist = random.randint(0,10)
            idx = random.randint(0, len(msa) - 1)
            hit_seqs, _, _ = local_db.search_seq(
                    chains[idx][1], chains[idx][0],
                    max_hamming_dist, max_hits=500,
                    region_label=region)
            hit_seqs.sort()

            gt_hit_idx = perform_exact_search(msa, region, idx,
                    labels, max_hamming_dist)

            gt_hits = [chains[idx][1] for idx in gt_hit_idx]
            gt_hits.sort()

            self.assertTrue(len(hit_seqs)==len(gt_hits))
            for hit_seq, gt_hit in zip(hit_seqs, gt_hits):
                self.assertTrue(hit_seq==gt_hit)


        # Now do the same for light chains.
        chains = [(a,s) for (a,s) in zip(annotations,
            seqs) if a[2]!="H"]
        codes, msa = sca.build_msa([t[1] for t in chains],
                [t[0] for t in chains])

        labels = sca.assign_cdr_labels(codes, "L", "imgt")
        possible_regions = ["cdr", "cdr3", "all"]
        # Randomly select a sequence and find its closest
        # matches using a random setting for distance and
        # region of interest.
        random.seed(123)

        for _ in range(250):
            region = random.choice(possible_regions)
            max_hamming_dist = random.randint(0,10)
            idx = random.randint(0, len(msa) - 1)
            hit_seqs, _, _ = local_db.search_seq(
                    chains[idx][1], chains[idx][0],
                    max_hamming_dist, max_hits=1000,
                    region_label=region)
            hit_seqs.sort()

            gt_hit_idx = perform_exact_search(msa, region, idx,
                    labels, max_hamming_dist)

            gt_hits = [chains[idx][1] for idx in gt_hit_idx]
            gt_hits.sort()

            self.assertTrue(len(hit_seqs)==len(gt_hits))
            for hit_seq, gt_hit in zip(hit_seqs, gt_hits):
                self.assertTrue(hit_seq==gt_hit)

        del local_db
        os.remove("TEMP_DB.db")



    def test_local_db_search_setup(self):
        """Check that the local db management tool is set
        up correctly."""
        current_dir = os.path.abspath(os.path.dirname(
            __file__))
        data_filepath = os.path.join(current_dir,
                "test_data", "addtnl_test_data.fasta.gz")

        build_database_from_fasta(data_filepath,
            "TEMP_DB.db", numbering_scheme="imgt",
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
        # expected. First check for heavy...


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
        ptable, _ = local_db.get_database_counts()
        gt_table = get_kmer_counts(msa, numbering, labels)

        for token, counts in gt_table.items():
            self.assertTrue(ptable[token]==counts)

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
        _, ptable = local_db.get_database_counts()
        gt_table = get_kmer_counts(msa, numbering, labels)

        for token, counts in gt_table.items():
            self.assertTrue(ptable[token]==counts)

        del numbering, msa, labels, sta

        del local_db
        os.remove("TEMP_DB.db")





def perform_exact_search(msa, region, idx,
        region_labels, max_hamming_dist):
    """Performs an exact search on an input msa, finding
    all sequences that have hamming distance <= specified."""
    hit_idx = []

    for i, seq in enumerate(msa):
        if region == "all":
            dist = len([a for a,b in zip(seq, msa[idx])
                if a != b])
        else:
            dist = len([a for a,b,code in zip(seq, msa[idx],
                        region_labels) if a != b and
                        code[:len(region)]== region])
        if dist <= max_hamming_dist:
            hit_idx.append(i)

    return hit_idx

def get_kmer_counts(msa, numbering,
        cdr_labels):
    """Extract kmers from the cdrs of an msa and
    count the number present of each."""
    kmer_count_table = {}
    for (cdr, next_fmwk) in [("cdr1", "fmwk2"),
            ("cdr2", "fmwk3"), ("cdr3", "fmwk4")]:
        cdr_start = cdr_labels.index(cdr)
        cdr_end = cdr_labels.index(next_fmwk)
        print(f"{cdr}, {cdr_start}, {cdr_end}")

        for i in range(cdr_start, cdr_end, 2):
            if numbering[i] not in kmer_count_table:
                kmer_count_table[numbering[i]] = [0]*21*21
            if i < cdr_end - 1:
                for m in msa:
                    kmer_code = AAMAP[m[i]] * 21 + AAMAP[m[i+1]]
                    kmer_count_table[numbering[i]][kmer_code] += 1
            else:
                for m in msa:
                    kmer_code = AAMAP[m[i]] * 21 + 20
                    kmer_count_table[numbering[i]][kmer_code] += 1


    return kmer_count_table



if __name__ == "__main__":
    unittest.main()
