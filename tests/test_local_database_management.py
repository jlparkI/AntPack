"""Test database management & search."""
import os
import random
import unittest
from antpack import (build_database_from_fasta,
        SingleChainAnnotator, LocalDBTool)
from antpack.utilities import read_fasta



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
            cdr_definition_scheme="north",
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

        labels = sca.assign_cdr_labels(codes, "H", "north")
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

        labels = sca.assign_cdr_labels(codes, "L", "north")
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
        AAMAP = {k:i for i,k in enumerate("ACDEFGHIKLMNPQRSTVWY-")}

        current_dir = os.path.abspath(os.path.dirname(
            __file__))
        data_filepath = os.path.join(current_dir,
                "test_data", "addtnl_test_data.fasta.gz")

        build_database_from_fasta(data_filepath,
            "TEMP_DB.db", numbering_scheme="imgt",
            cdr_definition_scheme="north",
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
        self.assertTrue('north' == metadata[3])
        self.assertTrue('' == metadata[4])

        # Next check that the position counts table matches
        # expected. First check for heavy...


        sca = SingleChainAnnotator(scheme="imgt")
        annotations = sca.analyze_seqs(seqs)

        heavy_chains = [(a,s) for (a,s) in zip(annotations,
            seqs) if a[2]=="H"]
        heavy_codes, msa = sca.build_msa(
                [t[1] for t in heavy_chains],
                [t[0] for t in heavy_chains])

        labels = sca.assign_cdr_labels(heavy_codes,
                "H", "north")
        cdr_labels = [(k,t) for (k,(l,t)) in enumerate(
            zip(labels, heavy_codes)) if l.startswith("cdr")]
        ptable, _ = local_db.get_database_counts()

        for i, (k,label) in enumerate(cdr_labels):
            gt_counts = [0]*21
            for m in msa:
                gt_counts[AAMAP[m[k]]] += 1
            self.assertTrue(ptable[label] == gt_counts)

        del heavy_codes, msa, labels


        # ...then light.
        light_chains = [(a,s) for (a,s) in zip(annotations,
            seqs) if a[2]!="H"]
        light_codes, msa = sca.build_msa(
                [t[1] for t in light_chains],
                [t[0] for t in light_chains])

        labels = sca.assign_cdr_labels(light_codes,
                "L", "north")
        cdr_labels = [(k,t) for (k,(l,t)) in enumerate(
            zip(labels, light_codes)) if l.startswith("cdr")]
        _, ptable = local_db.get_database_counts()

        for i, (k,label) in enumerate(cdr_labels):
            gt_counts = [0]*21
            for m in msa:
                gt_counts[AAMAP[m[k]]] += 1
            self.assertTrue(ptable[label] == gt_counts)

        del light_codes, msa, labels

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




if __name__ == "__main__":
    unittest.main()
