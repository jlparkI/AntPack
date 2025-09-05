"""Test database management & search."""
import os
import unittest
from antpack import (build_database_from_fasta,
        SingleChainAnnotator, LocalDBTool)
from antpack.utilities import read_fasta



class TestLocalDBManagement(unittest.TestCase):



    def test_local_db_search(self):
        """Check that searches yield expected results."""
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

        os.remove("TEMP_DB.db")
        os.remove("TEMP_DB.db-shm")
        os.remove("TEMP_DB.db-wal")



if __name__ == "__main__":
    unittest.main()
