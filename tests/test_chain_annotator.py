"""Tests basic functionality for the ChainAnnotator class.
Notice that this test relies on the correct functioning of
SingleChainAnnotator, so if test_single_chain_annotator does
not pass, this one will not either."""
import os
import random
import gzip
import unittest
from Bio import SeqIO
from antpack import ChainAnnotator, SingleChainAnnotator
from antpack.scoring_tools.scoring_constants import scoring_constants as SCCONST



class TestChainAnnotator(unittest.TestCase):


    def test_error_checking(self):
        """Check that sequences which have known issues are flagged
        as such, and that deliberately invalid inputs are recognized."""
        # Pass dummy sequences with errors.
        aligner = ChainAnnotator()

        results = aligner.analyze_seq("YaY", 0)
        self.assertTrue(results[0][3].startswith("Invalid sequence"))

        results = aligner.analyze_seq("YBW", 0)
        self.assertTrue(results[0][3].startswith("Invalid sequence"))

        results = aligner.analyze_seq("Y K", 0)
        self.assertTrue(results[0][3].startswith("Invalid sequence"))

        results = aligner.analyze_seq("Y-K", 0)
        self.assertTrue(results[0][3].startswith("Invalid sequence"))



    def test_sequence_extraction(self):
        """Take a set of test chains and number them using
        either SingleChainAnnotator or ChainAnnotator; results
        should be the same.

        Next, join pairs of these sequences and ensure that
        ChainAnnotator can still extract the correct subregions."""
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))
        with gzip.open("addtnl_test_data.fasta.gz", "rt") as fhandle:
            seqs = [str(seqrec.seq) for seqrec in
                    SeqIO.parse(fhandle, "fasta")]

        os.chdir(current_dir)

        sc_aligner = SingleChainAnnotator(scheme="imgt")
        chain_aligner = ChainAnnotator(scheme="imgt")
        init_alignments = []

        # First, an easy test -- make sure that ChainAnnotator
        # and SingleChainAnnotator give the same result on
        # sequences where only one chain is expected.
        for seq in seqs:
            anal_sc = sc_aligner.analyze_seq(seq)
            init_alignments.append(anal_sc)

            anal_chain = chain_aligner.analyze_seq(seq, 0.7)
            self.assertTrue(len(anal_chain) == 1)
            self.assertTrue(len(anal_chain[0][0]) == len(seq))

            self.assertTrue(anal_chain[0][2] == anal_sc[2])
            if anal_sc[-1] != "":
                continue
            self.assertTrue(anal_chain[0][0] == anal_sc[0])

        nbad1 = 0
        nbad2 = 0
        nbad3 = 0

        for i in range(0, len(seqs), 2):
            joined_seq = seqs[i] + seqs[i+1]
            annotations = chain_aligner.analyze_seq(joined_seq, 0.7)
            if init_alignments[i][-1] != "" or init_alignments[i+1][-1] != "":
                continue

            self.assertTrue(len(annotations)==2)
            self.assertTrue(len(annotations[0][0]) == len(joined_seq))
            self.assertTrue(len(annotations[1][0]) == len(joined_seq))

            if '1' not in annotations[0][0] or '1' not in annotations[1][0]:
                import pdb
                pdb.set_trace()
                nbad3 += 1
                continue

            if annotations[0][0].index('1') < annotations[1][0].index('1'):
                a1, a2 = annotations[0], annotations[1]
            else:
                a1, a2 = annotations[1], annotations[0]

            trs1, tra1, _, _ = chain_aligner.trim_alignment(joined_seq, a1)
            trs2, tra2, _, _ = chain_aligner.trim_alignment(joined_seq, a2)
            og_seq1, og_al1, _, _ = chain_aligner.trim_alignment(seqs[i],
                        init_alignments[i])
            og_seq2, og_al2, _, _ = chain_aligner.trim_alignment(seqs[i+1],
                        init_alignments[i+1])

            if (tra1 != og_al1):
                nbad2 += 1
            if (tra2 != og_al2):
                nbad2 += 1
            #self.assertTrue(tra1 == og_al1)
            #self.assertTrue(tra2 == og_al2)

        print("***********")
        print(nbad3)
        print(nbad2)
        print("***************")


if __name__ == "__main__":
    unittest.main()
