"""Tests basic functionality for the MultiChainAnnotator class.
Notice that this test relies on the correct functioning of
SingleChainAnnotator, so if test_single_chain_annotator does
not pass, this one will not either."""
import os
import random
import gzip
import unittest
from antpack import PairedChainAnnotator, SingleChainAnnotator
from antpack.scoring_tools.scoring_constants import scoring_constants as SCCONST



class TestMultiChainAnnotator(unittest.TestCase):


    def test_error_checking(self):
        """Check that sequences which have known issues are flagged
        as such, and that deliberately invalid inputs are recognized."""
        # Pass dummy sequences with errors.
        aligner = PairedChainAnnotator()

        results = aligner.analyze_seq("YaY")
        self.assertTrue(results[0][3].startswith("Invalid sequence"))
        self.assertTrue(results[1][3].startswith("Invalid sequence"))

        results = aligner.analyze_seq("YBW")
        self.assertTrue(results[0][3].startswith("Invalid sequence"))
        self.assertTrue(results[1][3].startswith("Invalid sequence"))

        results = aligner.analyze_seq("Y K")
        self.assertTrue(results[0][3].startswith("Invalid sequence"))
        self.assertTrue(results[1][3].startswith("Invalid sequence"))

        results = aligner.analyze_seq("Y-K")
        self.assertTrue(results[0][3].startswith("Invalid sequence"))
        self.assertTrue(results[1][3].startswith("Invalid sequence"))



    def test_sequence_extraction(self):
        """Take test heavy and light chains and number them using
        SingleChainAnnotator. Next, join pairs of chains with a random
        number of random amino acids and add a random number of
        random amino acids at the beginning and end. Compare the
        numbering assigned by SingleChainAnnotator to that assigned
        by MultiChainAnnotator (this of course assumes that
        SingleChainAnnotator is functioning correctly)."""
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))
        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs = [line.strip().split(",")[0] for line in fhandle]
        os.chdir(current_dir)

        sc_aligner = SingleChainAnnotator(scheme="imgt")
        m_aligner = PairedChainAnnotator(scheme="imgt")

        alignments = [sc_aligner.analyze_seq(seq) for seq in seqs]

        heavy_chains = [(seq, a) for a, seq in zip(alignments, seqs)
                if a[2] == "H"]
        light_chains = [(seq, a) for a, seq in zip(alignments, seqs)
                if a[2] != "H"]

        random.seed(0)

        for hc, lc in zip(heavy_chains, light_chains):
            prefix = [SCCONST.aa_list[random.randint(0,19)]
                for i in range(random.randint(0,25))]
            suffix = [SCCONST.aa_list[random.randint(0,19)]
                for i in range(random.randint(0,25))]
            joiner1 = [SCCONST.aa_list[random.randint(0,19)]
                for i in range(random.randint(0,25))]
            joiner2 = [SCCONST.aa_list[random.randint(0,19)]
                for i in range(random.randint(0,25))]

            merged_hc = "".join( ["".join(prefix), hc[0], "".join(joiner1) ] )
            merged_lc = "".join( ["".join(joiner2), lc[0], "".join(suffix) ] )
            merged_chain = merged_hc + merged_lc
            if merged_hc.endswith("SSA"):
                continue

            hc_align = sc_aligner.analyze_seq(merged_hc)
            lc_align = sc_aligner.analyze_seq(merged_lc)
            if hc_align[1] < 0.8 or lc_align[1] < 0.8:
                continue

            for i, lcpos in enumerate(lc_align[0]):
                if lcpos != "-":
                    break
            for j, lcpos in reversed(list(enumerate(lc_align[0]))):
                if lcpos != "-":
                    break

            trimmed_lc_seq = merged_lc[i:j+1]
            trimmed_lc_align = lc_align[0][i:j+1]
            if not trimmed_lc_align[0] == "1":
                continue

            for i, hcpos in enumerate(hc_align[0]):
                if hcpos != "-":
                    break
            for j, hcpos in reversed(list(enumerate(hc_align[0]))):
                if hcpos != "-":
                    break

            trimmed_hc_seq = merged_hc[i:j+1]
            trimmed_hc_align = hc_align[0][i:j+1]
            if not trimmed_hc_align[-1] == "128":
                continue

            mc_heavy, mc_light = m_aligner.analyze_seq(merged_chain)
            _, mchn, hstart, hend = m_aligner.trim_alignment(merged_chain, mc_heavy)
            _, mcln, lstart, lend = m_aligner.trim_alignment(merged_chain, mc_light)

            self.assertTrue(mcln == trimmed_lc_align)
            self.assertTrue(mchn == trimmed_hc_align)
            self.assertTrue(merged_chain[lstart:lend] ==
                    trimmed_lc_seq)
            self.assertTrue(merged_chain[hstart:hend] ==
                    trimmed_hc_seq)



if __name__ == "__main__":
    unittest.main()
