"""Tests basic functionality for the MultiChainAnnotator class.
Notice that this test relies on the correct functioning of
SingleChainAnnotator, so if test_single_chain_annotator does
not pass, this one will not either."""
import os
import random
import gzip
import unittest
from antpack import MultiChainAnnotator, SingleChainAnnotator
from antpack.scoring_tools.scoring_constants import scoring_constants as SCCONST


# A list of AAs which are unlikely to disrupt the alignment
# that can be used as filler.
JUNK_AAS = ["A", "N", "L", "P", "G"]



class TestMultiChainAnnotator(unittest.TestCase):


    def test_error_checking(self):
        """Check that sequences which have known issues are flagged
        as such, and that deliberately invalid inputs are recognized."""
        # Pass dummy sequences with errors.
        aligner = MultiChainAnnotator()

        results = aligner.analyze_seq("YaY")
        self.assertTrue(results[0][3].startswith("Invalid sequence"))
        results = aligner.analyze_seq("YBW")
        self.assertTrue(results[0][3].startswith("Invalid sequence"))
        results = aligner.analyze_seq("Y K")
        self.assertTrue(results[0][3].startswith("Invalid sequence"))
        results = aligner.analyze_seq("Y-K")
        self.assertTrue(results[0][3].startswith("Invalid sequence"))



    def test_sequence_extraction(self):
        """Take test heavy and light chains and number them using
        SingleChainAnnotator. Next, join pairs of chains with a random
        number of random amino acids and add a random number of
        random amino acids at the beginning and end. The heavy
        and light chains extracted by MultiChainAnnotator should
        match the ones that we put in, and the numbering
        assigned to them should match whatever SingleChainAnnotator
        assigned."""
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))
        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs = [line.strip().split(",")[0] for line in fhandle]
        os.chdir(current_dir)

        sc_aligner = SingleChainAnnotator(scheme="imgt")
        m_aligner = MultiChainAnnotator(scheme="imgt")

        alignments = [sc_aligner.analyze_seq(seq) for seq in seqs]

        heavy_chains = [(seq, a) for a, seq in zip(alignments, seqs)
                if a[2] == "H"]
        light_chains = [(seq, a) for a, seq in zip(alignments, seqs)
                if a[2] != "H"]

        random.seed(123)

        for hc, lc in zip(heavy_chains, light_chains):
            prefix = [SCCONST.aa_list[random.randint(0,19)]
                for i in range(random.randint(0,25))]
            suffix = [SCCONST.aa_list[random.randint(0,19)]
                for i in range(random.randint(0,25))]
            joiner = [SCCONST.aa_list[random.randint(0,19)]
                for i in range(random.randint(0,25))]
            merged_chain = "".join( ["".join(prefix),
                hc[0], "".join(joiner), lc[0], "".join(suffix)] )
            mc_analysis = m_aligner.analyze_seq(merged_chain)
            import pdb
            pdb.set_trace()

        import pdb
        pdb.set_trace()


if __name__ == "__main__":
    unittest.main()
