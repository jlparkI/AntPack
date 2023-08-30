"""Tests basic functionality for the SingleChainAnnotator class."""
import os
import unittest
from antpack import SingleChainAnnotator

class TestSingleChainAnnotator(unittest.TestCase):


    def test_error_checking(self):
        """Check that sequences which have known issues are flagged
        as such, and that deliberately invalid inputs are recognized."""

        # Pass inappropriate settings and make sure an error is raised.
        with self.assertRaises(ValueError):
            imgt = SingleChainAnnotator(species = "dinosaur")
        with self.assertRaises(ValueError):
            imgt = SingleChainAnnotator(species = ["woolly yak"])
        with self.assertRaises(ValueError):
            imgt = SingleChainAnnotator(species = ["human"], chains="VK")
        with self.assertRaises(ValueError):
            imgt = SingleChainAnnotator(species = ["human"], chains=["VK"])
        with self.assertRaises(ValueError):
            imgt = SingleChainAnnotator(species = ["rat"], chains=["A"])


        # Pass dummy sequences with errors.
        imgt = SingleChainAnnotator(species = ["all"], chains=["H", "K", "L"])
        with self.assertRaises(ValueError):
            imgt.analyze_online_seqs("YYY")

        results = imgt.analyze_online_seqs(["YaY"])
        self.assertTrue(results[0][3].startswith("Invalid sequence"))
        results = imgt.analyze_online_seqs(["YBW"])
        self.assertTrue(results[0][3].startswith("Invalid sequence"))
        results = imgt.analyze_online_seqs(["Y K"])
        self.assertTrue(results[0][3].startswith("Invalid sequence"))
        results = imgt.analyze_online_seqs(["Y-K"])
        self.assertTrue(results[0][3].startswith("Invalid sequence"))

        results = imgt.analyze_seq("Y-K")
        self.assertTrue(results[3].startswith("Invalid sequence"))
        results = imgt.analyze_seq("yAy")
        self.assertTrue(results[3].startswith("Invalid sequence"))


    def test_chain_recognition(self):
        """Ensure that single chain annotator can correctly recognize the
        input chain when supplied with something that could be L or H,
        and ensure it can correctly detect sequences with large deletions
        that remove one or more conserved residues."""
        imgt = SingleChainAnnotator(species = ["all"], chains=["H", "K", "L"])

        known_K = ("DIVMTQSPSSLTVTAGEKVTMSCKSSQSLLSSGNQKNYLTWYQQIPGQPPKLLIYWASTR"
                    "ESGVPDRFTGSGSGTDFTLTINSVQAEDLAVYYCQNDYTYPLTFGAGTKLELKRTV")
        known_L = ("QSALTQPASVSGSPGQSITISCTGTTSDVGTYNFVSWYQQHPGKAPKAIIFDVTNRPSGI"
                    "SNRFSGSKFGNTASLTISGLQAEDEADYYCAAYTVASTLLFGGGTKVTVLRQP")
        known_H = ("VHLQQSGAELMKPGASVKISCKASGYTFITYWIEWVKQRPGHGLEWIGDILPGSGSTNYN"
                    "ENFKGKATFTADSSSNTAYMQLSSLTSEDSAVYYCARSGYYGNSGFAYWGQGTLVTVSA")

        results = imgt.analyze_online_seqs([known_K, known_L, known_H])
        self.assertTrue(results[0][2] == "K")
        self.assertTrue(results[1][2] == "L")
        self.assertTrue(results[2][2] == "H")

        r1 = imgt.analyze_seq(known_K)
        r2 = imgt.analyze_seq(known_L)
        r3 = imgt.analyze_seq(known_H)
        self.assertTrue(r1[2] == "K")
        self.assertTrue(r2[2] == "L")
        self.assertTrue(r3[2] == "H")

        bad_chain = known_H[:100]
        results = imgt.analyze_online_seqs([bad_chain])
        self.assertTrue(results[0][3].startswith("Unexpected"))
        results = imgt.analyze_seq(bad_chain)
        self.assertTrue(results[3].startswith("Unexpected"))


    def test_performance(self):
        """Run a batch of test data (approximately 1600 sequences from the
        PDB) to ensure that numbering is consistent with numbering generated
        by another tool. There will occasionally be small differences in
        cases where there are multiple possible acceptable alignments,
        but in general we expect the numbering to be the same the vast
        majority of the time."""
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))
        with open("test_data.csv", "r") as fhandle:
            _ = fhandle.readline()
            seqs, comparator_numbering = [], []
            for line in fhandle:
                line_elements = line.strip().split(",")
                seqs.append(line_elements[0])
                comparator_numbering.append(line_elements[1].split("_"))

        imgt = SingleChainAnnotator(species = ["all"], chains=["H", "K", "L"])
        results = imgt.analyze_online_seqs(seqs)
        total_comparisons, num_correct = 0, 0
        for result, comparator in zip(results, comparator_numbering):
            if result[3] != '':
                continue
            total_comparisons += 1
            if result[0] == comparator:
                num_correct += 1
        print(f"Total comparisons: {total_comparisons}. Num matching: {num_correct}.")
        self.assertTrue(num_correct / total_comparisons > 0.97)


if __name__ == "__main__":
    unittest.main()
