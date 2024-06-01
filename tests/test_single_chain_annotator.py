"""Tests basic functionality for the SingleChainAnnotator class."""
import os
import re
import gzip
import copy
import random
import unittest
from antpack import SingleChainAnnotator

class TestSingleChainAnnotator(unittest.TestCase):


    def test_error_checking(self):
        """Check that sequences which have known issues are flagged
        as such, and that deliberately invalid inputs are recognized."""
        # Pass dummy sequences with errors.
        aligner = SingleChainAnnotator(chains=["H", "K", "L"])
        with self.assertRaises(ValueError):
            aligner.analyze_seqs("YYY")

        results = aligner.analyze_seqs(["YaY"])
        self.assertTrue(results[0][3].startswith("Invalid sequence"))
        results = aligner.analyze_seqs(["YBW"])
        self.assertTrue(results[0][3].startswith("Invalid sequence"))
        results = aligner.analyze_seqs(["Y K"])
        self.assertTrue(results[0][3].startswith("Invalid sequence"))
        results = aligner.analyze_seqs(["Y-K"])
        self.assertTrue(results[0][3].startswith("Invalid sequence"))

        results = aligner.analyze_seq("Y-K")
        self.assertTrue(results[3].startswith("Invalid sequence"))
        results = aligner.analyze_seq("yAy")
        self.assertTrue(results[3].startswith("Invalid sequence"))

        with self.assertRaises(RuntimeError):
            sorted_positions = aligner.sort_position_codes(["-", "1"], scheme="imgt")
        with self.assertRaises(RuntimeError):
            sorted_positions = aligner.sort_position_codes(["1", "2", "C3"], scheme="imgt")
        with self.assertRaises(RuntimeError):
            sorted_positions = aligner.sort_position_codes(["1", "2", "3"], scheme="alpha")


    def test_chain_recognition(self):
        """Ensure that single chain annotator can correctly recognize the
        input chain when supplied with something that could be L or H,
        and ensure it can correctly detect sequences with large deletions
        that remove one or more conserved residues."""

        known_K = ("DIVMTQSPSSLTVTAGEKVTMSCKSSQSLLSSGNQKNYLTWYQQIPGQPPKLLIYWASTR"
                    "ESGVPDRFTGSGSGTDFTLTINSVQAEDLAVYYCQNDYTYPLTFGAGTKLELKRTV")
        known_L = ("QSALTQPASVSGSPGQSITISCTGTTSDVGTYNFVSWYQQHPGKAPKAIIFDVTNRPSGI"
                    "SNRFSGSKFGNTASLTISGLQAEDEADYYCAAYTVASTLLFGGGTKVTVLRQP")
        known_H = ("VHLQQSGAELMKPGASVKISCKASGYTFITYWIEWVKQRPGHGLEWIGDILPGSGSTNYN"
                    "ENFKGKATFTADSSSNTAYMQLSSLTSEDSAVYYCARSGYYGNSGFAYWGQGTLVTVSA")

        for scheme in ["martin", "imgt", "kabat"]:
            aligner = SingleChainAnnotator(chains=["H", "K", "L"],
                            scheme = scheme)
            results = aligner.analyze_seqs([known_K, known_L, known_H])
            self.assertTrue(results[0][2] == "K")
            self.assertTrue(results[1][2] == "L")
            self.assertTrue(results[2][2] == "H")

            self.assertTrue(aligner.analyze_seq(known_K)[2] == "K")
            self.assertTrue(aligner.analyze_seq(known_L)[2] == "L")
            self.assertTrue(aligner.analyze_seq(known_H)[2] == "H")

            bad_chain = known_H[:100]
            self.assertTrue(aligner.analyze_seqs([bad_chain])[0][3].startswith("Unexpected"))
            self.assertTrue(aligner.analyze_seq(bad_chain)[3].startswith("Unexpected"))

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

        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs, martin_num, imgt_num, kabat_num = [], [], [], []
            for line in fhandle:
                line_elements = line.strip().split(",")
                seqs.append(line_elements[0])
                martin_num.append(line_elements[1].split("_"))
                imgt_num.append(line_elements[2].split("_"))
                kabat_num.append(line_elements[3].split("_"))

        os.chdir(current_dir)

        numberings = [martin_num, imgt_num, kabat_num]
        schemes = ["martin", "imgt", "kabat"]

        aligners = [SingleChainAnnotator(chains=["H", "K", "L"],
                        scheme=k) for k in schemes]

        for aligner, scheme, numbering in zip(aligners, schemes, numberings):
            total_comparisons, num_correct = compare_results(aligner.analyze_seqs(seqs),
                                    numbering, seqs, scheme)
            print(f"{scheme}: Total comparisons: {total_comparisons}. Num matching: {num_correct}.")
            self.assertTrue(num_correct / total_comparisons > 0.97)


    def test_region_labeling(self):
        """Ensure that the region labels assigned by the region labeling
        procedure correspond to our expectations, using a fairly
        inefficient procedure to determine ground-truth labeling."""
        regex = re.compile("^(?P<numbers>\d*)(?P<letters>\w*)$")

        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs = [line.strip().split(",")[0] for line in fhandle]

        imgt_labels = [(str(i), "fmwk1") for i in range(1,27)] + \
                [(str(i), "cdr1") for i in range(27,40)] + \
                [(str(i), "fmwk2") for i in range(39,56)] + \
                [(str(i), "cdr2") for i in range(56,66)] + \
                [(str(i), "fmwk3") for i in range(66,105)] + \
                [(str(i), "cdr3") for i in range(105,118)] + \
                [(str(i), "fmwk4") for i in range(118,129)]
        imgt_labels = {a:k for (a,k) in imgt_labels}

        martin_heavy = [(str(i), "fmwk1") for i in range(1,26)] + \
                [(str(i), "cdr1") for i in range(26,33)] + \
                [(str(i), "fmwk2") for i in range(33,52)] + \
                [(str(i), "cdr2") for i in range(52,57)] + \
                [(str(i), "fmwk3") for i in range(57,95)] + \
                [(str(i), "cdr3") for i in range(95,103)] + \
                [(str(i), "fmwk4") for i in range(103,114)]
        martin_heavy = {a:k for (a,k) in martin_heavy}
        martin_light = [(str(i), "fmwk1") for i in range(1,26)] + \
                [(str(i), "cdr1") for i in range(26,33)] + \
                [(str(i), "fmwk2") for i in range(33,50)] + \
                [(str(i), "cdr2") for i in range(50,53)] + \
                [(str(i), "fmwk3") for i in range(53,91)] + \
                [(str(i), "cdr3") for i in range(91,97)] + \
                [(str(i), "fmwk4") for i in range(97,108)]
        martin_light = {a:k for (a,k) in martin_light}

        kabat_heavy = [(str(i), "fmwk1") for i in range(1,31)] + \
                [(str(i), "cdr1") for i in range(31,36)] + \
                [(str(i), "fmwk2") for i in range(36,50)] + \
                [(str(i), "cdr2") for i in range(50,66)] + \
                [(str(i), "fmwk3") for i in range(66,95)] + \
                [(str(i), "cdr3") for i in range(95,103)] + \
                [(str(i), "fmwk4") for i in range(103,114)]
        kabat_heavy = {a:k for (a,k) in kabat_heavy}
        kabat_light = [(str(i), "fmwk1") for i in range(1,24)] + \
                [(str(i), "cdr1") for i in range(24,35)] + \
                [(str(i), "fmwk2") for i in range(35,50)] + \
                [(str(i), "cdr2") for i in range(50,57)] + \
                [(str(i), "fmwk3") for i in range(57,89)] + \
                [(str(i), "cdr3") for i in range(89,98)] + \
                [(str(i), "fmwk4") for i in range(98,108)]
        kabat_light = {a:k for (a,k) in kabat_light}



        scheme_labels = {"imgt":{"H":imgt_labels, "K":imgt_labels,
                "L":imgt_labels},
            "martin":{"H":martin_heavy, "K":martin_light,
                "L":martin_light},
            "kabat":{"H":kabat_heavy, "K":kabat_light,
                "L":kabat_light}
            }


        def get_gt_regions(numbering, label_map):
            gt_reg = []
            for n in numbering:
                if n == "-":
                    gt_reg.append("-")
                else:
                    gt_reg.append(label_map[regex.search(n).groups()[0]])
            return gt_reg

        os.chdir(current_dir)

        for scheme in ["imgt", "martin", "kabat"]:
            aligner = SingleChainAnnotator(chains=["H", "K", "L"],
                    scheme=scheme)
            num_err = 0

            for seq in seqs:
                numbering = aligner.analyze_seq(seq, get_region_labels=True)
                gt_regions = get_gt_regions(numbering[0],
                        scheme_labels[scheme][numbering[2]])
                if gt_regions != numbering[-1]:
                    num_err += 1
            self.assertTrue(num_err == 0)


    def test_position_code_sorting(self):
        """Checks the position code sorting function to make
        sure it is sorting positions correctly for different schemes."""
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs = [line.strip().split(",")[0] for line in fhandle]

        os.chdir(current_dir)
        random.seed(123)

        num_err = 0

        for scheme in ["martin", "imgt", "kabat"]:
            aligner = SingleChainAnnotator(chains=["H", "K", "L"], scheme=scheme)
            for seq in seqs:
                numbering = aligner.analyze_seq(seq)[0]
                numbering = [n for n in numbering if n != "-"]
                shuffled_numbering = copy.deepcopy(numbering)
                random.shuffle(shuffled_numbering)
                sorted_numbering = aligner.sort_position_codes(shuffled_numbering,
                        scheme)
                if numbering != sorted_numbering:
                    num_err += 1

        self.assertTrue(num_err == 0)





def compare_results(results, comparator_numbering, seqs, scheme):
    """Compares the numbering generated by AntPack with a comparator,
    and returns the number correct vs total comparisons. Also contains
    some (commented-out) code for writing non-matching results to
    a temporary file for closer inspection."""
    total_comparisons, num_correct = 0, 0
    #outhandle = open(f"temp_{scheme}.csv", "w+", encoding="utf-8")
    for result, comparator, seq in zip(results, comparator_numbering, seqs):
        if result[3] != '':
            continue
        total_comparisons += 1
        if result[0] == comparator:
            num_correct += 1
        #This code writes results to file in a format which is easy to look
        #at. Comment it out normally.
        #else:
        #    for i, resnum in enumerate(result[0]):
        #        if resnum != comparator[i]:
        #            result[0][i] = result[0][i] + "!!"
        #    outhandle.write(f"Sequence,{','.join(list(seq))}\n")
        #    outhandle.write(f"Antpack_result,{','.join(result[0])}\n")
        #    outhandle.write(f"Comparator_result,{','.join(comparator)}\n\n")

    #outhandle.close()
    return total_comparisons, num_correct


if __name__ == "__main__":
    unittest.main()
