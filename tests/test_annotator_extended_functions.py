"""Tests extended functionality shared by the SingleChainAnnotator and
PairedChainAnnotator classes, such as building distance matrices,
building and slicing MSAs, sorting position codes etc. We use
SingleChainAnnotator here since all functions shown are shared
by SingleChain and PairedChain annotators so only one set
of tests is needed."""
import os
import re
import random
import gzip
import copy
import unittest
import numpy as np
from antpack import SingleChainAnnotator

class TestAnnotatorExtendedFunctions(unittest.TestCase):


    def test_error_checking(self):
        """Check that sequences which have known issues are flagged
        as such, and that deliberately invalid inputs are recognized."""
        # Pass dummy sequences with errors.
        aligner = SingleChainAnnotator(chains=["H", "K", "L"])

        with self.assertRaises(RuntimeError):
            _ = aligner.sort_position_codes(["a", "1"])
        with self.assertRaises(RuntimeError):
            _ = aligner.sort_position_codes(["1", "2", "C3"])

        with self.assertRaises(RuntimeError):
            _ = aligner.build_msa(["AAAAAA", "AATTAAA"],
                    [(["1", "2", "3", "4", "5", "6"], 0.8, "T", ""),
                    (["1", "2", "3", "4", "5", "6", "7"],
                        0.8, "H", "")])
        with self.assertRaises(RuntimeError):
            _ = aligner.build_msa(["AAAAAA", "AATTAAA"],
                    [(["1", "2", "3", "4", "5", "6", "7"],
                        0.8, "H", "")])

        with self.assertRaises(RuntimeError):
            _ = aligner.trim_alignment("AAAAAA",
                (["1", "2", "3", "4", "5"], 0.2, "H", "") )

        with self.assertRaises(RuntimeError):
            _ = aligner.slice_msa(["AA--AA", "AATTAAA"],
                    ["1", "2", "3", "4", "5", "6"],
                    "H", "all")
        with self.assertRaises(RuntimeError):
            _ = aligner.slice_msa(["AA--AA", "AATTAA"],
                    ["1", "2", "3", "4", "5", "6"],
                    "H", "dummy")
        with self.assertRaises(RuntimeError):
            _ = aligner.slice_msa(["AA--AA", "AATTAA"],
                    ["1", "2", "3", "4", "5", "6"],
                    "H", "fmwk1", "hobgoblin")

        with self.assertRaises(RuntimeError):
            _ = aligner.build_distance_matrix(np.empty((2,2)),
                    ["AA--AA", "AATTAAA"],
                    ["1", "2", "3", "4", "5", "6"],
                    "H", "all")
        with self.assertRaises(RuntimeError):
            _ = aligner.build_distance_matrix(np.empty((2,2)),
                    ["AA--AA", "AATTAA"],
                    ["1", "2", "3", "4", "5", "6"],
                    "H", "UKNOWN")
        with self.assertRaises(RuntimeError):
            _ = aligner.build_distance_matrix(np.empty((2,2)),
                    ["AA--AA", "AATTAA"],
                    ["1", "2", "3", "4", "5", "6"],
                    "H", "fmwk1", "hobgoblin")
        with self.assertRaises(RuntimeError):
            _ = aligner.build_distance_matrix(np.empty((3,3)),
                    ["AA--AA", "AATTAA"],
                    ["1", "2", "3", "4", "5", "6"],
                    "H", "fmwk1")


    def test_alignment_trimming(self):
        """Make sure the alignment trimming procedure yields correct results.
        Since this procedure is the same for mabs and TCRs, we do not need
        to check separately."""
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs = [line.strip().split(",")[0] for line in fhandle]

        padded_seqs = []
        AAs = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N",
            "P", "Q", "R", "S", "T", "V", "W", "Y"]

        random.seed(123)

        for seq in seqs:
            left_padding = random.randint(0,10)
            right_padding = random.randint(0,10)
            padded_seq = "".join([random.choice(AAs) for i in range(left_padding)]) + \
                    seq + "".join([random.choice(AAs) for i in range(right_padding)])
            padded_seqs.append(padded_seq)

        os.chdir(current_dir)
        aligner = SingleChainAnnotator(chains=["H", "K", "L"],
                    scheme="imgt")
        aligner_results = aligner.analyze_seqs(padded_seqs)

        for sequence, alignment in zip(padded_seqs, aligner_results):
            numbering = alignment[0]
            exstart = next((i for i in range(len(numbering)) if numbering[i] != '-'), 0)
            exend = next((i for i in range(exstart + 1, len(numbering)) if numbering[i] == '-'),
                    len(numbering))
            trimmed_sequence = sequence[exstart:exend]
            trimmed_numbering = numbering[exstart:exend]

            eval_seq, eval_numbering, eval_start, eval_end = aligner.trim_alignment(sequence,
                    alignment)

            self.assertTrue(eval_seq == trimmed_sequence)
            self.assertTrue(eval_seq == sequence[eval_start:eval_end])
            self.assertTrue(eval_numbering == trimmed_numbering)



    def test_distance_matrix_construction(self):
        """Make sure the msa slicing procedure and distance
        matrix construction yield correct results."""
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs = [line.strip().split(",")[0] for line in fhandle]

        os.chdir(current_dir)
        aligner = SingleChainAnnotator(chains=["H", "K", "L"],
                    scheme="imgt")
        aligner_results = aligner.analyze_seqs(seqs)
        chain_dict = {k:{"chains":[], "annotations":[],
            "msa":None, "cdr_labels":None,
            "poscodes":None} for
            k in ["H", "L"]}

        for seq, result in zip(seqs, aligner_results):
            if result[2] == "H":
                chain_dict["H"]["chains"].append(seq)
                chain_dict["H"]["annotations"].append(result)
            else:
                chain_dict["L"]["chains"].append(seq)
                chain_dict["L"]["annotations"].append(result)


        # Since we test region labeling below, for now in
        # this test we will assume it is working correctly!
        for chain_type in ["H", "L"]:
            poscodes, msa = aligner.build_msa(
                    chain_dict[chain_type]["chains"],
                    chain_dict[chain_type]["annotations"])
            chain_dict[chain_type]["msa"] = msa
            chain_dict[chain_type]["poscodes"] = poscodes
            chain_dict[chain_type]["cdr_labels"] = \
                    aligner.assign_cdr_labels(
                            poscodes, chain_type)

        for region in ["fmwk", "cdr", "cdr3", "fmwk2", "all"]:
            for chain_type in ["H", "L"]:
                seqs = chain_dict[chain_type]["msa"]
                cdr_labels = chain_dict[chain_type]["cdr_labels"]
                poscodes = chain_dict[chain_type]["poscodes"]

                if region == "all":
                    python_slice = seqs
                else:
                    python_slice = [
                        "".join([l for l, c in zip(
                            seq, cdr_labels) if c[:len(region)] ==
                                region])
                            for seq in seqs]

                cpp_slice = aligner.slice_msa(
                        seqs, poscodes, chain_type, region)
                self.assertTrue(python_slice == cpp_slice)

                py_dist_mat = np.empty((len(seqs), len(seqs)))
                cpp_dist_mat = py_dist_mat.copy()
                # Compare the cpp distance matrix to a brutally
                # inefficient but simple and easy to troubleshoot
                # Python version.
                for i, s1 in enumerate(python_slice):
                    py_dist_mat[i,i] = 0
                    for j, s2 in enumerate(python_slice):
                        if j <= i:
                            continue
                        py_dist_mat[i,j] = len(
                                [a for a,b in zip(s1, s2)
                                    if a != b])
                        py_dist_mat[j,i] = py_dist_mat[i,j]

                aligner.build_distance_matrix(cpp_dist_mat,
                        seqs, poscodes, chain_type, region)
                self.assertTrue(np.allclose(cpp_dist_mat,
                    py_dist_mat))
                self.assertTrue(np.allclose(cpp_dist_mat.T,
                    cpp_dist_mat))

        # As a last somewhat silly check, make sure that
        # a distance matrix of a sequence vs itself is 0.
        cpp_dist_mat = np.empty((1,1))
        aligner.build_distance_matrix(cpp_dist_mat,
                seqs[:1], poscodes, chain_type, region)
        self.assertTrue(np.allclose(cpp_dist_mat,
            np.zeros((1,1))))





    def test_region_labeling(self):
        """Ensure that the region labels assigned by the region labeling
        procedure correspond to our expectations, using a fairly
        inefficient procedure to determine ground-truth labeling."""
        regex = re.compile(r"^(?P<numbers>\d*)(?P<letters>\w*)$")

        project_path = os.path.abspath(os.path.dirname(__file__))
        os.chdir(os.path.join(project_path, "test_data"))

        for receptor, testfile in [("mab", "test_data.csv.gz"),
                ("tcr", "tcr_test_data.csv.gz")]:
            with gzip.open(testfile, "rt") as fhandle:
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



            scheme_labels = {"imgt":{k:imgt_labels for k in
                ["H", "K", "L", "A", "B", "D", "G"]},
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

            for scheme in ["imgt", "martin", "kabat"]:
                if receptor == "tcr":
                    if scheme != "imgt":
                        continue
                    aligner = SingleChainAnnotator(
                        chains=["A", "B", "D", "G"],
                        scheme=scheme)
                else:
                    aligner = SingleChainAnnotator(chains=["H", "K", "L"],
                        scheme=scheme)
                num_err = 0

                for seq in seqs:
                    numbering = aligner.analyze_seq(seq)
                    labels = aligner.assign_cdr_labels(numbering[0],
                            numbering[2])

                    gt_regions = get_gt_regions(numbering[0],
                            scheme_labels[scheme][numbering[2]])
                    if gt_regions != labels:
                        num_err += 1
                self.assertTrue(num_err == 0)



    def test_cross_scheme_cdr_assignment(self):
        """For antibodies, there is the possibility of numbering
        using one scheme then assigning cdr labels using another.
        This module tests whether this procedure works as
        expected, by comparing it to a slow but easily
        debugged ground-truth assignment. Since there are
        many such possible conversions, we test the most
        popular rather than exhaustively testing all of
        them."""
        regex = re.compile(r"^(?P<numbers>\d*)(?P<letters>\w*)$")

        project_path = os.path.abspath(os.path.dirname(__file__))
        os.chdir(os.path.join(project_path, "test_data"))

        # In the following, the key is
        # numbering-scheme_chain_cdr-scheme
        scheme_labels = {
                "kabat_H_imgt":dict(
                [(str(i), "fmwk1") for i in range(1,26)] + \
                [(str(i), "cdr1") for i in range(26,36)] + \
                [(str(i), "fmwk2") for i in range(36,51)] + \
                [(str(i), "cdr2") for i in range(51,58)] + \
                [(str(i), "fmwk3") for i in range(58,93)] + \
                [(str(i), "cdr3") for i in range(93,103)] + \
                [(str(i), "fmwk4") for i in range(103,114)]
                ),
                "kabat_L_imgt":dict(
                [(str(i), "fmwk1") for i in range(1,27)] + \
                [(str(i), "cdr1") for i in range(27,33)] + \
                [(str(i), "fmwk2") for i in range(33,50)] + \
                [(str(i), "cdr2") for i in range(50,53)] + \
                [(str(i), "fmwk3") for i in range(53,89)] + \
                [(str(i), "cdr3") for i in range(89,98)] + \
                [(str(i), "fmwk4") for i in range(98,108)]
                ),
                "imgt_H_kabat":dict(
                [(str(i), "fmwk1") for i in range(1,32)] + \
                [(str(i), "cdr1") for i in range(32,41)] + \
                [(str(i), "fmwk2") for i in range(41,55)] + \
                [(str(i), "cdr2") for i in range(55,75)] + \
                [(str(i), "fmwk3") for i in range(75,107)] + \
                [(str(i), "cdr3") for i in range(107,118)] + \
                [(str(i), "fmwk4") for i in range(118,129)]
                ),
                "imgt_L_kabat":dict(
                [(str(i), "fmwk1") for i in range(1,24)] + \
                [(str(i), "cdr1") for i in range(24,41)] + \
                [(str(i), "fmwk2") for i in range(41,56)] + \
                [(str(i), "cdr2") for i in range(56,70)] + \
                [(str(i), "fmwk3") for i in range(70,105)] + \
                [(str(i), "cdr3") for i in range(105,118)] + \
                [(str(i), "fmwk4") for i in range(118,129)]
                ),
                "imgt_H_north":dict(
                [(str(i), "fmwk1") for i in range(1,24)] + \
                [(str(i), "cdr1") for i in range(24,41)] + \
                [(str(i), "fmwk2") for i in range(41,55)] + \
                [(str(i), "cdr2") for i in range(55,67)] + \
                [(str(i), "fmwk3") for i in range(67,105)] + \
                [(str(i), "cdr3") for i in range(105,118)] + \
                [(str(i), "fmwk4") for i in range(118,129)]
                ),
                "imgt_L_north":dict(
                [(str(i), "fmwk1") for i in range(1,24)] + \
                [(str(i), "cdr1") for i in range(24,41)] + \
                [(str(i), "fmwk2") for i in range(41,55)] + \
                [(str(i), "cdr2") for i in range(55,70)] + \
                [(str(i), "fmwk3") for i in range(70,105)] + \
                [(str(i), "cdr3") for i in range(105,118)] + \
                [(str(i), "fmwk4") for i in range(118,129)]
                ),
                "kabat_H_north":dict(
                [(str(i), "fmwk1") for i in range(1,23)] + \
                [(str(i), "cdr1") for i in range(23,36)] + \
                [(str(i), "fmwk2") for i in range(36,50)] + \
                [(str(i), "cdr2") for i in range(50,59)] + \
                [(str(i), "fmwk3") for i in range(59,93)] + \
                [(str(i), "cdr3") for i in range(93,103)] + \
                [(str(i), "fmwk4") for i in range(103,114)]
                ),
                "kabat_L_north":dict(
                [(str(i), "fmwk1") for i in range(1,24)] + \
                [(str(i), "cdr1") for i in range(24,35)] + \
                [(str(i), "fmwk2") for i in range(35,49)] + \
                [(str(i), "cdr2") for i in range(49,57)] + \
                [(str(i), "fmwk3") for i in range(57,89)] + \
                [(str(i), "cdr3") for i in range(89,98)] + \
                [(str(i), "fmwk4") for i in range(98,108)]
                ),
                "imgt_H_aho":dict(
                [(str(i), "fmwk1") for i in range(1,25)] + \
                [(str(i), "cdr1") for i in range(25,39)] + \
                [(str(i), "fmwk2") for i in range(39,56)] + \
                [(str(i), "cdr2") for i in range(56,76)] + \
                [(str(i), "fmwk3") for i in range(76,107)] + \
                [(str(i), "cdr3") for i in range(107,117)] + \
                [(str(i), "fmwk4") for i in range(117,129)]
                ),
                "aho_H_imgt":dict(
                [(str(i), "fmwk1") for i in range(1,27)] + \
                [(str(i), "cdr1") for i in range(27,41)] + \
                [(str(i), "fmwk2") for i in range(41,58)] + \
                [(str(i), "cdr2") for i in range(58,68)] + \
                [(str(i), "fmwk3") for i in range(68,107)] + \
                [(str(i), "cdr3") for i in range(107,139)] + \
                [(str(i), "fmwk4") for i in range(139,150)]
                ),
            }

        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs = [line.strip().split(",")[0] for line in fhandle]

            def get_gt_regions(numbering, label_map):
                gt_reg = []
                for n in numbering:
                    if n == "-":
                        gt_reg.append("-")
                    else:
                        gt_reg.append(label_map[regex.search(n).groups()[0]])
                return gt_reg

            for scheme in ["imgt", "kabat", "aho", "martin"]:
                for cdr_scheme in ["aho", "imgt", "north",
                        "martin", "kabat"]:
                    if scheme == cdr_scheme:
                        continue
                    if f"{scheme}_H_{cdr_scheme}" not in scheme_labels:
                        continue

                    aligner = SingleChainAnnotator(chains=["H", "K", "L"],
                            scheme=scheme)
                    num_err = 0

                    for seq in seqs:
                        numbering = aligner.analyze_seq(seq)
                        labels = aligner.assign_cdr_labels(numbering[0],
                            numbering[2], scheme=cdr_scheme)

                        if scheme in ("aho", "imgt") and \
                                cdr_scheme in ("aho", "imgt"):
                            label_key = f"{scheme}_H_{cdr_scheme}"
                        else:
                            if numbering[2] in ("K", "L"):
                                label_key = f"{scheme}_L_{cdr_scheme}"
                            else:
                                label_key = f"{scheme}_H_{cdr_scheme}"

                        gt_regions = get_gt_regions(numbering[0],
                                scheme_labels[label_key])
                        if gt_regions != labels:
                            num_err += 1
                    self.assertTrue(num_err == 0)



    def test_position_code_sorting(self):
        """Checks the position code sorting function to make
        sure it is sorting positions correctly for different schemes.
        This procedure is the same for mabs and tcrs, so no need
        to test separately."""
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs = [line.strip().split(",")[0] for line in fhandle]

        os.chdir(current_dir)
        random.seed(123)

        num_err = 0

        for scheme in ["martin", "imgt", "kabat", "aho"]:
            aligner = SingleChainAnnotator(chains=["H", "K", "L"], scheme=scheme)
            for seq in seqs:
                numbering = aligner.analyze_seq(seq)[0]
                numbering = [n for n in numbering if n != "-"]
                shuffled_numbering = copy.deepcopy(numbering)
                random.shuffle(shuffled_numbering)
                sorted_numbering = aligner.sort_position_codes(shuffled_numbering)
                if numbering != sorted_numbering:
                    num_err += 1

        self.assertTrue(num_err == 0)


if __name__ == "__main__":
    unittest.main()
