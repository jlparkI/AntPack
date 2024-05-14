"""Tests basic functionality for the VJGeneTool."""
import os
import copy
import gzip
import unittest
from antpack import VJGeneTool

class TestVJGeneTool(unittest.TestCase):


    def test_error_checking(self):
        """Check that sequences which have known issues are flagged
        as such, and that deliberately invalid inputs are recognized."""
        # Pass dummy sequences with errors.
        vj_tool = VJGeneTool()
        with self.assertRaises(RuntimeError):
            vj_tool.assign_numbered_sequence("AYAYAYA",
                    ["1", "2", "3"], "H")
        with self.assertRaises(RuntimeError):
            vj_tool.assign_numbered_sequence("AYAYAYA",
                ["1", "2", "3", "4", "5", "6", "7"], "Z")
        with self.assertRaises(RuntimeError):
            vj_tool.assign_numbered_sequence("BYBYBYY",
                ["1", "2", "3", "4", "5", "6", "7"], "H")

        results = vj_tool.assign_sequence("AYAYAYA")
        self.assertTrue(results[0] is None)
        self.assertTrue(results[1] is None)



    def test_percent_ident_calc(self):
        """Double checks the percent identity calculation
        done by the cpp extension using a simple stupid
        Python version, just to make sure it's correct."""
        vj_tool = VJGeneTool()

        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            test_seqs = [line.strip().split(",")[0] for line in fhandle]

        seq_dict = {}

        for chain in ['H', 'K', 'L']:
            for recep in [f"IG{chain}V", f"IG{chain}J"]:
                seqs, names = vj_tool.vj_gene_matchups["human"][recep].getSeqLists()
                seq_dict[recep] = {}
                seq_dict[recep]["seqs"] = seqs
                seq_dict[recep]["names"] = names

        os.chdir(current_dir)
        for seq in test_seqs[:100]:
            numbering, _, chain, _ = vj_tool.default_aligner.analyze_seq(seq)
            fmt_seq = vj_tool._prep_sequence(seq, numbering)

            for recep in [f"IG{chain}V", f"IG{chain}J"]:
                gpred, id_pred = vj_tool.vj_gene_matchups["human"][recep].vjMatch(fmt_seq)

                best_pid, matchnum = 0, 0
                for i, template in enumerate(seq_dict[recep]["seqs"]):
                    matches = 0

                    for qletter, tletter in zip(fmt_seq, template):
                        if qletter == tletter:
                            matches += 1

                    true_pid = float(matches) / float(len(template.replace("-", "")))
                    if true_pid > best_pid:
                        matchnum = i
                        best_pid = copy.deepcopy(true_pid)

                gtrue = seq_dict[recep]["names"][matchnum]

                self.assertTrue(gtrue == gpred)
                self.assertTrue(best_pid == id_pred)


if __name__ == "__main__":
    unittest.main()
