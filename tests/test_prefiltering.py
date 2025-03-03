"""Tests basic functionality for prefiltering routines."""
import os
import unittest
import gzip
import numpy as np
from antpack.numbering_tools.cterm_finder import PyCTermFinder
from antpack.numbering_tools.cterm_finder import _load_nterm_kmers

KMER_WIDTH = 9
KMER_WINDOW = 14


class TestPrefiltering(unittest.TestCase):


    def test_pycterm_routine(self):
        """Compare scoring for N and C terminals with the
        wrapped C++ code with an inefficient but simple Python
        routine to ensure results match."""
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs = [line.strip().split(",")[0] for line
                    in fhandle]

        os.chdir(current_dir)

        aa_to_idx = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N",
            "P", "Q", "R", "S", "T", "V", "W", "Y"]
        aa_to_idx = {a:i for i,a in enumerate(aa_to_idx)}

        pyct_finder = PyCTermFinder()
        kmer_dict = _load_nterm_kmers()
        cterm_to_aa = np.zeros((10,21,3))

        for k, chain in enumerate(["H", "K", "L"]):
            cterm_to_aa[:,:,k] = np.load(
                    os.path.join(project_path, "..", "src",
                "antpack", "numbering_tools", "consensus_data",
                "mabs", f"CTERMFINDER_CONSENSUS_{chain}.npy"))

        for seq in seqs:
            nterm_gt_scores, cterm_gt_scores = \
                    np.zeros((3), dtype=np.int32), \
                    np.zeros((3))

            kmer_presence = np.zeros((len(seq), 3),
                        dtype=np.int32)
            for i, _ in enumerate(seq[:-KMER_WIDTH]):
                if seq[i:i+KMER_WIDTH] in kmer_dict:
                    kmer_presence[i,
                            kmer_dict[seq[i:i+KMER_WIDTH]]] += 1

            for i, _ in enumerate(seq[:-5]):
                cterm_scores = np.zeros((3))
                len_to_traverse = min(cterm_to_aa.shape[0],
                        len(seq) - i)
                for j in range(len_to_traverse):
                    for k in range(3):
                        cterm_scores[k] += cterm_to_aa[j,
                                aa_to_idx[seq[i+j]], k]
                for k in range(3):
                    if cterm_scores[k] > cterm_gt_scores[k]:
                        cterm_gt_scores[k] = cterm_scores[k]

            for i in range(kmer_presence.shape[0] -
                    KMER_WINDOW):
                for k in range(3):
                    score = kmer_presence[i:i+KMER_WINDOW,k].sum()
                    if score > nterm_gt_scores[k]:
                        nterm_gt_scores[k] = score

            nterm_test_scores, cterm_test_scores = \
                    np.zeros((3), dtype=np.int32), \
                    np.zeros((3))
            nterm_test_positions, cterm_test_positions = \
                    np.zeros((3), dtype=np.int32), \
                    np.zeros((3), dtype=np.int32)
            pyct_finder.pyfind_c_terminals(seq, cterm_test_scores,
                    cterm_test_positions, nterm_test_scores,
                    nterm_test_positions)
            self.assertTrue(np.allclose(nterm_gt_scores,
                nterm_test_scores))
            self.assertTrue(np.allclose(cterm_gt_scores,
                cterm_test_scores))


if __name__ == "__main__":
    unittest.main()
