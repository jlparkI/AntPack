"""Tests basic functionality for prefiltering routines."""
import os
import pytest
import numpy as np
from antpack.numbering_tools.cterm_finder import PyCTermFinder
from antpack.numbering_tools.cterm_finder import _load_nterm_kmers

from ..conftest import (standard_aa_list,
            get_test_base_filepath, load_mab_nmbr_test_data)


def test_pycterm_routine(load_mab_nmbr_test_data, get_test_base_filepath,
                         standard_aa_list, standard_window_settings):
    """Compare scoring for N and C terminals with the
    wrapped C++ code with an inefficient but simple Python
    routine to ensure results match."""
    seqs, _ = load_mab_nmbr_test_data
    aa_to_idx = standard_aa_list
    aa_to_idx = {a:i for i,a in enumerate(aa_to_idx)}

    pyct_finder = PyCTermFinder()
    kmer_dict = _load_nterm_kmers()
    cterm_to_aa = np.zeros((10,21,3))

    for k, chain in enumerate(["H", "K", "L"]):
        cterm_to_aa[:,:,k] = np.load(
                os.path.join(get_test_base_filepath, "src",
            "antpack", "numbering_tools", "consensus_data",
            "mabs", f"CTERMFINDER_CONSENSUS_{chain}.npy"))

    for seq in seqs:
        nterm_gt_scores, cterm_gt_scores = \
                np.zeros((3), dtype=np.int32), \
                np.zeros((3))

        kmer_presence = np.zeros((len(seq), 3),
                    dtype=np.int32)
        kmer_width, kmer_window = standard_window_settings

        for i, _ in enumerate(seq[:-kmer_width]):
            if seq[i:i+kmer_width] in kmer_dict:
                kmer_presence[i,
                        kmer_dict[seq[i:i+kmer_width]]] += 1

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

        for i in range(kmer_presence.shape[0] - kmer_width):
            for k in range(3):
                score = kmer_presence[i:i+kmer_window,k].sum()
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
        assert np.allclose(nterm_gt_scores,
            nterm_test_scores)
        assert np.allclose(cterm_gt_scores,
            cterm_test_scores)


@pytest.fixture
def standard_window_settings():
    """Returns standard kmer window settings for the test."""
    width = 9
    window_size = 14
    return width, window_size
