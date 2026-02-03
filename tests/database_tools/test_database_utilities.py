"""Test database utility functions that the end user or
Python functions may need to access."""
import os
import sys
import gzip
import struct
import sqlite3
import pytest
import numpy as np
from antpack import (number_imgt_imgt_cdr, SingleChainAnnotator)
from antpack.utilities import read_fasta
from .database_utilities import (get_vgene_code,
                setup_canonical_numbering, int_to_bin)
from ..conftest import (get_test_data_filepath, get_test_base_filepath)


def test_mab_renumbering(load_mab_nmbr_test_data):
    """Run a batch of test data (approximately 1600 sequences from the
    PDB) to ensure that numbering isolated cdrs is consistent with
    numbering the full sequence and extracting the cdr."""
    seqs, numberings = load_mab_nmbr_test_data
    numbering = numberings["imgt"]

    aligner = SingleChainAnnotator(chains=["H", "K", "L"],
        scheme="imgt")
    for seq, nmbr in zip(seqs, numbering):
        # Chain here is irrelevant since we are using IMGT.
        labels = aligner.assign_cdr_labels(nmbr, "H", "imgt")
        for i in range(1, 4):
            extracted_cdr, extracted_tokens = [], []
            for s, c, n in zip(seq, labels, nmbr):
                if c == f"cdr{i}":
                    extracted_cdr.append(s)
                    extracted_tokens.append(n)

            assigned_tokens = number_imgt_imgt_cdr(
                    ''.join(extracted_cdr), i)
            assert assigned_tokens == extracted_tokens
