"""Fixtures for providing larger external datasets
to all test functions."""
import os
import gzip
import pytest
import numpy as np
from Bio import SeqIO
from antpack import SingleChainAnnotator


@pytest.fixture(scope="module")
def get_test_base_filepath():
    """Gets the base filepath for the project regardless
    of which test calls."""
    return os.path.join( os.path.abspath(os.path.dirname(__file__)),
                        ".." )


@pytest.fixture(scope="module")
def get_test_data_filepath():
    """Gets the base filepath for test data regardless
    of which test calls."""
    return os.path.join( os.path.abspath(os.path.dirname(__file__)),
                        "test_data" )

@pytest.fixture(scope="module")
def standard_aa_list():
    """Returns a standard list of amino acids sorted alphabetically."""
    return ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N",
        "P", "Q", "R", "S", "T", "V", "W", "Y", "-"]


@pytest.fixture(scope="module")
def load_mab_nmbr_test_data(get_test_data_filepath):
    """Loads test data for ensuring that numbering is accurate
    for different schemes."""
    with gzip.open(os.path.join(get_test_data_filepath,
                                "test_data.csv.gz"), "rt") as fhandle:
        _ = fhandle.readline()
        nmbr = {scheme:[] for scheme in ["martin", "imgt",
                                              "kabat", "aho"]}
        seqs = []
        for line in fhandle:
            line_elements = line.strip().split(",")
            seqs.append(line_elements[0])
            nmbr["martin"].append(line_elements[1].split("_"))
            nmbr["imgt"].append(line_elements[2].split("_"))
            nmbr["kabat"].append(line_elements[3].split("_"))
            nmbr["aho"].append(line_elements[4].split("_"))
    return seqs, nmbr


@pytest.fixture(scope="module")
def load_tcr_nmbr_test_data(get_test_data_filepath):
    """Loads test data for ensuring that numbering is accurate
    for TCRs."""
    with gzip.open(os.path.join(get_test_data_filepath,
                                "tcr_test_data.csv.gz"), "rt") as fhandle:
        _ = fhandle.readline()
        seqs, nmbr = [], []
        for line in fhandle:
            line_elements = line.strip().split(",")
            seqs.append(line_elements[0])
            nmbr.append(line_elements[1].split("_"))
    return seqs, nmbr


@pytest.fixture
def load_mab_emclust_test_data(get_test_data_filepath, standard_aa_list):
    """Load some test data primarily used for testing EM
    clustering."""
    data_path = os.path.join(get_test_data_filepath,
            "addtnl_test_data.fasta.gz")
    with gzip.open(data_path, "rt") as fhandle:
        seqrecs = list(SeqIO.parse(fhandle, "fasta"))

    # For these purposes, just use the light chains.
    sc = SingleChainAnnotator()
    light_seqs = [str(seqrec.seq) for seqrec in seqrecs if "light" in seqrec.name]
    light_alignments = sc.analyze_seqs(light_seqs)
    lpos_codes, lmsa = sc.build_msa(light_seqs, light_alignments)

    # Duplicate each light seq a bunch of times to make it easier to test
    # multithreading.
    lmsa = 5 * lmsa

    # Encode the light seqs using a simple stupid Python routine
    # for comparison with the model's build-in encoder.
    aa_codes = {k:i for i,k in enumerate(standard_aa_list)}
    encoded_data = np.empty((len(lmsa), len(lmsa[0])),
            dtype=np.uint8)

    for i, seq in enumerate(lmsa):
        for j, letter in enumerate(seq):
            encoded_data[i,j] = aa_codes[letter]

    region_labels = sc.assign_cdr_labels(lpos_codes, "L")
    return lmsa, lpos_codes, encoded_data, region_labels


@pytest.fixture
def load_non_mab_test_data(get_test_data_filepath, standard_aa_list):
    """Loads some non antibody non tcr test data."""
    data_path = os.path.join(get_test_data_filepath,
            "non_antibody_test_data", "converted_seqs.txt.gz")
    with gzip.open(data_path, "rt") as fhandle:
        seqs = [l.strip() for l in fhandle]

    aa_codes = {k:i for i,k in enumerate(standard_aa_list)}
    encoded_data = np.empty((len(seqs), len(seqs[0])),
            dtype=np.uint8)

    for i, seq in enumerate(seqs):
        for j, letter in enumerate(seq):
            encoded_data[i,j] = aa_codes[letter]
    return seqs, encoded_data
