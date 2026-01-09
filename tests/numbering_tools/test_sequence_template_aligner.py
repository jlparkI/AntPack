"""Tests SequenceTemplateAligners to ensure that
new sequences can be added to existing msas in
a straightforward way, and that the aligner used
for database search is functioning correctly."""
import os
import gzip
import pytest
from Bio import SeqIO
from antpack import SingleChainAnnotator as SCA
from antpack.antpack_cpp_ext import (SequenceTemplateAligner)

from ..conftest import get_test_data_filepath



def test_template_aligner_error_checking(load_antibody_data):
    """Test that the template aligner raises an exception if initialized
    with unacceptable parameters."""
    seqrecs = load_antibody_data
    _, position_codes, _ = get_chain_data("L", seqrecs)
    with pytest.raises(RuntimeError):
        _ = SequenceTemplateAligner(position_codes, "L", "")
    with pytest.raises(RuntimeError):
        _ = SequenceTemplateAligner(position_codes, "Z", "imgt")
    with pytest.raises(RuntimeError):
        _ = SequenceTemplateAligner(position_codes, "H", "imgt", "gollum")
    with pytest.raises(RuntimeError):
        _ = SequenceTemplateAligner(position_codes + ["Z"], "H", "imgt", "")



def test_correct_alignment_and_slicing(load_antibody_data):
    """Tests that the template aligner can load and slice sequences
    correctly."""
    sc = SCA()
    seqrecs = load_antibody_data

    for chain_type in ["H", "L"]:
        seqs, position_codes, region_labels = get_chain_data(chain_type, seqrecs)
        unaligned_seqs = [s.replace('-', '') for s in seqs]
        unaligned_seq_nmbr = sc.analyze_seqs(unaligned_seqs)

        sta = SequenceTemplateAligner(position_codes,
                    chain_type, "imgt", "")
        realigned_seqs = [sta.align_sequence(s, u[0], False) for u, s in zip(
            unaligned_seq_nmbr, unaligned_seqs)]
        for r, s in zip(realigned_seqs, seqs):
            assert r==s

        for region in ["all", "fmwk", "cdr2"]:
            realigned_seqs = [sta.align_and_slice_sequence(s, u[0], region, False)
                    for u, s in zip(unaligned_seq_nmbr, unaligned_seqs)]
            if region == "all":
                sliced_seqs = seqs
            else:
                sliced_seqs = [''.join([l for token, l in zip(region_labels, s)
                    if token.startswith(region)])
                        for s in seqs]
            for r, s in zip(realigned_seqs, sliced_seqs):
                assert r==s

            if region == "all":
                expected_region_size = len(seqs[0])
            else:
                expected_region_size = len([token for token in region_labels
                    if token.startswith(region)])
            assert expected_region_size==sta.get_region_size(region)

        assert sta.get_template_length() == len(position_codes)



def test_forward_reverse_numbering(load_antibody_data):
    """Check to make sure that the forward and reverse numbering
    returned by the template aligner makes sense."""
    sc = SCA()
    seqrecs = load_antibody_data

    for chain_type in ["H", "L"]:
        seqs, position_codes, _ = get_chain_data(chain_type, seqrecs)
        unaligned_seqs = [s.replace('-', '') for s in seqs]
        unaligned_seq_nmbr = sc.analyze_seqs(unaligned_seqs)

        pcode_dict = {p:i for i, p in enumerate(position_codes)}
        sta = SequenceTemplateAligner(position_codes,
                    chain_type, "imgt", "")

        for unaligned_seq, nmbr in zip(
                unaligned_seqs, unaligned_seq_nmbr):
            gt_fwd_numbering = [pcode_dict[p] for p in nmbr[0]]
            assert (gt_fwd_numbering ==
                    sta.retrieve_alignment_forward_numbering(
                        unaligned_seq, nmbr[0], False))

            reverse_dict = {p:i for i, p in enumerate(nmbr[0])}
            gt_rev_numbering = []
            for p in position_codes:
                if p in reverse_dict:
                    gt_rev_numbering.append(reverse_dict[p])
                else:
                    gt_rev_numbering.append(-1)
            assert (gt_rev_numbering ==
                    sta.retrieve_alignment_back_numbering(
                        unaligned_seq, nmbr[0], False))



@pytest.fixture
def load_antibody_data(get_test_data_filepath):
    """Loads some saved test data specific for antibodies."""
    data_path = os.path.join(get_test_data_filepath, "addtnl_test_data.fasta.gz")
    with gzip.open(data_path, "rt") as fhandle:
        seqrecs = list(SeqIO.parse(fhandle, "fasta"))
    return seqrecs


def get_chain_data(chain_type, seqrecs):
    """Extracts data specific for a given chain from
    the fixture-generated data."""
    sc = SCA()
    if chain_type == "L":
        light_seqs = [str(seqrec.seq) for seqrec in seqrecs if "light" in seqrec.name]
        light_alignments = sc.analyze_seqs(light_seqs)
        lpos_codes, lmsa = sc.build_msa(light_seqs, light_alignments)
        region_labels = sc.assign_cdr_labels(lpos_codes, "L")
        return lmsa, lpos_codes, region_labels
    heavy_seqs = [str(seqrec.seq) for seqrec in seqrecs if "heavy" in seqrec.name]
    heavy_alignments = sc.analyze_seqs(heavy_seqs)
    hpos_codes, hmsa = sc.build_msa(heavy_seqs, heavy_alignments)
    region_labels = sc.assign_cdr_labels(hpos_codes, "H")
    return hmsa, hpos_codes, region_labels
