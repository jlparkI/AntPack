"""Tests the SequenceTemplateAligner to ensure that
new sequences can be added to existing msas in a straightforward
way."""
import os
import gzip
import unittest
from Bio import SeqIO
import numpy as np
from antpack import SingleChainAnnotator as SCA
from antpack.antpack_cpp_ext import SequenceTemplateAligner


class TestSequenceTemplateAligner(unittest.TestCase):


    def test_error_checking(self):
        """Ensure that absurd inputs generate expected
        errors."""
        _, position_codes, _ = load_antibody_test_data()
        with self.assertRaises(RuntimeError):
            sta = SequenceTemplateAligner(position_codes,
                    "L", "")

        with self.assertRaises(RuntimeError):
            sta = SequenceTemplateAligner(position_codes,
                    "Z", "imgt")

        with self.assertRaises(RuntimeError):
            sta = SequenceTemplateAligner(position_codes,
                    "H", "imgt", "gollum")

        with self.assertRaises(RuntimeError):
            sta = SequenceTemplateAligner(position_codes + ['Z'],
                    "H", "imgt", "")


    def test_alignment_and_slicing(self):
        """Check to make sure the template aligner correctly
        aligns and slices input sequences. We assume here that
        cdr label assignment is correct since this is tested
        under test_annotator_extended_functions."""
        sc = SCA()
        seqs, position_codes, region_labels = load_antibody_test_data()
        unaligned_seqs = [s.replace('-', '') for s in seqs]
        unaligned_seq_nmbr = sc.analyze_seqs(unaligned_seqs)

        sta = SequenceTemplateAligner(position_codes,
                    "L", "imgt", "")
        realigned_seqs = [sta.align_sequence(s, u[0], False) for u, s in zip(
            unaligned_seq_nmbr, unaligned_seqs)]
        for r, s in zip(realigned_seqs, seqs):
            self.assertTrue(r==s)

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
                self.assertTrue(r==s)

            if region == "all":
                expected_region_size = len(seqs[0])
            else:
                expected_region_size = len([token for token in region_labels
                    if token.startswith(region)])
            self.assertTrue(expected_region_size==sta.get_region_size(region))

        self.assertTrue(sta.get_template_length() ==
                len(position_codes))




def load_antibody_test_data():
    """Loads some saved test data specific for antibodies."""
    start_dir = os.path.abspath(os.path.dirname(__file__))

    data_path = os.path.join(start_dir, "test_data",
            "addtnl_test_data.fasta.gz")
    with gzip.open(data_path, "rt") as fhandle:
        seqrecs = list(SeqIO.parse(fhandle, "fasta"))

    # For purposes of this test, just use the light chains.
    sc = SCA()
    light_seqs = [str(seqrec.seq) for seqrec in seqrecs if "light" in seqrec.name]
    light_alignments = sc.analyze_seqs(light_seqs)
    lpos_codes, lmsa = sc.build_msa(light_seqs, light_alignments)

    region_labels = sc.assign_cdr_labels(lpos_codes, "L")
    return lmsa, lpos_codes, region_labels


if __name__ == "__main__":
    unittest.main()
