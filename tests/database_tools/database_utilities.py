"""Utilities shared by all database tests."""
from antpack.antpack_cpp_ext import (SequenceTemplateAligner,
        return_imgt_canonical_numbering_cpp)


def setup_canonical_numbering(numbering_scheme,
        cdr_scheme, chain_type):
    """Sets up a template aligner for the numbering
    scheme we have chosen."""
    numbering = return_imgt_canonical_numbering_cpp()

    sta = SequenceTemplateAligner(numbering,
                        chain_type, numbering_scheme,
                        cdr_scheme)
    return sta, numbering




def get_vgene_code(vgene, species):
    """Converts the input vgene and species to a code."""
    species_map = {"human":0, "mouse":1, "alpaca":2, "rabbit":3}
    if species not in species_map:
        species_code = 255
    else:
        species_code = species_map[species]
    
    codemap = {"A":0, "G":1, "H":2, "L":3, "K":4, "B":5, "D":6}
    if vgene[2] not in codemap:
        raise RuntimeError("Incorrect chain code found in test data.")
    chain_code = codemap[vgene[2]]

    extracted_nums = []
    read_now = False
    current_num = ""

    for i in range(4, len(vgene)):
        if vgene[i] == "*":
            break
        if vgene[i].isnumeric():
            if not read_now:
                read_now = True
                current_num = ""
            current_num += vgene[i]
        elif read_now:
            read_now = False
            if len(current_num) > 0:
                extracted_nums.append(int(current_num))
            current_num = ""

    if len(current_num) > 0:
        extracted_nums.append(int(current_num))
    if len(extracted_nums)==1:
        extracted_nums.append(255)
    elif len(extracted_nums) < 1:
        raise RuntimeError("Invalid vgene found in test data.")

    return (chain_code, species_code, extracted_nums[0],
            extracted_nums[1])



def int_to_bin(int_val, num_bytes):
    """Converts an integer value to a binary
    string of length num bytes."""
    binstring = format(int_val, 'b')
    return '0'*(num_bytes*8 - len(binstring)) + binstring
