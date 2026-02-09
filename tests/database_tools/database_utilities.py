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




def get_vgene_code(input_vgene, species):
    """Converts the input vgene and species to a code."""
    vgene = input_vgene.split("_")[0]

    species_map = {"human":1, "mouse":2, "alpaca":3, "rabbit":4}
    if species not in species_map:
        species_code = 0
    else:
        species_code = species_map[species]

    codemap = {"A":1, "G":2, "H":3, "L":4, "K":5, "B":6, "D":7}
    if vgene[2] not in codemap:
        raise RuntimeError("Incorrect chain code found in test data.")
    chain_code = codemap[vgene[2]]

    family, gene, allele = "", "", ""

    if "S" in vgene:
        family = vgene[4:vgene.find("S")+1]
        if "*" in vgene:
            gene = vgene[vgene.find("S")+1:vgene.find("*")]
            allele = vgene[vgene.find("*")+1:]
        else:
            gene = vgene[vgene.find("S")+1:]
    elif "-" in vgene:
        family = vgene[4:vgene.find("-")]
        if "*" in vgene:
            gene = vgene[vgene.find("-")+1:vgene.find("*")]
            allele = vgene[vgene.find("*")+1:]
        else:
            gene = vgene[vgene.find("-")+1:]
    elif "*" in vgene:
        family = vgene[4:vgene.find("*")]
        allele = vgene[vgene.find("*")+1:]
    else:
        family = vgene[4:]


    return ([chain_code, min(extract_numeric(family), 255),
             min(extract_numeric(gene), 255), species_code],
            family, gene, allele)




def extract_numeric(input_string):
    """Extracts the first numeric portion of an input
    string."""
    read_now = False
    extracted_num = ""

    for letter in input_string:
        if letter.isnumeric():
            if not read_now:
                read_now = True
            extracted_num += letter
        elif read_now:
            break

    if len(extracted_num) == 0:
        return 0
    return int(extracted_num)






def int_to_bin(int_val, num_bytes):
    """Converts an integer value to a binary
    string of length num bytes."""
    binstring = format(int_val, 'b')
    return '0'*(num_bytes*8 - len(binstring)) + binstring
