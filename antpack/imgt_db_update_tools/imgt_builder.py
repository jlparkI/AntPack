"""
Contains the functions needed to download selected IMGT files and convert them
into updated profiles for numbering. The package contains a standard set of these
profiles, but this tool enables the user to build updated profiles (if the IMGT db has
for example been updated). There are two kinds of profiles which are generated:
HMM profiles for scoring sequences using HMMer AND consensus profiles used for the PSSM-
based alignment in the C++ extension. We prefer the latter and may eventually remove
the HMMer-based scoring if it turns out to be less useful.

The parent sequences are found here:

https://www.imgt.org/vquest/refseqh.html     
"""

import os
import requests
from .imgt_html_parser import IMGT_HTML_Parser
from .alignment_formatting_tools import build_stockholm_alignments
from ..constants import hmmbuild_constants as hmbc




def build_consensus_alignment(output_path, muscle_fpath, cleanup = True):
    """Constructs consensus alignment schemes for all species and chains.

    output_path (str): Filename of a folder where the output HMM
        will be constructed. Should be a directory.
    muscle_fpath (str): The filepath to an installed instance of the muscle sequence
        aligner.
    cleanup (bool): If True, intermediate files generated while building the consensus
        sequence are deleted. If False they are retained.
    """
    current_dir = os.getcwd()
    try:
        os.chdir(output_path)
        os.chdir(current_dir)
    except Exception as exc:
        raise ValueError("Invalid output file path supplied.") from exc

    #Note that currently we generate for all species and chain types, although
    #we may relax this in future.
    genes, selected_species = validate_inputs("ALL", "ALL")
    cleanup_files = []

    for gene in genes:
        for target_species in selected_species:
            if gene in hmbc.genes["TR"] and target_species not in hmbc.allowed_species["TR"]:
                continue
            if gene[0] in "KL" and target_species == "Vicugna+pacos":
                continue
            print(f"Now working on {target_species}, {gene}")
            filename = os.path.join(output_path,
                    f"{target_species.replace('+', '_')}_{gene}.fasta" )
            cleanup_files.append(filename)
            if os.path.isfile(filename):
                print(f"{target_species}, {gene} already downloaded")
                continue

            if retrieve_fasta(target_species, gene, output_path):
                print(f"******\nWARNING!!! Error retrieving {target_species}, {gene}. "
                "This species & gene will not be represented in the consensus.")
            else:
                print(f"Downloaded {target_species}, {gene}")

    if os.path.isfile(os.path.join(output_path, "ALL_ALIGNED.stockholm")):
        os.remove(os.path.join(output_path, "ALL_ALIGNED.stockholm"))
    build_stockholm_alignments(output_path, selected_species, muscle_fpath)
    if cleanup:
        for cleanup_file in cleanup_files:
            os.remove(cleanup_file)
        os.remove(os.path.join(output_path, "all_js_aligned.fasta"))
        os.remove(os.path.join(output_path, "all_js.fasta"))



def validate_inputs(species, receptor_type):
    """Checks user specified options for validity, and returns a list
    of urls to retrieve.

    Raises:
        ValueError: A ValueError is raised if user selections are
            not valid.
    """
    if receptor_type in hmbc.genes:
        genes = hmbc.genes[receptor_type]
    else:
        raise ValueError("Unrecognized receptor type supplied. "
                    "Should be one of 'ALL', 'IG', 'TR'.")

    if isinstance(species, str):
        if species == "ALL":
            selected_species = hmbc.allowed_species[receptor_type]
        else:
            raise ValueError("Unexpected species supplied.")
    elif isinstance(species, list):
        for species_name in species:
            if species_name not in hmbc.allowed_species[receptor_type]:
                raise ValueError(f"For {receptor_type}, only species "
                    f"{hmbc.allowed_species[receptor_type]} are allowed.")
        selected_species = species

    return genes, selected_species



def retrieve_fasta(species, gene_type, target_dir):
    """ Retrieve the fasta sequences for a species and gene type from IMGT."""
    url = hmbc.urls[gene_type]%species
    resp = requests.get(url, timeout=300)
    if resp.status_code >= 300 or resp.status_code < 200:
        return 1
    parser = IMGT_HTML_Parser()
    sequences = parser.parse_imgt_html(resp.text)
    if sequences:
        filename = os.path.join(target_dir, f"{species.replace('+', '_')}_{gene_type}.fasta" )
        with open(filename, "w+", encoding="utf-8") as fhandle:
            for name, sequence in sequences:
                fhandle.write(f">{name}\n{sequence}\n")
    else:
        return 1
    return 0
