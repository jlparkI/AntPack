"""Contains tools needed to update the VJ gene database using
the latest result from IMGT. Note that this module uses Biopython,
but it is run when the package is built, not when the end user
installs, so Biopython is not required as a dependency."""
import sys
import os
import gzip
from datetime import date
from Bio import SeqIO
import wget
from . import utility_constants as UTC
from .special_functions import translate_imgt_nt



def update_vj_databases():
    """Downloads the reference sequences from IMGT and writes
    updated database files with the translated sequences. Each
    file indicates the last date it was updated."""
    current_dir = os.getcwd()
    start_path = os.path.abspath(os.path.dirname(__file__))

    # IMGT descriptions contain specific information in fields separated by
    # '|'. This list for readability provides the names of the fields in the
    # order they appear.
    imgt_fields = ["accession", "gene_name", "species", "functionality",
            "exon_region_label", "start_end", "num_nucleotides",
            "codon_start", "nucleotides_added_5'", "nucleotides_added_3'",
            "all_added_nucleotides", "num_amino_acids", "num_characters",
            "partial", "reverse"]

    # We save the created databases to the vj tools directory so they
    # can be retrieved easily by the appropriate tool.
    os.chdir(os.path.join(start_path, "..", "vj_tools", "consensus_data"))

    today = date.today()

    # Remove all outdated dbs.
    for f in os.listdir():
        if f.endswith(".fa.gz"):
            os.remove(f)


    for receptor_type, species_list in UTC.allowed_species.items():
        for species in species_list:
            for gene in UTC.genes[receptor_type]:
                target_url = f"https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/{species}/{receptor_type}/{receptor_type}{gene}.fasta"
                fasta_file = wget.download(target_url)
                ofname = f"{UTC.latin_to_common[species]}_{receptor_type}{gene}_{today}.fa.gz"

                if gene.endswith("V") and receptor_type == "IG":
                    gene_formatter = format_v_genes
                elif gene.endswith("J") and receptor_type == "IG":
                    gene_formatter = format_j_genes
                else:
                    raise RuntimeError("Only IG V and IG J genes are currently supported.")

                num_retained = 0

                with gzip.open(ofname, "wt", encoding="utf-8") as output_handle:
                    with open(fasta_file, "r", encoding="utf-8") as in_handle:
                        for seqrec in SeqIO.parse(in_handle, "fasta"):
                            nt_sequence = str(seqrec.seq)
                            description = seqrec.description
                            description_fields = dict(zip(imgt_fields,
                                [d.strip() for d in description.split("|")[:-1] ] ))

                            codon_start = 1
                            if description_fields["codon_start"] != "":
                                codon_start = int(description_fields["codon_start"])
                            elif gene.endswith("J"):
                                codon_start = check_problem_seqs(description_fields["gene_name"],
                                    species)

                            aa_seq = translate_imgt_nt(nt_sequence, description, codon_start)

                            # Skip entries that are partial, are not F for functional,
                            # that have no accession number, or contain "X" or "*".
                            if description_fields["accession"] == "":
                                continue
                            if description_fields["partial"] != "":
                                if 'partial' not in description_fields['partial']:
                                    print(f"Partial on {description_fields['partial']}")
                                continue
                            # ORF indicates open reading frame, P indicates pseudogene,
                            # (F) and [F] indicate this is not germline. Don't use
                            # any of these.
                            if description_fields["functionality"] != "F":
                                continue
                            if "X" in aa_seq or "*" in aa_seq:
                                print(aa_seq)
                                continue

                            fmt_seq = gene_formatter(aa_seq, description_fields["gene_name"])

                            num_retained += 1
                            _ = output_handle.write(f">{description_fields['gene_name']}\n{fmt_seq}\n")

                os.remove(fasta_file)
                print(f"For {ofname}, found {num_retained} valid sequences.")
    os.chdir(current_dir)


def format_v_genes(aa_seq, gene_name):
    """Takes as input an aa sequence and the corresponding
    description and formats the V gene, returning the formatted
    copy. This is fairly straightforward since V genes from IMGT
    are already numbered / aligned, so we merely need to trim at
    107 (three after the conserved cysteine at 104) and add blanks
    to bring the length up to 128. We also check that the conserved
    cysteines are present -- if not, there's something seriously
    wrong."""
    if aa_seq[103] != "C" or aa_seq[22] != "C":
        raise RuntimeError("Expected cysteines not present at expected locations "
                f"for {gene_name}")
    if len(aa_seq) < 100:
        raise RuntimeError(f"{gene_name} is unusually short for an aligned v-gene.")

    formatted_vgene = aa_seq[:108].ljust(128, "-").replace(".", "-")

    if len(formatted_vgene) != 128:
        raise RuntimeError(f"Error when formatting vgene {gene_name}")

    return formatted_vgene


def format_j_genes(aa_seq, gene_name):
    """Takes as input an aa sequence and the corresponding description
    and formats the J gene, returning the formatted copy. This is less
    straightforward since J genes are not aligned by default. The BEST
    way to do this would be to run a simple alignment against a template
    with a large gap opening penalty so that gaps are inserted only at
    the beginning and ending. For now, to save time, since we currently support
    only mouse and human and there is limited diversity in these species,
    we just add gaps based on the known length of the j-genes in IMGT for
    these species. This very simplistic procedure should be replaced on
    the next db update."""
    receptor = gene_name[:4]
    if receptor in ["IGKJ", "IGLJ"]:
        formatted_jgene = aa_seq.replace(".", "-") + "-"
    elif receptor == "IGHJ":
        formatted_jgene = aa_seq.replace(".", "-")[-14:]
    else:
        raise RuntimeError("Currently unsupported receptor supplied.")

    formatted_jgene = formatted_jgene.rjust(128, "-")
    return formatted_jgene


def check_problem_seqs(gene_name, species):
    """There are some known problem_seqs where a codon start SHOULD be specified
    but is not (this information is missing in the IMGT db but should be specified).
    We try to look these up and 'patch' these here with a dictionary mapping
    gene name to codon start. Definitely not ideal...but this is a shortcoming
    of the IMGT db."""
    known_problem_seqs = {"Homo_sapiens":{
        "IGHJ6*03":3
        }}
    if gene_name not in known_problem_seqs[species]:
        raise RuntimeError(f"Unexpected error encountered for gene {gene_name}. "
                "Check the IMGT db.")
    return known_problem_seqs[species][gene_name]
