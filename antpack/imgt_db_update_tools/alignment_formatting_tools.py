"""Contains tools needed to convert a set of fasta files with v, j etc. into
alignments that can be fed to HMMER to build an HMM profile. We use the
muscle software for alignment of j genes with a very high gapopen penalty (-10),
then trim the alignment before IMGT position 116. All v are combined pairwise with
all j to create all possible recombinants.
"""
import os
import sys
from copy import deepcopy
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from ..constants import allowed_inputs
from ..constants import hmmbuild_constants as hmbc


def read_fasta_input(input_file, read_all=False, region_name=""):
    """Reads an aligned or unaligned fasta file and saves key fields plus
    the sequence info while discarding the unnecessary information provided by
    IMGT."""
    imgt_fields =  ["accession_number",
    "allele",  
    "species",  
    "functionality",  
    "region",  
    "start_and_end_positions_IMGT/LIGM-DB_accession_number", 
    "number_of_nucleotides", 
    "codon_start", 
    "number_nucleotides_added_in_5'_compared_IMGT/LIGM-DB", 
    "number_nucleotides_added_in_3'_compared_IMGT/LIGM-DB", 
    "number_nucleotides_to_correct_sequencing_errors", 
    "number_of_amino_acids", 
    "number_of_characters", 
    "partial",  
    "reverse"]

    extended_aas = deepcopy(allowed_inputs.allowed_amino_acids)
    extended_aas.add(".")

    records = {}
    try:
        handle = open(input_file, "r", encoding="utf-8")
        handle.close()
    except IOError:
        print('Warning file', input_file, 'could not be found')
        return records

    region=""
    with open(input_file, "r", encoding="utf-8") as fhandle:
        for record in SeqIO.parse(fhandle, "fasta"):
            fields = dict(list(zip( imgt_fields, record.description.strip(">").split("|"))) )
            sequence = str(record.seq)
            for required_key in ["accession_number", "functionality", "partial", "reverse",
                    "region", "species", "allele"]:
                if required_key not in fields:
                    raise ValueError(f"Missing data for file {input_file}; "
                            "required field {required_key} not found.")

            if fields['accession_number'] == 'None':
                continue
            if fields["functionality"]=="F" and not fields["partial"].strip() and \
                    not fields["reverse"].strip():
                if fields["allele"].split("*")[-1].strip()!="01" and not read_all:
                    continue
                if len([aa for aa in sequence if aa not in extended_aas]) > 0:
                    continue

                if fields["region"] != region_name:
                    if region_name.startswith("C"):
                        if len(sequence) < 100:
                            continue
                    if region:
                        if fields["region"] != region:
                            raise ValueError(f"For file {input_file}, some of the "
                                "entries have different regions.")

                region=fields["region"]
                records[ (fields["species"], fields[ "allele" ] ) ] = sequence
    return records

def write_fasta(sequences, output_filename):
    """This is duplicate to Biopython's SeqIO.write but is retained here
    since correctly formatting the name and description would add some
    slight additional complexity if using SeqIO.write. This function writes
    the input sequence dictionary to file."""
    with open(output_filename, "w+", encoding="utf-8") as fhandle:
        for (species, chain_type) in sequences:
            for (species2, allele) in sequences[(species, chain_type)]:
                fhandle.write(f">{species}|{chain_type}|{species2}|{allele}\n")
                fhandle.write(f"{sequences[(species, chain_type)][(species2, allele)]}\n")



def format_j_genes(jalignments, target_dir, muscle_fpath):
    """Applies special formatting for the j genes, for which it is more
    difficult to build an alignment. Using a large gap open penalty
    prevents muscle from opening gaps in multiple locations."""

    reference = ("WFAYWGQGTLVTVSA", 4  , 19 )
    #                 seq           start  end

    input_file = os.path.join(target_dir, "all_js.fasta" )
    write_fasta(jalignments, input_file)

    al_filename = os.path.join(target_dir, "all_js_aligned.fasta" )
    muscle_cline = MuscleCommandline(muscle_fpath, input=input_file, out=al_filename, gapopen=-10.0)
    _ = muscle_cline()
    #No error checking here, since MuscleCommandline will automatically throw an exception if muscle
    #returns an error.
    with open(al_filename, "r", encoding="utf-8") as fhandle:
        aligned = [(s.description.strip(">"), str(s.seq)) for s in SeqIO.parse(fhandle, "fasta")]
    new_jalignments = {}

    # Find the reference sequence and what we need to do to map
    for name, sequence in aligned:
        if name == "Mus|H|Mus musculus|IGHJ3*01":
            ref_aligned = sequence
            break
    start = ref_aligned.index( reference[0] )
    if start > reference[1]:
        START = start+1-reference[1]
    else:
        START = 0
    END = start + 15


    for name, sequence in aligned:
        species, chain_type, id1, id2 =  name.strip(">").split("|")

        if (species, chain_type) not in new_jalignments:   
            new_jalignments[(species, chain_type)] = {}
        # We take the last 13 of the new alignment and pad into 20 long string 
        new_jalignments[(species, chain_type)][ (id1, id2) ] = sequence[START: END][-14:].rjust(20).replace(" ", ".")
    return new_jalignments

def format_v_genes(valignments):
    """
    Take upto and including imgt position 108 in the alignment. Pad with gaps on the right side
    """

    new_valignments = {}
    for (species, chain_type) in valignments:
        new_valignments[(species, chain_type)] = {}
        for seq in valignments[(species, chain_type)]:
            sequence = valignments[(species, chain_type)][seq]
            common_species_name = hmbc.latin_to_common[species]
            #Unfortunately the IMGT database seems to have some alignment issues that
            #necessitate some manual fixes.
            if common_species_name == "rhesus":
                if chain_type == "H":
                    sequence = fix_rhesus_heavy(sequence)
                elif chain_type == "L":
                    sequence = fix_rhesus_lambda(sequence)
                elif chain_type == "K":
                    sequence = fix_rhesus_kappa(sequence)
            elif common_species_name == "mouse":
                if chain_type == "A":
                    sequence = fix_mouse_alpha(sequence)
                if chain_type == "D":
                    sequence = fix_mouse_delta(sequence)
            elif common_species_name == "rat":
                if chain_type == "H":
                    sequence = fix_rat_heavy(sequence)
            new_valignments[(species, chain_type)][seq] = \
                    sequence[:108].ljust( 108 ).replace(" ",".")
            if new_valignments[(species, chain_type)][seq][103] != "C" or \
                    new_valignments[(species, chain_type)][ seq ][22] != "C":
                sys.stderr.write(f"Warning - {chain_type},{species} doesn't "
                        "have CYS at position 23 and/or position 104.\n")

    return new_valignments



def fix_rhesus_heavy(sequence):
    """Rhesus heavy chains have gaps inserted in an inappropriate place in
    the IMGT data, which causes them to not follow the proper IMGT alignment.
    This problem is not unique to rhesus heavy chains and is in fact observed for
    some other chains (see below). For now we are fixing these manually
    since the problem lies with the IMGT data not with this script."""
    return sequence[:15] + sequence[16:27] + sequence[28:]

def fix_rhesus_kappa(sequence):
    """Rhesus kappa chains have insertions. This function fixes these insertions to preserve the
    alignment."""
    return sequence[:20] + sequence[21:]

def fix_rhesus_lambda(sequence):
    """Rhesus lambda chains have insertions. This function fixes these to preserve the alignment."""
    return sequence[:20]+sequence[21:51]+ sequence[53:]

def fix_mouse_delta(sequence):
    """Mouse delta chains have insertions. This function fixes these to preserve the alignment."""
    # Check in here because not all are bad...recheck again in the format v genes just to make sure.
    if sequence[103] != "C" or sequence[22] != "C":
        return sequence[ : 8 ] + sequence[ 9:85 ] + sequence[86:]
    return sequence

def fix_mouse_alpha(sequence):
    """Mouse alpha chains have insertions. This function fixes these to preserve the alignment."""
    return sequence[:8]+sequence[9:85]+sequence[86:]

def fix_rat_heavy(sequence):
    """Rat heavy chains have insertions. This function fixes these to preserve the alignment."""
    if sequence[-5] == ".":
        return sequence[:33] + sequence[34:-5] + sequence[-4:]
    return sequence[:33] + sequence[34:-6] + sequence[-5:]


def make_vj_alignments(vsequences, jsequences, fpath):
    """Generates all possible combinations of the V and J genes to form a consensus sequence."""
    for (species, chain_type) in vsequences:
        if (species, chain_type) not in jsequences:
            continue
        combined_sequences = {}
        for (vspecies, vallele), vsequence in vsequences[(species, chain_type)].items():
            for (_, jallele), jsequence in jsequences[(species, chain_type)].items():
                combined_sequences[f"{vspecies}_{vallele}_{jallele}".replace(" ", "_")] = \
                        vsequence + jsequence
        append_to_stockholm(combined_sequences,
                f"{hmbc.latin_to_common[species]}_{chain_type}",
                os.path.join(fpath, "ALL_ALIGNED.stockholm"), "a+")


def append_to_stockholm(sequences, idnum, outfile, mode = "a+"):
    """Appends to a stockholm file (or creates it if it does not yet exist)
    a new batch of sequences."""
    with open(outfile, mode, encoding="utf-8") as fhandle:
        fhandle.write("# STOCKHOLM 1.0\n")
        fhandle.write(f"#=GF ID {idnum}\n")

        pad_length = max(list(map(len, list(sequences.keys()))))+1
        for seq_id, seq in sequences.items():
            fhandle.write(f"{seq_id.replace(' ', '_').ljust(pad_length)}"
                    f" {seq.replace('.', '-')}\n")
        fhandle.write("#=GC RF".ljust(pad_length))
        fhandle.write(f"{''.join(['x' for s in seq])}\n")
        fhandle.write("//\n")





def build_stockholm_alignments(target_dir, selected_species, muscle_fpath):
    """Read in raw v and j sequence data from IMGT; combine it to form
    germline sequences, then use these to build stockholm alignment
    files for each chain. These stockholm alignment files are in turn
    used to form either HMMer profiles (this will likely be deprecated
    soon) and / or consensus alignments.
    """
    valignments, jalignments = {},{}
    print("Formatting IGs")
    for species in selected_species:
        for chain_type in ['H', 'K', 'L']:
            #Remove the + required in species names for IMGT
            species_c = species.replace("+", "_")
            if not os.path.isfile(os.path.join(target_dir, f"{species_c}_{chain_type}V.fasta")):
                continue

            print(f"Now working on {species_c}, {chain_type}")
            valignments[(species_c, chain_type)]  = read_fasta_input(os.path.join(target_dir,
                            f"{species_c}_{chain_type}V.fasta"), region_name = "V-REGION" )
            jalignments[(species_c, chain_type)]  = read_fasta_input(os.path.join(target_dir,
                            f"{species_c}_{chain_type}J.fasta"), region_name = "J-REGION")

    print("Formatting TRs")
    for species in selected_species:
        for chain_type in ['A', 'B', 'G', 'D']:
            #Remove the + required in species names for IMGT
            species_c = species.replace("+", "_")
            if not os.path.isfile( os.path.join(target_dir, f"{species_c}_{chain_type}V.fasta")):
                continue

            print(f"Now working on {species_c}, {chain_type}")
            valignments[(species_c, chain_type)]  = read_fasta_input(os.path.join(target_dir,
                            f"{species_c}_{chain_type}V.fasta"))
            jalignments[(species_c, chain_type)]  = read_fasta_input(os.path.join(target_dir,
                            f"{species_c}_{chain_type}J.fasta"))

    valignments = format_v_genes(valignments)
    jalignments = format_j_genes(jalignments, target_dir, muscle_fpath)
    make_vj_alignments(valignments, jalignments, target_dir)
