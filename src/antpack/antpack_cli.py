"""Contains tools needed to run the command line interface."""
import os
import sys
import argparse
from antpack import SingleChainAnnotator, PairedChainAnnotator, VJGeneTool

from .cli_tools.biofile_read_tools import read_fasta


class ReconfigParser(argparse.ArgumentParser):
    """Reconfigure argparse's parser so an automatic error message
    is generated if no args supplied."""
    def error(self, message):
        self.print_help()
        sys.exit(2)


def gen_arg_parser():
    """Build the command line arg parser."""
    parser = ReconfigParser(description="Run AntPack from the command line "
            "to number and assign VJ genes for input mAbs and TCRs. "
            "This command line tool writes its output to a csv or fasta. "
            "It is appropriate for a quick generic analysis on datasets from a couple to "
            "a few hundred thousand sequences. If you need to work with a larger "
            "dataset or do a more customized analysis, use the Python API instead, which "
            "is a little faster and much more flexible.")
    parser.add_argument("input", nargs = 1, help=
            "Input filepath. Must be a fasta file. It can be gzipped -- AntPack "
            "will assume it is gzipped if the path ends in .gz.")
    parser.add_argument("output", nargs = 1, help=
            "Output filepath. Two output files are created, one with '_heavy' "
            "appended for heavy chains and one with '_light' appended for light "
            "chains.")
    parser.add_argument("species", nargs = 1, help=
            "Species. This is used to find the closest matching germline genes (by percent "
            "identity). If multiple genes have the same percent identity they are "
            "all reported. Valid options are human, rabbit, alpaca, mouse, unknown. "
            "For TCRs, valid options are human, mouse, unknown. If fasta output is "
            "selected, no VJ gene id is done and this argument is ignored.")
    parser.add_argument("scheme", nargs = 1, help=
            "One of aho, imgt, kabat or martin. For "
            "TCRs, only imgt is accepted.")
    parser.add_argument("--paired", action="store_true", help=
            "This tool will normally assume there is one variable region "
            "per input sequence. If this flag is supplied, it will instead "
            "assume each input may contain a heavy chain, a light chain, or "
            "both and try to extract whatever is present. Use this if at "
            "least some of your sequences are paired chains.")
    parser.add_argument("--tcrs", action="store_true", help=
            "This tool will normally assume input sequences are mAbs. "
            "If this flag is supplied it will instead assume all input "
            "sequences are TCRs.")
    parser.add_argument("--evalue", action="store_true", help=
            "If this flag is supplied, VJ genes are assigned by evalue "
            "rather than by identity.")
    parser.add_argument("--fasta", action="store_true", help=
            "Output is normally written to a csv. If this flag is supplied, "
            "output is instead written to a fasta file and no VJ gene matching "
            "is done. This is faster but also less informative (no error messages "
            "are output into the fasta, no VJ genes are assigned etc.")
    parser.add_argument("--offline", action="store_true", help=
            "This tool will normally load all sequences to memory before "
            "processing them. This is faster but undesirable if there are "
            "a large number of sequences. If this flag is supplied, temporary "
            "files are created in your current working directory and processing "
            "is performed using on-disk storage. This is slower but has minimal "
            "memory usage.")
    return parser






def process_fasta_online(cli_args):
    """Processes a fasta file by first loading all sequences
    into memory. This is faster and simpler but obviously
    greatly increases memory consumption if the number of
    sequences in the input file is large."""
    vj_tool = VJGeneTool(scheme = cli_args.scheme[0])
    vj_mode = "identity"
    if cli_args.evalue:
        vj_mode = "evalue"

    if cli_args.tcrs:
        if cli_args.species[0] not in ["human", "mouse", "unknown"]:
            raise RuntimeError("Invalid species supplied for TCRs; "
                    "please check directions.")
        if cli_args.paired:
            sc_tool = PairedChainAnnotator(receptor_type="tcr",
                    scheme=cli_args.scheme[0])
        else:
            sc_tool = SingleChainAnnotator(chains=["A", "B", "D", "G"],
                    scheme=cli_args.scheme[0])
    else:
        if cli_args.species[0] not in ["human", "alpaca", "rabbit",
                "mouse", "unknown"]:
            raise RuntimeError("Invalid species supplied; "
                    "please check directions.")
        if cli_args.paired:
            sc_tool = PairedChainAnnotator(scheme=cli_args.scheme[0])
        else:
            sc_tool = SingleChainAnnotator(scheme=cli_args.scheme[0])

    seqrecs = list(read_fasta(cli_args.input[0]))
    seqs = [seqrec[1] for seqrec in seqrecs]
    seqinfo = [seqrec[0] for seqrec in seqrecs]
    seq_annotations = [sc_tool.analyze_seq(s) for s in seqs]

    output_dict = {
        k:{"annotations":[], "seqs":[],
        "seqinfo":[], "v_genes":[], "j_genes":[],
        "v_scores":[], "j_scores":[], "vj_species":[],
        "msa":None}
        for k in ["heavy", "light"]
        }

    if cli_args.paired:
        for i, chain in enumerate(["heavy", "light"]):
            output_dict[chain]["annotations"] = [s[i] for s in seq_annotations]
            output_dict[chain]["seqs"] = seqs
            output_dict[chain]["seqinfo"] = seqinfo
            output_dict[chain]["msa"] = sc_tool.build_msa(seqs,
                    output_dict[chain]["annotations"], True)

        del seq_annotations

    else:
        for i, (seq, annot) in enumerate(zip(seqs, seq_annotations)):
            if annot[2] not in ("H", "A", "G"):
                output_dict["light"]["annotations"].append(annot)
                output_dict["light"]["seqs"].append(seq)
                output_dict["light"]["seqinfo"].append(seqinfo[i])
            else:
                output_dict["heavy"]["annotations"].append(annot)
                output_dict["heavy"]["seqs"].append(seq)
                output_dict["heavy"]["seqinfo"].append(seqinfo[i])

        for chain in ["heavy", "light"]:
            if len(output_dict[chain]["seqs"]) > 0:
                output_dict[chain]["msa"] = sc_tool.build_msa(output_dict[chain]["seqs"],
                        output_dict[chain]["annotations"], True)

    if not cli_args.fasta:
        for chain in ["heavy", "light"]:
            for seq, annot in zip(output_dict[chain]["seqs"], output_dict[chain]["annotations"]):
                v_gene, j_gene, v_ident, j_ident, vj_species = \
                        vj_tool.assign_vj_genes(annot, seq,
                        cli_args.species[0], vj_mode)
                output_dict[chain]["vj_species"].append(vj_species)
                output_dict[chain]["v_genes"].append(v_gene)
                output_dict[chain]["j_genes"].append(j_gene)
                output_dict[chain]["v_scores"].append(v_ident)
                output_dict[chain]["j_scores"].append(j_ident)

        for chain, cdict in output_dict.items():
            if cdict["msa"] is None:
                continue
            with open(cli_args.output[0] + f"_{chain}.csv", "w+",
                    encoding="utf-8") as fhandle:
                _ = fhandle.write("Sequence_info,percent_identity,vj_species,"
                            "vj_mode,v_gene,v_score,j_gene,j_score,")
                _ = fhandle.write(f"{','.join(cdict['msa'][0])},error_message\n")

                for i, (seq_id, annotation, msa_row) in enumerate(zip(cdict["seqinfo"],
                    cdict["annotations"], cdict["msa"][1])):
                    percent_identity = 100 * round(annotation[1], 3)
                    _ = fhandle.write(f"{seq_id},{percent_identity},")
                    _ = fhandle.write(",".join([cdict["vj_species"][i], vj_mode,
                        cdict["v_genes"][i], str(cdict["v_scores"][i]),
                        cdict["j_genes"][i], str(cdict["j_scores"][i]) ] + list(msa_row) ))
                    _ = fhandle.write(f",{cdict['annotations'][i][3]}\n")

    else:
        for chain, cdict in output_dict.items():
            if cdict["msa"] is None:
                continue
            with open(cli_args.output[0] + f"_{chain}.fasta", "w+",
                    encoding="utf-8") as fhandle:
                for seq_id, msa_row in zip(cdict["seqinfo"], cdict["msa"][1]):
                    _ = fhandle.write(f">{seq_id}\n{''.join(msa_row)}\n")



def process_fasta_offline(cli_args):
    """Processes a fasta file by writing results to a temporary
    file. This uses minimal memory but obviously is slower so
    is only preferable if the number of sequences in the input
    file is large."""
    vj_tool = VJGeneTool(scheme = cli_args.scheme[0])
    vj_mode = "identity"
    if cli_args.evalue:
        vj_mode = "evalue"

    if cli_args.tcrs:
        if cli_args.species[0] not in ["human", "mouse", "unknown"]:
            raise RuntimeError("Invalid species supplied for TCRs; "
                    "please check directions.")
        if cli_args.paired:
            sc_tool = PairedChainAnnotator(receptor_type="tcr",
                    scheme=cli_args.scheme[0])
        else:
            sc_tool = SingleChainAnnotator(chains=["A", "B", "D", "G"],
                    scheme=cli_args.scheme[0])
    else:
        if cli_args.species[0] not in ["human", "alpaca", "rabbit",
                "mouse", "unknown"]:
            raise RuntimeError("Invalid species supplied; "
                    "please check directions.")
        if cli_args.paired:
            sc_tool = PairedChainAnnotator(scheme=cli_args.scheme[0])
        else:
            sc_tool = SingleChainAnnotator(scheme=cli_args.scheme[0])

    try:
        with open("ANTPACK_light_TEMP_FILE.txt", "w+",
                encoding="utf-8") as fhandle:
            pass
        os.remove("ANTPACK_light_TEMP_FILE.txt")
    except:
        raise RuntimeError("Could not open a temporary file in "
                "current working directory!")

    heavy_codes, light_codes = set(), set()

    with open("ANTPACK_heavy_TEMP_FILE.txt", "w+",
            encoding="utf-8") as heavy_handle:
        with open("ANTPACK_light_TEMP_FILE.txt", "w+",
                encoding="utf-8") as light_handle:
            if cli_args.paired:
                for seqrec in read_fasta(cli_args.input[0]):
                    seqinfo, seq = seqrec[0], seqrec[1]
                    heavy_ann, light_ann = sc_tool.analyze_seq(seq)
                    heavy_codes.update(heavy_ann[0])
                    light_codes.update(light_ann[0])

                    trimmed_h_seq, trimmed_h_align, _, _ = sc_tool.trim_alignment(
                                seq, heavy_ann)
                    trimmed_l_seq, trimmed_l_align, _, _ = sc_tool.trim_alignment(
                                seq, light_ann)

                    hv_gene, hj_gene, hv_ident, hj_ident = "", "", "", ""
                    lv_gene, lj_gene, lv_ident, lj_ident = "", "", "", ""
                    hvj_species, lvj_species = "", ""

                    if not cli_args.fasta:
                        hv_gene, hj_gene, hv_ident, hj_ident, hvj_species = \
                                vj_tool.assign_vj_genes(heavy_ann, seq,
                                cli_args.species[0], vj_mode)
                        lv_gene, lj_gene, lv_ident, lj_ident, lvj_species = \
                                vj_tool.assign_vj_genes(light_ann, seq,
                                cli_args.species[0], vj_mode)

                    heavy_handle.write(f"{seqinfo},{heavy_ann[1]},{hvj_species},"
                                f"{vj_mode},{hv_gene},{hv_ident},{hj_gene},{hj_ident},"
                                f"\t\t\t\t{'_'.join(trimmed_h_align)},"
                                f"{trimmed_h_seq},{heavy_ann[-1]}\n")
                    light_handle.write(f"{seqinfo},{light_ann[1]},{lvj_species},"
                                f"{vj_mode},{lv_gene},{lv_ident},{lj_gene},{lj_ident},"
                                f"\t\t\t\t{'_'.join(trimmed_l_align)},"
                                f"{trimmed_l_seq},{light_ann[-1]}\n")

            else:
                for seqrec in read_fasta(cli_args.input[0]):
                    seqinfo, seq = seqrec[0], seqrec[1]
                    annotation = sc_tool.analyze_seq(seq)
                    if annotation[2] in ('H', 'A', 'G'):
                        output_handle = heavy_handle
                        heavy_codes.update(annotation[0])
                    else:
                        output_handle = light_handle
                        light_codes.update(annotation[0])

                    trimmed_seq, trimmed_align, _, _ = sc_tool.trim_alignment(
                                seq, annotation)
                    v_gene, j_gene, v_ident, j_ident = "", "", "", ""
                    vj_species = ""

                    if not cli_args.fasta:
                        v_gene, j_gene, v_ident, j_ident, vj_species = \
                                vj_tool.assign_vj_genes(annotation, seq,
                                cli_args.species[0], vj_mode)

                    output_handle.write(f"{seqinfo},{annotation[1]},{vj_species},"
                                f"{vj_mode},{v_gene},{v_ident},{j_gene},{j_ident},"
                                f"\t\t\t\t{'_'.join(trimmed_align)},"
                                f"{trimmed_seq},{annotation[-1]}\n")


    for i in range(1,129):
        if str(i) not in heavy_codes:
            heavy_codes.add(str(i))
        if str(i) not in light_codes:
            light_codes.add(str(i))

    all_codes = {"heavy":sc_tool.sort_position_codes(list(heavy_codes)),
            "light":sc_tool.sort_position_codes(list(light_codes))}
    code_key = {chain:{k:i for i,k in enumerate(codes)}
            for chain, codes in all_codes.items()}

    for chain, chain_key in code_key.items():
        if len(all_codes[chain]) == 0:
            continue

        with open(f"ANTPACK_{chain}_TEMP_FILE.txt", "r",
                encoding="utf-8") as input_handle:

            if not cli_args.fasta:
                with open(cli_args.output[0] + f"_{chain}.csv", "w+",
                    encoding="utf-8") as output_handle:
                    _ = output_handle.write("Sequence_info,percent_identity,vj_species,"
                            "vj_mode,v_gene,v_score,j_gene,j_score,")
                    _ = output_handle.write(f"{','.join(all_codes[chain])},error_message\n")

                    for row in input_handle:
                        output_handle.write(row.split("\t\t\t\t")[0])
                        gapped_sequence = ['-' for c in chain_key]
                        numbering, seq, err = row.strip().split("\t\t\t\t")[1].split(',')
                        for nmbr, letter in zip(numbering.split("_"), seq):
                            if nmbr == '-':
                                continue
                            try:
                                gapped_sequence[chain_key[nmbr]] = letter
                            except:
                                import pdb
                                pdb.set_trace()

                        output_handle.write(f"{','.join(gapped_sequence)},{err}\n")


            else:
                with open(cli_args.output[0] + f"_{chain}.fasta", "w+",
                    encoding="utf-8") as output_handle:
                    for row in input_handle:
                        output_handle.write(f">{row.split(',')[0]}\n")
                        gapped_sequence = ['-' for c in chain_key]
                        numbering, seq, err = row.strip().split("\t\t\t\t")[1].split(',')
                        for nmbr, letter in zip(numbering.split("_"), seq):
                            if nmbr == '-':
                                continue
                            gapped_sequence[chain_key[nmbr]] = letter

                        output_handle.write(f"{''.join(gapped_sequence)}\n")

        os.remove(f"ANTPACK_{chain}_TEMP_FILE.txt")



def run_cli_interface():
    """Runs the command line interface using the user's specified
    options."""
    parser = gen_arg_parser()
    if len(sys.argv) < 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    cli_args = parser.parse_args()

    if cli_args.offline:
        process_fasta_offline(cli_args)
    else:
        process_fasta_online(cli_args)
