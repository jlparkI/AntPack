"""Contains tools needed to run the command line interface."""
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
            "is slightly faster and much more flexible.")
    parser.add_argument("input", nargs = 1, help=
            "Input filepath. Must be a fasta file. The location of this file "
            "will also be used for temporary files created during the sequence "
            "processing.")
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
    parser.add_argument("--fasta", action="store_true", help=
            "Output is normally written to a csv. If this flag is supplied, "
            "output is instead written to a fasta file and no VJ gene matching "
            "is done. This is faster but also less informative (no error messages "
            "are output into the fasta, no VJ genes are assigned etc.")
    return parser






def process_fasta_online(cli_args):
    """Processes a fasta file by first loading all sequences
    into memory. This is faster and simpler but obviously
    greatly increases memory consumption if the number of
    sequences in the input file is large."""
    vj_tool = VJGeneTool(scheme = cli_args.scheme[0])

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
        "v_scores":[], "j_scores":[], "msa":None}
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
                v_gene, j_gene, v_ident, j_ident = vj_tool.assign_vj_genes(annot, seq,
                        cli_args.species[0], "identity")
                output_dict[chain]["v_genes"].append(v_gene)
                output_dict[chain]["j_genes"].append(j_gene)
                output_dict[chain]["v_scores"].append(v_ident)
                output_dict[chain]["j_scores"].append(j_ident)

    if not cli_args.fasta:
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
                    _ = fhandle.write(",".join([cli_args.vj[0], cli_args.vj[1],
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


def process_fasta_offline(cli_args, chains):
    """Processes a fasta file by loading sequences into
    memory one at a time. This is slower than processing online
    but avoids memory issues if the dataset is large."""
    seqrecs = list(read_fasta(cli_args.input[0]))
    seqs = [seqrec[1] for seqrec in seqrecs]
    seqinfo = [seqrec[0] for seqrec in seqrecs]
    vj_tool = VJGeneTool(scheme = cli_args.scheme[0])

    output_dict = {
        k:{"annotations":[], "seqs":[],
        "seqinfo":[], "v_genes":[], "j_genes":[],
        "v_scores":[], "j_scores":[], "msa":None}
        for k in ["heavy", "light"]
        }

    if cli_args.paired:
        sc_tool = PairedChainAnnotator(cli_args.scheme[0])
        seq_annotations = [sc_tool.analyze_seq(s) for s in seqs]
        for i, chain in enumerate(["heavy", "light"]):
            output_dict[chain]["annotations"] = [s[i] for s in seq_annotations]
            output_dict[chain]["seqs"] = seqs
            output_dict[chain]["seqinfo"] = seqinfo
            output_dict[chain]["msa"] = sc_tool.build_msa(seqs,
                    output_dict[chain]["annotations"], True)

        del seq_annotations

    else:
        sc_tool = SingleChainAnnotator(chains, cli_args.scheme[0])
        annotations = sc_tool.analyze_seqs(seqs)

        for i, (seq, annot) in enumerate(zip(seqs, annotations)):
            if annot[2] in ("K", "L"):
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

    if cli_args.vj is not None:
        for chain in ["heavy", "light"]:
            for seq, annot in zip(output_dict[chain]["seqs"], output_dict[chain]["annotations"]):
                v_gene, j_gene, v_ident, j_ident = vj_tool.assign_vj_genes(annot, seq,
                        cli_args.vj[0], cli_args.vj[1])
                output_dict[chain]["v_genes"].append(v_gene)
                output_dict[chain]["j_genes"].append(j_gene)
                output_dict[chain]["v_scores"].append(v_ident)
                output_dict[chain]["j_scores"].append(j_ident)

    for chain, cdict in output_dict.items():
        if cdict["msa"] is None:
            continue
        with open(cli_args.output[0] + f"_{chain}.csv", "w+",
                encoding="utf-8") as fhandle:
            if cli_args.vj is not None:
                _ = fhandle.write("Sequence_info,percent_identity,vj_species,"
                        "vj_mode,v_gene,v_score,j_gene,j_score,")
            else:
                _ = fhandle.write("Sequence_info,percent_identity,")
            _ = fhandle.write(f"{','.join(cdict['msa'][0])},error_message\n")

            if not cli_args.vj:
                for i, (seq_id, annotation, msa_row) in enumerate(zip(cdict["seqinfo"],
                    cdict["annotations"], cdict["msa"][1])):
                    percent_identity = 100 * round(annotation[1], 3)
                    _ = fhandle.write(f"{seq_id},{percent_identity},")
                    _ = fhandle.write(",".join(list(msa_row)))
                    _ = fhandle.write(f",{cdict['annotations'][i][3]}\n")

            else:
                for i, (seq_id, annotation, msa_row) in enumerate(zip(cdict["seqinfo"],
                    cdict["annotations"], cdict["msa"][1])):
                    percent_identity = 100 * round(annotation[1], 3)
                    _ = fhandle.write(f"{seq_id},{percent_identity},")
                    _ = fhandle.write(",".join([cli_args.vj[0], cli_args.vj[1],
                        cdict["v_genes"][i], str(cdict["v_scores"][i]),
                        cdict["j_genes"][i], str(cdict["j_scores"][i]) ] +
                        list(msa_row) ))
                    _ = fhandle.write(f",{cdict['annotations'][i][3]}\n")





def run_cli_interface():
    """Runs the command line interface using the user's specified
    options."""
    parser = gen_arg_parser()
    if len(sys.argv) < 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    process_fasta_online(parser.parse_args())
