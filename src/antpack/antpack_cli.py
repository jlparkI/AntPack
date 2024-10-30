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
            "to number and provide VJ assignments for input AA sequences. "
            "This command line tool first reads all input sequences into memory, "
            "and writes them to a csv in a format that makes them easy to review "
            "and compare (appropriate for small-moderate datasets). If you need to work "
            "with a larger dataset or do a more customized or extensive analysis, "
            "use the Python API instead.")
    parser.add_argument("input", nargs = 1, help=
            "Input filepath. Must be in fasta format.")
    parser.add_argument("output", nargs = 1, help=
            "Output filepath. Two output files are created, one output "
            "file with '_heavy.csv' appended to the output path for "
            "heavy chains (if any) and one with '_light.csv' appended "
            "to the output path for light chains (if any).")
    parser.add_argument("scheme", nargs = 1, help=
            "The numbering scheme. One of aho, imgt, kabat or martin.")
    parser.add_argument("--paired", action="store_true", help=
            "AntPack will normally assume there is one variable region "
            "per input sequence and try to extract it. If this "
            "flag is supplied, it will instead assume each input "
            "sequence contains a heavy chain and a light chain and "
            "try to extract both. The --paired flag can still be "
            "supplied if some of the sequences are single-chain.")
    parser.add_argument("--vj", nargs=2, help=
            "If this flag is supplied, AntPack will find the closest "
            "VJ genes for the specified species and write these "
            "into the output. There are two arguments. The first is "
            "the species which must be either human or mouse. The "
            "second is the mode which must be identity or evalue. ",
            metavar=("species", "mode"))
    parser.add_argument("--chains", nargs=1, help=
            "AntPack will normally look for an H, K or L chain in each "
            "input sequence. If desired, you can instead restrict it to "
            "search for a comma-separated list of specific chains (e.g. H or "
            "K,L) by using this argument. If you supply the argument "
            "--paired this argument is ignored.",
            metavar=("chains"))
    return parser






def process_fasta_online(cli_args, chains):
    """Processes a fasta file by first loading all sequences
    into memory. This is faster and simpler but obviously
    greatly increases memory consumption if the number of
    sequences in the input file is large. For a large
    file, it is strongly recommended to use the Python
    API instead."""
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

    args = parser.parse_args()
    if args.chains:
        chains = args.chains.split(",")
    else:
        chains = ['H', 'K', 'L']

    process_fasta_online(args, chains)
