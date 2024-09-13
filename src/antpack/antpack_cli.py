"""Contains tools needed to run the command line interface."""
import os
import sys
import argparse



class ReconfigParser(argparse.ArgumentParser):
    """Reconfigure argparse's parser so an automatic error message
    is generated if no args supplied."""
    def error(self, message):
        self.print_help()
        sys.exit(2)


def gen_arg_parser():
    """Build the command line arg parser."""
    parser = ReconfigParser(description="Run AntPack from the command line "
            "to number and provide VJ assignments for input AA sequences.")
    parser.add_argument("input", nargs = 1, help=
            "Input filepath. Must be in fasta format.")
    parser.add_argument("output", nargs = 1, help=
            "Output filepath. No file extension is needed since "
            "an appropriate extension will be added depending on "
            "additional options.")
    parser.add_argument("--csv", action="store_true", help=
            "The output is normally provided in fasta format. If this "
            "flag is supplied, the output will be supplied as a csv "
            "with one chain per row.")
    parser.add_argument("--paired", action="store_true", help=
            "AntPack will normally assume there is one variable region "
            "per input sequence and try to extract it. If this "
            "flag is supplied, it will instead assume each input "
            "sequence contains a heavy chain and a light chain.")
    parser.add_argument("--vj", nargs=1, help=
            "AntPack will not by default assign VJ genes. If this "
            "flag is supplied, however, it will find the closest "
            "VJ genes for the specified species and write these "
            "into the output. The specified species must be one of "
            "human, mouse. This flag is ONLY accepted for csv "
            "output; if the --csv flag is not supplied this flag "
            "is ignored.")
    parser.add_argument("--low_memory", action="store_true", help=
            "AntPack will normally read all of the sequences from "
            "your input into memory before analyzing them, which is "
            "fine if the number of sequences is small. If it is large "
            "supply this flag, which causes AntPack to read only one "
            "sequence into memory at a time and save results into a "
            "temporary file. This approach is slower "
            "but has minimal memory consumption.")
    return parser


def run_cli_interface():
    """Runs the command line interface using the user's specified
    options."""
    parser = gen_arg_parser()
    if len(sys.argv) < 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
