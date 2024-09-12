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
            "Input file. Must be in fasta format.")
    return parser


def run_cli_interface():
    """Runs the command line interface using the user's specified
    options."""
    parser = gen_arg_parser()
    if len(sys.argv) < 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    print(args.input)
