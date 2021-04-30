"""pytrapment main module create entrapment dbs."""

import argparse
import os
import sys
import time

from pytrapment import __version__ as xv
from pytrapment import entrapment


def arg_parser():  # pragma: not covered
    """
    Parse the arguments from the CLI.

    Returns:
        arguments, from parse_args
    """
    description = """
    pytrapment is a convenient tool to create entrapment databases for protemic MS analysis.
    Entrapment databases are created by sampling for each source fasta a protein from the
    entrapment database. This is done while conserving some important properties of the proteins
    in the source file.
    
    Use --help to see the command line arguments.

    Current Version: {}
    """.format(xv)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--in_fasta",
                        help="Input protein fasta file.",
                        required=True, action="store", dest="in_fasta")

    parser.add_argument("-e", "--entrapment",
                        help="Entrapment proteins.",
                        required=True, action="store", dest="entrapment")

    parser.add_argument("-o", "--out_dir",
                        help="Directory to store the results",
                        required=True, action="store", dest="out_dir")
    return parser


def pytrapment_runner(fasta_host, fasta_trap):
    """
    Execute pytrapment.

    Args:
        in_fasta: str, location of the input fasta file.
        entrapment_fasta: str, location of the input entrapment file.

    Returns:
        None
    """
    start_time = time.time()
    entrapment.get_nearest_neighbor_proteins(fasta_host, fasta_trap)
    end_time = time.time()
    print(f"Took {(end_time-start_time)/60.} minutes")


if __name__ == "__main__":  # pragma: no cover
    """Run pytrapment main function."""
    parser = arg_parser()
    try:
        args = parser.parse_args(sys.argv[1:])
    except TypeError:
        parser.print_usage()

    # create dir if not there
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    # call function
    pytrapment_runner(args.in_fasta, args.entrapment, args.outdir)
