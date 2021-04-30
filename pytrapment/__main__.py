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
    Entrapment databases are created by sampling for each host fasta protein a similar protein
    from the entrapment database.

    Use --help to see the command line arguments.

    Current Version: {}
    """.format(xv)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--fasta_host",
                        help="Input protein fasta file.",
                        required=True, action="store", dest="fasta_host")

    parser.add_argument("-t", "--fasta_trap",
                        help="Entrapment proteins.",
                        required=True, action="store", dest="fasta_trap")

    parser.add_argument("-o", "--out_dir",
                        help="Directory to store the results",
                        required=True, action="store", dest="out_dir")
    return parser


def main():
    """
    Execute pytrapment.

    Returns:
        None
    """
    """Run pytrapment main function."""
    parser = arg_parser()
    try:
        args = parser.parse_args(sys.argv[1:])
    except TypeError:
        parser.print_usage()

    # create dir if not there
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    start_time = time.time()
    entrapment.get_nearest_neighbor_proteins(args.fasta_host, args.fasta_trap)
    end_time = time.time()
    print(f"Took {(end_time-start_time)/60.} minutes")


if __name__ == "__main__":  # pragma: no cover
    """Run pytrapment main function."""
    main()
