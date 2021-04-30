"""pytrapment main module create entrapment dbs."""

import argparse
import os
import sys
import time
from datetime import date

from pyteomics import fasta

from pytrapment import __version__ as xv
from pytrapment import entrapment, qc


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


def main():  # pragma: no cover
    """
    Execute pytrapment.

    Returns:
        None
    """
    today = date.today().strftime("%Y%m%d")

    parser = arg_parser()
    try:
        args = parser.parse_args(sys.argv[1:])
    except TypeError:
        parser.print_usage()

    # create dir if not there
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    start_time = time.time()
    fasta_df = entrapment.get_nearest_neighbor_proteins(args.fasta_host, args.fasta_trap)
    _ = fasta.write(zip(fasta_df.index, fasta_df.sequence),
                    os.path.join(args.out_dir, f"entrapment_{today}.fasta"),
                    file_mode="w")

    end_time = time.time()
    print(f"Took {(end_time-start_time)/60.} minutes")

    # doing qc
    host_peptides = entrapment.digest_protein_df(fasta_df[fasta_df["db_type"] == "host"])
    trap_peptides = entrapment.digest_protein_df(fasta_df[fasta_df["db_type"] == "trap"])

    features_df_host = qc.compute_sequence_features(host_peptides)
    features_df_host["Type"] = "host"

    features_df_trap = qc.compute_sequence_features(trap_peptides)
    features_df_trap["Type"] = "trap"

    qc.qc_peptides(features_df_host, features_df_trap, args.out_dir)


if __name__ == "__main__":  # pragma: no cover
    """Run pytrapment main function."""
    main()
