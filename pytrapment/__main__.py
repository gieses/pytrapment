"""xiRT main module to run the training and prediction."""

import argparse
import logging
import os
import pickle
import sys
import time
from datetime import datetime
from pytrapment import __version__ as xv
logger = logging.getLogger(__name__)


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


def pytrapment_runner(in_fasta, entrapment_fasta):
    """
    Execute pytrapment.

    Args:
        in_fasta: str, location of the input fasta file.
        entrapment_fasta: str, location of the input entrapment file.

    Returns:
        None
    """
    start_time = time.time()
    logger.info("pytrapment fasta: {}".format(in_fasta))
    logger.info("pytrapment entrapment: {}".format(entrapment_fasta))
    pass

def main():  # pragma: no cover
    """Run xiRT main function."""
    parser = arg_parser()
    try:
        args = parser.parse_args(sys.argv[1:])
    except TypeError:
        parser.print_usage()

    # create dir if not there
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    # create logger
    logger = logging.getLogger('xirt')
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to debug
    ch = logging.FileHandler(os.path.join(args.out_dir, "pytrapment_logger.log"), "w")
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(ch)

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.DEBUG)
    sh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(sh)

    logger.info("command line call:")
    logger.info("pytrapment -i {} -e {} -o {}".format(args.in_fasta, args.entrapment,
                                                            args.outdir))
    logger.info("Init logging file.")
    logger.info("Starting Time: {}".format(datetime.now().strftime("%H:%M:%S")))
    logger.info("Starting xiRT.")
    logger.info("Using xiRT version: {}".format(xv))

    # call function
    pytrapment_runner(args.in_fasta, args.entrapment, args.outdir)


if __name__ == "__main__":  # pragma: no cover
    main()
