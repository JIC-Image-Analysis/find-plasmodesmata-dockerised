"""Analyse all images in an input directory."""

import os
import os.path
import argparse
import logging

from plasmodesmata_analysis import (
    get_microscopy_collection,
    __version__,
    log_settings,
)
from analyse_image import analyse_all_series

# Setup logging with a stream handler.
logger = logging.getLogger(os.path.basename(__file__))
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
logger.addHandler(ch)


def analyse_dir(args):
    """Analyse all images in an input directory."""
    for fname in os.listdir(args.input_dir):
        if not fname.lower().endswith(".lif"):
            continue
        logger.info("Analysing image: {}".format(fname))

        def get_dir_name(fname):
            no_suffix_list = fname.split(".")[0:-1]
            return ".".join(no_suffix_list)
        dir_name = get_dir_name(fname)
        specific_out_dir = os.path.join(args.output_dir, dir_name)

        # Skip analysis of image if output directory exists.
        if os.path.isdir(specific_out_dir):
            logger.info("Directory exists: {}".format(specific_out_dir))
            logger.info("Skipping: {}".format(fname))
            continue
        os.mkdir(specific_out_dir)

        fpath = os.path.join(args.input_dir, fname)
        microscopy_collection = get_microscopy_collection(fpath)
        analyse_all_series(microscopy_collection, specific_out_dir, args.threshold,
                           args.min_voxel, args.max_voxel)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("input_dir", help="input directory")
    parser.add_argument("output_dir", help="output directory")
    parser.add_argument("-t", "--threshold",
                        default=15000, type=int,
                        help="abs threshold (default=15000)")
    parser.add_argument("--min-voxel", default=2, type=int,
                        help="Minimum voxel volume (default=2)")
    parser.add_argument("--max-voxel", default=50, type=int,
                        help="Maximum voxel volume (default=50)")
    args = parser.parse_args()

    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)

    # Create file handle logger.
    fh = logging.FileHandler(os.path.join(args.output_dir, "log"), mode="w")
    fh.setLevel(logging.DEBUG)
    format_ = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    formatter = logging.Formatter(format_)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    log_settings(__version__, args)

    analyse_dir(args)

if __name__ == "__main__":
    main()
