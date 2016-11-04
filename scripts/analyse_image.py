"""Analyse all series."""

import os
import os.path
import argparse
import logging

from jicbioimage.core.io import AutoName

from plasmodesmata_analysis import (
    get_microscopy_collection,
    plasmodesmata_analysis,
    __version__,
    log_settings,
)


# Setup logging with a stream handler.
logger = logging.getLogger(os.path.basename(__file__))
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
logger.addHandler(ch)


def analyse_all_series(microscopy_collection, output_dir, threshold, min_voxel,
                       max_voxel):
    """Analyse all series in input microscopy file."""
    for s in microscopy_collection.series:
        sub_dir = os.path.join(output_dir, str(s))
        if not os.path.isdir(sub_dir):
            os.mkdir(sub_dir)

        AutoName.directory = sub_dir

        logger.info("Analysing series: {}".format(s))
        plasmodesmata_analysis(microscopy_collection, s, threshold, min_voxel,
                               max_voxel)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("input_file", help="path to raw microscopy data")
    parser.add_argument("output_dir", help="output directory")
    parser.add_argument("-t", "--threshold",
                        default=15000, type=int,
                        help="abs threshold (default=15000)")
    parser.add_argument("--min-voxel", default=2, type=int,
                        help="Minimum voxel volume (default=2)")
    parser.add_argument("--max-voxel", default=50, type=int,
                        help="Maximum voxel volume (default=50)")
    args = parser.parse_args()

    dir_name = os.path.basename(args.input_file).split(".")[0]
    specific_out_dir = os.path.join(args.output_dir, dir_name)

    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)
    if not os.path.isdir(specific_out_dir):
        os.mkdir(specific_out_dir)
    if not os.path.isfile(args.input_file):
        parser.error("No such microscopy file: {}".format(args.input_file))

    # Create file handle logger.
    fh = logging.FileHandler(os.path.join(specific_out_dir, "log"), mode="w")
    fh.setLevel(logging.DEBUG)
    format_ = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    formatter = logging.Formatter(format_)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    log_settings(logger, __version__, args)

    microscopy_collection = get_microscopy_collection(args.input_file)
    analyse_all_series(microscopy_collection, specific_out_dir, args.threshold,
                       args.min_voxel, args.max_voxel)

if __name__ == "__main__":
    main()
