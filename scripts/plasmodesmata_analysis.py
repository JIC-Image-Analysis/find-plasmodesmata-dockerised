"""Analyse the plasmodesmata in a 3D image.

Segmentation of the image is done using an absolute threshold.

Small unwanted regions, such as background pixels, remaining from the
thresholding are filtered out based on a minimal requirement of a voxel size.

Large unwanted regions, such as stomata, remaining from the thresholding
are filtered out based on a maximum allowed voxel size.

"""

import os.path
import argparse
import warnings

import numpy as np

from jicbioimage.core.transform import transformation
from jicbioimage.core.io import DataManager, FileBackend, AutoName, AutoWrite
from jicbioimage.core.util.array import normalise
from jicbioimage.core.util.color import pretty_color
from jicbioimage.segment import SegmentedImage, connected_components
from jicbioimage.illustrate import AnnotatedImage


__version__ = "0.7.1"

AutoWrite.on = False
HERE = os.path.dirname(os.path.realpath(__file__))

# Suppress spurious scikit-image warnings.
warnings.filterwarnings("ignore", module="skimage.io._io")


def get_microscopy_collection(input_file):
    """Return microscopy collection from input file."""
    data_dir = os.path.abspath(os.path.join(HERE, "..", "data"))
    if not os.path.isdir(data_dir):
        os.mkdir(data_dir)
    backend_dir = os.path.join(data_dir, 'unpacked')
    file_backend = FileBackend(backend_dir)
    data_manager = DataManager(file_backend)
    microscopy_collection = data_manager.load(input_file)
    return microscopy_collection


@transformation
def threshold_abs(image, threshold):
    """Return thresholded image from an absolute cutoff."""
    return image > threshold


def segment3D(microscopy_collection, series, threshold):
    """Return segmented plasmodesmata in 3D."""
    zstack = microscopy_collection.zstack_array(s=series)
    segmentation = np.zeros(zstack.shape, dtype=bool)
    ziterator = microscopy_collection.zstack_proxy_iterator(s=series)
    for z, proxy_image in enumerate(ziterator):
        image = proxy_image.image
        image = threshold_abs(image, threshold)
        segmentation[:, :, z] = image
    return connected_components(segmentation, background=0)


def filter_large(segmentation3D, max_voxel):
    """Remove large regions."""
    removed = segmentation3D.copy()
    for i in segmentation3D.identifiers:
        region = segmentation3D.region_by_identifier(i)
        if region.area > max_voxel:
            segmentation3D[region] = 0
        else:
            removed[region] = 0
    return segmentation3D, removed


def filter_small(segmentation3D, min_voxel):
    """Remove small regions."""
    removed = segmentation3D.copy()
    for i in segmentation3D.identifiers:
        region = segmentation3D.region_by_identifier(i)
        if region.area < min_voxel:
            segmentation3D[region] = 0
        else:
            removed[region] = 0
    return segmentation3D, removed


@transformation
def annotate(image, segmentation):
    """Return annotated image."""
    uint8_normalised = normalise(image) * 255
    annotation = AnnotatedImage.from_grayscale(uint8_normalised)
    for i in segmentation.identifiers:
        region = segmentation.region_by_identifier(i)
        annotation.mask_region(region.dilate(1).border,
                               color=pretty_color(i))
    return annotation


def annotate3D(microscopy_collection, series, segmentation3D, name):
    for z in microscopy_collection.zslices(series):
        image = microscopy_collection.image(s=series, z=z)
        zslice = segmentation3D[:, :, z]
        segmentation = SegmentedImage.from_array(zslice)
        annotation = annotate(image, segmentation)
        fname = "z{:03d}_{}.png".format(z, name)
        fpath = os.path.join(AutoName.directory, fname)
        with open(fpath, "wb") as fh:
            fh.write(annotation.png())


def write_csv(segmentation3D, intensity, fname):
    """Write out a csv file with information about each spot."""
    header = ["id", "rgb", "voxels", "sum", "min", "max", "mean"]
    row = '{id:d},"{rgb}",{voxels:d},{sum:d},{min:d},{max:d},{mean:.3f}\n'
    with open(fname, "w") as fh:
        fh.write("{}\n".format(",".join(header)))
        for i in segmentation3D.identifiers:
            region = segmentation3D.region_by_identifier(i)
            values = intensity[region]
            data = dict(id=i, rgb=str(pretty_color(i)), voxels=region.area,
                        sum=int(np.sum(values)), min=int(np.min(values)),
                        max=int(np.max(values)), mean=float(np.mean(values)))
            fh.write(row.format(**data))


def plasmodesmata_analysis(microscopy_collection, series, threshold,
                           min_voxel, max_voxel):
    """Analyse the plasmodesmata in a 3D image.

    Segmentation of the image is done using an absolute threshold.

    Large unwanted regions, such as stomata, remaining from the thresholding
    are filtered out based on a maximum allowed voxel size.
    """
    segmentation = segment3D(microscopy_collection, series, threshold)

    # Filter out small and large regions.
    segmentation, small_removed = filter_small(segmentation, min_voxel)
    segmentation, large_removed = filter_large(segmentation, max_voxel)

    # Create annotated images.
    annotate3D(microscopy_collection, series, segmentation, "plasmodesmata")
    annotate3D(microscopy_collection, series, small_removed, "small_removed")
    annotate3D(microscopy_collection, series, large_removed, "large_removed")

    # Write out data to CSV files.
    csv_fn = os.path.join(AutoName.directory, "plasmodesmata.csv")
    write_csv(segmentation, microscopy_collection.zstack_array(s=series),
              csv_fn)
    csv_fn = os.path.join(AutoName.directory, "small.removed.csv")
    write_csv(small_removed, microscopy_collection.zstack_array(s=series),
              csv_fn)
    csv_fn = os.path.join(AutoName.directory, "large.removed.csv")
    write_csv(large_removed, microscopy_collection.zstack_array(s=series),
              csv_fn)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_file", help="path to raw microscopy data")
    parser.add_argument("series", type=int,
                        help="zero based index of series to analyse")
    parser.add_argument("output_dir", help="output directory")
    parser.add_argument("-t", "--threshold",
                        default=15000, type=int,
                        help="abs threshold (default=15000)")
    parser.add_argument("--min-voxel", default=2, type=int,
                        help="Minimum voxel volume (default=2)")
    parser.add_argument("--max-voxel", default=50, type=int,
                        help="Maximum voxel volume (default=50)")
    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        parser.error("No such file: {}".format(args.input_file))

    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)
    AutoName.directory = args.output_dir

    microscopy_collection = get_microscopy_collection(args.input_file)
    plasmodesmata_analysis(microscopy_collection, args.series, args.threshold,
                           args.min_voxel, args.max_voxel)


if __name__ == "__main__":
    main()
