"""Functional test."""

import os
import os.path

from jicbioimage.core.image import MicroscopyCollection
from jicbioimage.core.io import (
    AutoName,
    DataManager,
    FileBackend,
    _md5_hexdigest_from_file,
)
from plasmodesmata_analysis import plasmodesmata_analysis


def test_plasmodesmata_analysis():
    output_dir = "/output/tmp"
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    try:
        AutoName.directory = output_dir

        backend_dir = "/backend"
        test_image = "/data/test.tif"
        file_backend = FileBackend(backend_dir)
        data_manager = DataManager(file_backend)

        data_manager.load(test_image)
        md5_hex = _md5_hexdigest_from_file(test_image)
        manifest_path = os.path.join(backend_dir, md5_hex, "manifest.json")

        microscopy_collection = MicroscopyCollection()
        microscopy_collection.parse_manifest(manifest_path)

        plasmodesmata_analysis(microscopy_collection, 0, 15000, 2, 50)

        expected_data_dir = "/scripts/test_data/expected"
        for fname in os.listdir(expected_data_dir):
            expected_fpath = os.path.join(expected_data_dir, fname)
            result_fpath = os.path.join(output_dir, fname)
            expected_md5 = _md5_hexdigest_from_file(expected_fpath)
            result_md5 = _md5_hexdigest_from_file(result_fpath)
            assert expected_md5 == result_md5
    finally:
        for fname in os.listdir(output_dir):
            fpath = os.path.join(output_dir, fname)
            os.unlink(fpath)
        os.rmdir(output_dir)



if __name__ == "__main__":
    test_plasmodesmata_analysis()
