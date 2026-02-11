"""Module for getting cancer hotspots data.

Downloads Site: https://www.cancerhotspots.org/#/download
Data URL: https://www.cancerhotspots.org/files/hotspots_v2.xls
"""

import json
import logging
import shutil
from pathlib import Path

import boto3
from botocore.config import Config

from evidence import DATA_DIR_PATH
from evidence.data_sources.base import DownloadableDataSource
from evidence.schemas import Response, SourceDataType, SourceMeta, Sources

_logger = logging.getLogger(__name__)


class CancerHotspots(DownloadableDataSource):
    """Class for Cancer Hotspots Data Access."""

    def __init__(
        self,
        data_url: str = "https://www.cancerhotspots.org/files/hotspots_v2.xls",
        src_dir_path: Path = DATA_DIR_PATH / "cancer_hotspots",
        transformed_data_path: Path | None = None,
        ignore_transformed_data: bool = False,
    ) -> None:
        """Initialize Cancer Hotspots class

        :param data_url: URL to data file
        :param src_dir_path: Path to cancer hotspots data directory
        :param transformed_data_path: Path to transformed cancer hotspots file
        :param ignore_transformed_data: `True` if only bare init is needed. This is
            intended for developers when using the CLI to transform cancer hotspots
            data. Ignores path set in `_transformed_data_path`. `False` will load
            transformed data from s3
        """
        super().__init__(data_url, src_dir_path, ignore_transformed_data)

        self.source_meta = SourceMeta(label=Sources.CANCER_HOTSPOTS, version="2")
        transformed_data_path = self.get_transformed_data_path(
            transformed_data_path, SourceDataType.CANCER_HOTSPOTS
        )
        if transformed_data_path:
            with transformed_data_path.open("r") as f:
                self.transformed_data = json.load(f)
        else:
            self.transformed_data = {}

    def download_s3_data(
        self, src_data_type: SourceDataType = SourceDataType.CANCER_HOTSPOTS
    ) -> None:
        """Download Cancer Hotspots data from public s3 bucket if it
        does not already exist in data directory and set the corresponding data path

        :param src_data_type: The data type contained in the transformed data file
        :return: Path to transformed data file
        """
        data_path = None
        _logger.info("Retrieving transformed cancer hotspots data from s3 bucket...")
        s3 = boto3.resource("s3", config=Config(region_name="us-east-2"))
        prefix = f"evidence_normalization/cancer_hotspots/{src_data_type.value}_"
        bucket = sorted(
            s3.Bucket("vicc-normalizers").objects.filter(Prefix=prefix).all(),
            key=lambda o: o.key,
        )
        if len(bucket) > 0:
            obj = bucket.pop().Object()
            obj_s3_path = obj.key
            zip_fn = obj_s3_path.split("/")[-1]
            fn = zip_fn[:-4]
            transformed_data_path = self.src_dir_path / fn
            if not transformed_data_path.exists():
                zip_path = self.src_dir_path / zip_fn
                with zip_path.open("wb") as f:
                    obj.download_fileobj(f)
                shutil.unpack_archive(zip_path, self.src_dir_path)
                Path.unlink(zip_path)
                _logger.info("Successfully downloaded transformed Cancer Hotspots data")
            else:
                _logger.info("Latest transformed Cancer Hotspots data already exists")

            data_path = transformed_data_path
        else:
            _logger.warning(
                "Could not find transformed Cancer Hotspots data in vicc-normalizers "
                "s3 bucket"
            )
        return data_path

    def mutation_hotspots(self, vrs_variation_id: str) -> Response:
        """Get cancer hotspot data for a variant

        :param vrs_variation_id: The VRS digest for the variation
        :return: Mutation hotspots data for variation
        """
        return self.format_response(
            Response(
                data=self.transformed_data.get(vrs_variation_id, {}),
                source_meta_=self.source_meta,
            )
        )
