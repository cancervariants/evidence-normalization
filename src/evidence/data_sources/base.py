"""Module for the base data source class"""

import hashlib
import json
import logging
from pathlib import Path

from evidence.schemas import Response, SourceDataType

_logger = logging.getLogger(__name__)


class DataSource:
    """A base class for data sources"""

    @staticmethod
    def format_response(resp: Response) -> Response:
        """Add `id` to resp object if data exists

        :param Response resp: Response object
        :return: Response object with `id` field added if data exists
        """
        if resp.data:
            blob = json.dumps(
                resp.model_dump(), sort_keys=True, separators=(",", ":"), indent=None
            ).encode("utf-8")
            digest = hashlib.md5(blob)  # noqa: S324
            resp.id = f"normalize.evidence:{digest.hexdigest()}"
        return resp


class DownloadableDataSource(DataSource):
    """A base class for sources that use downloadable data"""

    def __init__(
        self, data_url: str, src_dir_path: Path, ignore_transformed_data: bool
    ) -> None:
        """Initialize DownloadableDataSource class

        :param str data_url: URL to data file
        :param Path src_dir_path: Path to source data directory
        :param bool ignored_transformed_data: `True` if only bare init is needed. This
            is intended for developers when using the CLI to transform source data.
            `False` will load the transformed data from s3
        """
        self.data_url = data_url
        self.src_dir_path = src_dir_path
        self.src_dir_path.mkdir(exist_ok=True, parents=True)
        self.ignore_transformed_data = ignore_transformed_data

    def download_s3_data(self, src_data_type: SourceDataType) -> Path:
        """Download data from public s3 bucket if it does not already exist in data
        directory and set the corresponding data path

        :param SourceDataType src_data_type: The data type contained in the
            transformed data file
        """
        raise NotImplementedError

    def get_transformed_data_path(
        self, transformed_data_path: Path, src_data_type: SourceDataType
    ) -> Path | None:
        """Get transformed data path for source

        :param Path transformed_data_path: The path to the transformed data file
        :param SourceDataType src_data_type: The data type contained in the
            transformed data file
        :return: Path to transformed data file
        """
        data_path = None
        if not self.ignore_transformed_data:
            if transformed_data_path:
                if transformed_data_path.exists():
                    data_path = transformed_data_path
                else:
                    _logger.error(
                        "The supplied path at %s does not exist.", transformed_data_path
                    )
            else:
                data_path = self.download_s3_data(src_data_type)
        return data_path
