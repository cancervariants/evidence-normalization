"""Module for accessing python client for cBioPortal."""

import csv
import logging
import shutil
from pathlib import Path

import boto3

from evidence import DATA_DIR_PATH
from evidence.data_sources.base import DownloadableDataSource
from evidence.schemas import Response, SourceDataType, SourceMeta, Sources

_logger = logging.getLogger(__name__)


class CBioPortal(DownloadableDataSource):
    """cBioPortal class."""

    def __init__(
        self,
        data_url: str = "https://cbioportal-datahub.s3.amazonaws.com/msk_impact_2017.tar.gz",
        src_dir_path: Path = DATA_DIR_PATH / "cbioportal",
        transformed_mutations_data_path: Path | None = None,
        transformed_case_lists_data_path: Path | None = None,
        ignore_transformed_data: bool = False,
    ) -> None:
        """Initialize cbioportal class

        :param str data_url: URL to data file
        :param Path src_dir_path: Path to cbioportal data directory
        :param Optional[Path] transformed_mutations_data_path: Path to transformed
            cbioportal mutations file
        :param Optional[Path] transformed_case_lists_data_path: Path to transformed
            cbioportal case_lists file
        :param bool ignore_transformed_data: `True` if only bare init is needed. This
            is intended for developers when using the CLI to transform cbioportal data.
            Ignores path set in `transformed_mutations_data_path` and
            `transformed_case_lists_data_path`. `False` will load transformed
            data from s3
        """
        super().__init__(data_url, src_dir_path, ignore_transformed_data)
        self.source_meta = SourceMeta(
            label=Sources.CBIOPORTAL, version="msk_impact_2017"
        )
        self.transformed_mutations_data_path = self.get_transformed_data_path(
            transformed_mutations_data_path, SourceDataType.CBIOPORTAL_MUTATIONS
        )
        self.transformed_case_lists_data_path = self.get_transformed_data_path(
            transformed_case_lists_data_path, SourceDataType.CBIOPORTAL_CASE_LISTS
        )

    def download_s3_data(
        self, src_data_type: SourceDataType = SourceDataType.CBIOPORTAL_MUTATIONS
    ) -> Path:
        """Download MSK Impact 2017 mutation and case_lists data from public s3 bucket
        if it does not already exist in data directory and set the corresponding
        data path

        :param SourceDataType src_data_type: The data type contained in the
            transformed data file
        :return: Path to transformed data file
        """
        is_mutations = src_data_type == SourceDataType.CBIOPORTAL_MUTATIONS
        data_type = "mutations" if is_mutations else "case_lists"
        _logger.info("Retrieving transformed %s data from s3 bucket...", data_type)
        s3 = boto3.client("s3")
        zip_fn = (
            "msk_impact_2017_mutations.csv.zip"
            if is_mutations
            else "msk_impact_2017_case_lists.csv.zip"
        )
        zip_path = self.src_dir_path / zip_fn
        with zip_path.open("wb") as f:
            s3.download_fileobj(
                "vicc-normalizers", f"evidence_normalization/cbioportal/{zip_fn}", f
            )
        shutil.unpack_archive(zip_path, self.src_dir_path)
        Path.unlink(zip_path)
        _logger.info(
            "Successfully downloaded transformed cBioPortal %s data", data_type
        )
        return self.src_dir_path / zip_fn[:-4]

    def cancer_types_summary(self, hgnc_symbol: str) -> Response:
        """Get cancer types with gene mutations data

        :param str hgnc_symbol: HGNC symbol
        :return: Cancer types summary for gene
        """
        hgnc_symbol = hgnc_symbol.upper()

        mutation_sample_ids = set()
        with self.transformed_mutations_data_path.open() as f:
            data = csv.reader(f)
            headers = next(data)
            for row in data:
                if row[headers.index("Hugo_Symbol")] == hgnc_symbol:
                    sample_id = row[headers.index("Tumor_Sample_Barcode")]
                    mutation_sample_ids.add(sample_id)

        if not mutation_sample_ids:
            return self.format_response(
                Response(data={}, source_meta_=self.source_meta)
            )

        tumor_type_totals = {}
        with self.transformed_case_lists_data_path.open() as f:
            data = csv.reader(f)
            headers = next(data)
            for row in data:
                case_list_name = row[headers.index("case_list_name")]
                if ":" in case_list_name:
                    tumor_type = case_list_name.split(": ")[-1]
                    sample_ids = row[headers.index("case_list_ids")].split("\t")
                    tumor_type_totals[tumor_type] = {
                        "count": 0,
                        "total": len(sample_ids),
                    }
                    for sample_id in sample_ids:
                        if sample_id in mutation_sample_ids:
                            tumor_type_totals[tumor_type]["count"] += 1
                    tumor_type_totals[tumor_type]["percent_altered"] = (
                        tumor_type_totals[tumor_type]["count"]
                        / tumor_type_totals[tumor_type]["total"]
                    ) * 100
        return self.format_response(
            Response(data=tumor_type_totals, source_meta_=self.source_meta)
        )
