"""Module for accessing python client for cBioPortal."""
from typing import Optional
from pathlib import Path
import shutil
from os import remove
import csv

import boto3

from evidence import DATA_DIR_PATH, logger
from evidence.data_sources.base import DataSource
from evidence.schemas import SourceMeta, Response, Sources


class CBioPortal(DataSource):
    """cBioPortal class."""

    def __init__(
        self, data_url: str = "https://cbioportal-datahub.s3.amazonaws.com/msk_impact_2017.tar.gz",  # noqa: E501
        src_dir_path: Path = DATA_DIR_PATH / "cbioportal",
        transformed_mutations_data_path: Optional[Path] = None,
        transformed_case_lists_data_path: Optional[Path] = None,
        ignore_transformed_data: bool = False
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
            data from s3 and load transformed excel sheet data.
        """
        self.data_url = data_url
        self.src_dir_path = src_dir_path
        self.src_dir_path.mkdir(exist_ok=True, parents=True)
        self.source_meta = SourceMeta(
            label=Sources.CBIOPORTAL,
            version="msk_impact_2017"
        )
        self.transformed_mutations_data_path = None
        self.transformed_case_lists_data_path = None

        if not ignore_transformed_data:
            if transformed_mutations_data_path:
                if transformed_mutations_data_path.exists():
                    self.transformed_mutations_data_path = \
                        transformed_mutations_data_path
                else:
                    logger.error(f"The supplied path at `transformed_mutations_data_"
                                 f"path`, {transformed_mutations_data_path}, for "
                                 f"cBioPortal does not exist.")
            else:
                self.get_transformed_data_path(is_mutations=True)

            if not self.transformed_mutations_data_path:
                raise FileNotFoundError(
                    "Unable to retrieve path for transformed cBioPortal mutations data")

            if transformed_case_lists_data_path:
                if transformed_case_lists_data_path.exists():
                    self.transformed_case_lists_data_path = \
                        transformed_case_lists_data_path
                else:
                    logger.error(f"The supplied path at `transformed_case_lists_data_"
                                 f"path`, {transformed_case_lists_data_path}, for "
                                 f"cBioPortal does not exist.")
            else:
                self.get_transformed_data_path(is_mutations=False)

            if not self.transformed_case_lists_data_path:
                raise FileNotFoundError(
                    "Unable to retrieve path for transformed cBioPortal case_lists"
                    " data")

    def get_transformed_data_path(self, is_mutations: bool = True) -> None:
        """Download MSK Impact 2017 mutation and case_lists data from public s3 bucket
        if it does not already exist in data directory and set the corresponding
        data path

        :param bool is_mutations: `True` if getting mutations data. `False` if
            getting case_lists data
        """
        data_type = "mutations" if is_mutations else "case_lists"
        logger.info(f"Retrieving transformed {data_type} data from s3 bucket...")
        s3 = boto3.client("s3")
        zip_fn = "msk_impact_2017_mutations.csv.zip" if is_mutations else "msk_impact_2017_case_lists.csv.zip"  # noqa: E501
        zip_path = self.src_dir_path / zip_fn
        with open(zip_path, "wb") as f:
            s3.download_fileobj("vicc-normalizers",
                                f"evidence_normalization/cbioportal/{zip_fn}", f)
        shutil.unpack_archive(zip_path, self.src_dir_path)
        remove(zip_path)
        logger.info(f"Successfully downloaded transformed cBioPortal {data_type} data")
        data_path = self.src_dir_path / zip_fn[:-4]
        if is_mutations:
            self.transformed_mutations_data_path = data_path
        else:
            self.transformed_case_lists_data_path = data_path

    def cancer_types_summary(self, hgnc_symbol: str) -> Response:
        """Get cancer types with gene mutations data

        :param str hgnc_symbol: HGNC symbol
        :return: Cancer types summary for gene
        """
        hgnc_symbol = hgnc_symbol.upper()

        mutation_sample_ids = set()
        with open(self.transformed_mutations_data_path) as f:
            data = csv.reader(f)
            headers = next(data)
            for row in data:
                if row[headers.index("Hugo_Symbol")] == hgnc_symbol:
                    sample_id = row[headers.index("Tumor_Sample_Barcode")]
                    mutation_sample_ids.add(sample_id)

        if not mutation_sample_ids:
            return self.format_response(Response(data=dict(),
                                                 source_meta_=self.source_meta))

        tumor_type_totals = dict()
        with open(self.transformed_case_lists_data_path) as f:
            data = csv.reader(f)
            headers = next(data)
            for row in data:
                case_list_name = row[headers.index("case_list_name")]
                if ":" in case_list_name:
                    tumor_type = case_list_name.split(": ")[-1]
                    sample_ids = row[headers.index("case_list_ids")].split("\t")
                    tumor_type_totals[tumor_type] = {
                        "count": 0,
                        "total": len(sample_ids)
                    }
                    for sample_id in sample_ids:
                        if sample_id in mutation_sample_ids:
                            tumor_type_totals[tumor_type]["count"] += 1
                    tumor_type_totals[tumor_type]["percent_altered"] = (tumor_type_totals[tumor_type]["count"] / tumor_type_totals[tumor_type]["total"]) * 100  # noqa: E501
        return self.format_response(Response(data=tumor_type_totals,
                                             source_meta_=self.source_meta))
