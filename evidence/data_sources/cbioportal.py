"""Module for accessing python client for cBioPortal."""
from typing import Optional
from pathlib import Path
import shutil
from os import remove

import boto3
import xlrd

from evidence import DATA_DIR_PATH, logger
from evidence.schemas import SourceMeta, Response, Sources


class CBioPortal:
    """cBioPortal class."""

    def __init__(
        self, data_url: str = "https://cbioportal-datahub.s3.amazonaws.com/msk_impact_2017.tar.gz",  # noqa: E501
        src_dir_path: Path = DATA_DIR_PATH / "cbioportal",
        transformed_data_path: Optional[Path] = None,
        ignore_transformed_data: bool = False
    ) -> None:
        """Initialize cbioportal class

        :param str data_url: URL to data file
        :param Path src_dir_path: Path to cbioportal data directory
        :param Optional[Path] transformed_data_path: Path to transformed cbioportal file
        :param bool ignore_transformed_data: `True` if only bare init is needed. This
            is intended for developers when using the CLI to transform cbioportal data.
            Ignores path set in `transformed_data_path`. `False` will load transformed
            data from s3 and load transformed excel sheet data.
        """
        self.data_url = data_url
        self.src_dir_path = src_dir_path
        self.src_dir_path.mkdir(exist_ok=True, parents=True)
        self.source_meta = SourceMeta(
            label=Sources.CBIOPORTAL,
            version="msk_impact_2017"
        )

        if not ignore_transformed_data:
            self.transformed_data_path = None
            if not transformed_data_path:
                self.get_transformed_data_path()
            else:
                if transformed_data_path.exists():
                    self.transformed_data_path = transformed_data_path
                else:
                    logger.error(f"The supplied path at `transformed_data_path`, "
                                 f"{transformed_data_path}, for cBioPortal does not "
                                 f"exist.")

            if not self.transformed_data_path:
                raise FileNotFoundError(
                    "Unable to retrieve path for transformed cBioPortal data")

            wb = xlrd.open_workbook(self.transformed_data_path)
            self.mutations = wb.sheet_by_name("mutations")
            self.mutations_headers = self.mutations.row_values(0)
            self.case_lists = wb.sheet_by_name("case_lists")
            self.case_lists_headers = self.case_lists.row_values(0)

    def get_transformed_data_path(self) -> None:
        """Download MSK Impact 2017 data from public s3 bucket if it does not already
        exist in data dir and set transformed_data_path
        """
        logger.info("Retrieving normalized data from s3 bucket...")
        s3 = boto3.client("s3")
        zip_fn = "msk_impact_2017.xls.zip"
        zip_path = self.src_dir_path / zip_fn
        with open(zip_path, "wb") as f:
            s3.download_fileobj("vicc-normalizers",
                                f"evidence_normalization/cbioportal/{zip_fn}", f)
        shutil.unpack_archive(zip_path, self.src_dir_path)
        remove(zip_path)
        logger.info("Successfully downloaded transformed cBioPortal data")
        self.transformed_data_path = self.src_dir_path / "msk_impact_2017.xls"

    def cancer_types_summary(self, hgnc_symbol: str) -> Response:
        """Get cancer types with gene mutations data

        :param str hgnc_symbol: HGNC symbol
        :return: Cancer types summary for gene
        """
        hgnc_symbol = hgnc_symbol.upper()
        mutation_sample_ids = set()
        for row_idx in range(1, self.mutations.nrows):
            row = self.mutations.row_values(row_idx)
            if row[self.mutations_headers.index("Hugo_Symbol")] == hgnc_symbol:
                sample_id = row[self.mutations_headers.index("Tumor_Sample_Barcode")]
                mutation_sample_ids.add(sample_id)

        tumor_type_totals = dict()
        for row_idx in range(1, self.case_lists.nrows):
            row = self.case_lists.row_values(row_idx)
            case_list_name = row[self.case_lists_headers.index("case_list_name")]
            if ":" in case_list_name:
                tumor_type = case_list_name.split(": ")[-1]
                sample_ids = \
                    row[self.case_lists_headers.index("case_list_ids")].split("\t")
                tumor_type_totals[tumor_type] = {
                    "count": 0,
                    "total": len(sample_ids)
                }
                for sample_id in sample_ids:
                    if sample_id in mutation_sample_ids:
                        tumor_type_totals[tumor_type]["count"] += 1
                tumor_type_totals[tumor_type]["percent_altered"] = (tumor_type_totals[tumor_type]["count"] / tumor_type_totals[tumor_type]["total"]) * 100  # noqa: E501

        return Response(
            data=tumor_type_totals,
            source_meta_=self.source_meta
        )
