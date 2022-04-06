"""Module for ETL cancer hotspots data"""
import logging
from timeit import default_timer as timer
from datetime import datetime
from typing import Optional
from pathlib import Path

import requests
import pandas as pd
from variation.query import QueryHandler

from evidence import DATA_DIR_PATH
from evidence.data_sources import CancerHotspots


class CancerHotspotsETLException(Exception):
    """Exceptions for Cancer Hotspots ETL"""

    pass


logger = logging.getLogger("evidence.etl.cancer_hotspots")
logger.setLevel(logging.DEBUG)


class CancerHotspotsETL(CancerHotspots):
    """Class for Cancer Hotspots ETL methods."""

    def __init__(
        self, data_url: str = "https://www.cancerhotspots.org/files/hotspots_v2.xls",
        src_dir_path: Path = DATA_DIR_PATH / "cancer_hotspots",
        snv_transformed_data_path: Optional[Path] = None,
        indel_transformed_data_path: Optional[Path] = None,
        ignore_transformed_data: bool = True
    ) -> None:
        """Initialize the CancerHotspotsETL class

        :param str data_url: URL to data file
        :param Path src_dir_path: Path to cancer hotspots data directory
        :param Optional[Path] snv_transformed_data_path: Path to transformed cancer
            hotspots SNV file
        :param Optional[Path] indel_transformed_data_path: Path to transformed cancer
            hotspots INDEL file
        :param bool ignore_transformed_data: `True` if only bare init is needed. This
            is intended for developers when using the CLI to normalize cancer hotspots
            data. Ignores paths set in `snv_transformed_data_path` and
            `indel_transformed_data_path`. `False` will load transformed data from s3
        """
        super().__init__(
            data_url, src_dir_path, snv_transformed_data_path,
            indel_transformed_data_path, ignore_transformed_data
        )
        fn = self.data_url.split("/")[-1]
        self.data_path = self.src_dir_path / fn
        self.og_snv_sheet_name = "SNV-hotspots"
        self.og_indel_sheet_name = "INDEL-hotspots"
        self.new_snv_sheet_name = "snv_hotspots"
        self.new_indel_sheet_name = "indel_hotspots"

    def download_data(self) -> None:
        """Download Cancer Hotspots data."""
        if not self.data_path.exists():
            r = requests.get(self.data_url)
            if r.status_code == 200:
                with open(self.data_path, "wb") as f:
                    f.write(r.content)
            else:
                logger.error(f"Unable to download Cancer Hotspots data. "
                             f"Received status code: {r.status_code}")

    async def add_vrs_identifier_to_data(self) -> None:
        """Normalize variations in cancer hotspots SNV sheet and adds `vrs_identifier`
        column to dataframe. Run manually each time variation-normalizer
        or Cancer Hotspots releases a new version.
        """
        self.download_data()
        if not self.data_path.exists():
            raise CancerHotspotsETLException("Downloading Cancer Hotspots data"
                                             " was unsuccessful")

        snv_hotspots = pd.read_excel(self.data_path,
                                     sheet_name=self.og_snv_sheet_name)
        indel_hotspots = pd.read_excel(self.data_path,
                                       sheet_name=self.og_indel_sheet_name)
        variation_normalizer = QueryHandler()

        logger.info("Normalizing Cancer Hotspots data...")
        start = timer()
        await self.get_transformed_data(
            snv_hotspots, variation_normalizer, self.new_snv_sheet_name)
        await self.get_transformed_data(
            indel_hotspots, variation_normalizer, self.new_indel_sheet_name)
        end = timer()
        logger.info(f"transformed Cancer Hotspots data in {(end-start):.5f} s")

        today = datetime.strftime(datetime.today(), "%Y%m%d")
        snv_transformed_data_path = self.src_dir_path / f"hotspots_snv_v{self.source_meta.version}_{today}.csv"  # noqa: E501
        indel_transformed_data_path = self.src_dir_path / f"hotspots_indel_v{self.source_meta.version}_{today}.csv"  # noqa: E501
        snv_hotspots.to_csv(snv_transformed_data_path)
        indel_hotspots.to_csv(indel_transformed_data_path)
        logger.info("Successfully transformed Cancer Hotspots data.")

    async def get_transformed_data(
        self, df: pd.DataFrame, variation_normalizer: QueryHandler, df_name: str
    ) -> None:
        """Normalize variant and add vrs_identifier column to df

        :param pd.DataFrame df: Dataframe to transform
        :param QueryHandler variation_normalizer: Variation Normalizer handler
        :param str df_name: Name of df.
            Must be either `snv_hotspots` or `indel_hotspots`
        """
        df["vrs_identifier"] = None
        for i, row in df.iterrows():
            if df_name == self.new_snv_sheet_name:
                variation = f"{row['Hugo_Symbol']} {row['ref']}{row['Amino_Acid_Position']}{row['Variant_Amino_Acid'].split(':')[0]}"  # noqa: E501
            else:
                variation = f"{row['Hugo_Symbol']} {row['Variant_Amino_Acid'].split(':')[0]}"  # noqa: E501
            try:
                norm_vd = await variation_normalizer.normalize(variation)
            except Exception as e:
                logger.warning(f"variation-normalizer unable to normalize {variation}: {e}")  # noqa: E501
            else:
                if norm_vd:
                    norm_vd = norm_vd.dict()
                    if norm_vd["variation"]["type"] != "Text":
                        df.at[i, "vrs_identifier"] = norm_vd["variation"]["id"]
                    else:
                        logger.warning(f"variation-normalizer unable to normalize: {variation}")  # noqa: E501
                else:
                    logger.warning(f"variation-normalizer unable to normalize: {variation}")  # noqa: E501
