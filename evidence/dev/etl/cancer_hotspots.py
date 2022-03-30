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
        normalized_data_path: Optional[Path] = None,
        ignore_normalized_data: bool = True
    ) -> None:
        """Initialize the CancerHotspotsETL class"""
        super().__init__(
            data_url, src_dir_path, normalized_data_path, ignore_normalized_data
        )

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

    def add_vrs_identifier_to_data(self) -> None:
        """Normalize variations in cancer hotspots SNV sheet and adds `vrs_identifier`
        column to dataframe. Run manually each time variation-normalizer
        or Cancer Hotspots releases a new version.

        :param CancerHotspots self: Cancer Hotspots data source
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
        self.get_normalized_data(
            snv_hotspots, variation_normalizer, self.new_snv_sheet_name)
        self.get_normalized_data(
            indel_hotspots, variation_normalizer, self.new_indel_sheet_name)
        end = timer()
        logger.info(f"Normalized Cancer Hotspots data in {(end-start):.5f} s")

        today = datetime.strftime(datetime.today(), "%Y%m%d")
        self.normalized_data_path = \
            self.src_dir_path / f"normalized_hotspots_v{self.source_meta.version}_{today}.xls"  # noqa: E501
        with pd.ExcelWriter(self.normalized_data_path) as writer:
            snv_hotspots.to_excel(
                writer, sheet_name=self.og_snv_sheet_name, index=False)
            indel_hotspots.to_excel(
                writer, sheet_name=self.og_indel_sheet_name, index=False)
        logger.info(f"Successfully normalized Cancer Hotspots data. "
                    f"Normalized data can be found at: {self.normalized_data_path}")

    def get_normalized_data(self, df: pd.DataFrame,
                            variation_normalizer: QueryHandler, df_name: str) -> None:
        """Normalize variant and add vrs_identifier column to df

        :param CancerHotspots self: Cancer Hotspots data source
        :param pd.DataFrame df: Dataframe to normalize
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
                norm_vd = variation_normalizer.normalize(variation)
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
