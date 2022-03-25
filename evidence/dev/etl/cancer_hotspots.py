"""Module for ETL cancer hotspots data"""
import logging
from timeit import default_timer as timer
from datetime import datetime

import requests
import pandas as pd
from variation.query import QueryHandler

from evidence.data_sources import CancerHotspots


class CancerHotspotsETLException(Exception):
    """Exceptions for Cancer Hotspots ETL"""

    pass


logger = logging.getLogger("evidence.etl.cancer_hotspots")
logger.setLevel(logging.DEBUG)


class CancerHotspotsETL():
    """Class for Cancer Hotspots ETL methods."""

    def __init__(self) -> None:
        """Initialize the CancerHotspotsETL class"""
        self.cancer_hotspots = CancerHotspots(ignore_normalized_data=True)

    def download_data(self) -> None:
        """Download Cancer Hotspots data."""
        if not self.cancer_hotspots.data_path.exists():
            r = requests.get(self.cancer_hotspots.data_url)
            if r.status_code == 200:
                with open(self.cancer_hotspots.data_path, "wb") as f:
                    f.write(r.content)

    def add_vrs_identifier_to_data(self) -> None:
        """Normalize variations in cancer hotspots SNV sheet and adds `vrs_identifier`
        column to dataframe. Run manually each time variation-normalizer
        or Cancer Hotspots releases a new version.

        :param CancerHotspots self.cancer_hotspots: Cancer Hotspots data source
        """
        self.download_data()
        if not self.cancer_hotspots.data_path.exists():
            raise CancerHotspotsETLException("Downloading Cancer Hotspots data"
                                             " was unsuccessful")

        snv_hotspots = pd.read_excel(self.cancer_hotspots.data_path,
                                     sheet_name=self.cancer_hotspots.og_snv_sheet_name)
        indel_hotspots = pd.read_excel(self.cancer_hotspots.data_path,
                                       sheet_name=self.cancer_hotspots.og_indel_sheet_name)  # noqa: E501
        variation_normalizer = QueryHandler()

        logger.info("Normalizing Cancer Hotspots data...")
        start = timer()
        self.get_normalized_data(
            snv_hotspots, variation_normalizer, self.cancer_hotspots.new_snv_sheet_name)
        self.get_normalized_data(
            indel_hotspots, variation_normalizer, self.cancer_hotspots.new_indel_sheet_name)  # noqa: E501
        end = timer()
        logger.info(f"Normalized Cancer Hotspots data in {(end-start):.5f} s")

        today = datetime.strftime(datetime.today(), "%Y%m%d")
        self.cancer_hotspots.normalized_data_path = \
            self.cancer_hotspots.src_dir_path / f"normalized_hotspots_v{self.cancer_hotspots.source_meta.version}_{today}.xls"  # noqa: E501
        with pd.ExcelWriter(self.cancer_hotspots.normalized_data_path) as writer:
            snv_hotspots.to_excel(
                writer, sheet_name=self.cancer_hotspots.og_snv_sheet_name, index=False)
            indel_hotspots.to_excel(
                writer, sheet_name=self.cancer_hotspots.og_indel_sheet_name, index=False)  # noqa: E501
        logger.info(f"Successfully normalized Cancer Hotspots data. "
                    f"Normalized data can be found at: {self.cancer_hotspots.normalized_data_path}")  # noqa: E501

    def get_normalized_data(self, df: pd.DataFrame,
                            variation_normalizer: QueryHandler, df_name: str) -> None:
        """Normalize variant and add vrs_identifier column to df

        :param CancerHotspots self.cancer_hotspots: Cancer Hotspots data source
        :param pd.DataFrame df: Dataframe to normalize
        :param QueryHandler variation_normalizer: Variation Normalizer handler
        :param str df_name: Name of df.
            Must be either `snv_hotspots` or `indel_hotspots`
        """
        df["vrs_identifier"] = None
        for i, row in df.iterrows():
            if df_name == self.cancer_hotspots.new_snv_sheet_name:
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
