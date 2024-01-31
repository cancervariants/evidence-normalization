"""Module for ETL cancer hotspots data"""
import json
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


class CancerHotspotsETLError(Exception):
    """Exceptions for Cancer Hotspots ETL"""


logger = logging.getLogger(__name__)


class CancerHotspotsETL(CancerHotspots):
    """Class for Cancer Hotspots ETL methods."""

    def __init__(
        self, data_url: str = "https://www.cancerhotspots.org/files/hotspots_v2.xls",
        src_dir_path: Path = DATA_DIR_PATH / "cancer_hotspots",
        transformed_data_path: Optional[Path] = DATA_DIR_PATH / "cancer_hotspots",
        ignore_transformed_data: bool = True
    ) -> None:
        """Initialize the CancerHotspotsETL class

        :param data_url: URL to data file
        :param src_dir_path: Path to cancer hotspots data directory
        :param transformed_data_path: Path to transformed cancer hotspots file
        :param ignore_transformed_data: `True` if only bare init is needed. This is
            intended for developers when using the CLI to transform cancer hotspots
            data. Ignores path set in `_transformed_data_path`. `False` will load
            transformed data from s3
        """
        super().__init__(
            data_url, src_dir_path, transformed_data_path, ignore_transformed_data
        )
        fn = self.data_url.split("/")[-1]
        self.data_path = self.src_dir_path / fn
        self.transformed_data = {}  # vrs_id: hotspot data

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
        """Normalize variations in cancer hotspots and updates `transformed_data`

        Run manually each time variation-normalizer or Cancer Hotspots releases a new
        version.
        """
        self.download_data()
        if not self.data_path.exists():
            raise CancerHotspotsETLError(
                "Downloading Cancer Hotspots data was unsuccessful"
            )

        snv_hotspots = pd.read_excel(self.data_path, sheet_name="SNV-hotspots")
        indel_hotspots = pd.read_excel(self.data_path, sheet_name="INDEL-hotspots")
        variation_normalizer = QueryHandler()

        logger.info("Normalizing Cancer Hotspots data...")
        start = timer()
        await self.get_transformed_data(
            snv_hotspots, variation_normalizer, is_snv=True)
        await self.get_transformed_data(
            indel_hotspots, variation_normalizer, is_snv=False)
        end = timer()

        logger.info("Transformed Cancer Hotspots data in %.*f s", 2, end - start)

        today = datetime.strftime(datetime.today(), "%Y%m%d")
        transformed_data_path = self.src_dir_path / f"cancer_hotspots_{today}.json"
        with transformed_data_path.open("w") as f:
            json.dump(self.transformed_data, f)

        logger.info("Successfully transformed Cancer Hotspots data.")

    async def get_transformed_data(
        self, df: pd.DataFrame, variation_normalizer: QueryHandler, is_snv: bool
    ) -> None:
        """Normalize variant and updates `transformed_data`

        :param df: Dataframe to transform
        :param variation_normalizer: Variation Normalizer handler
        :param is_snv: `True` if SNV data, else INDEL
        """
        for _, row in df.iterrows():
            hugo_symbol = row["Hugo_Symbol"]
            alt = row["Variant_Amino_Acid"]
            pos = row["Amino_Acid_Position"]

            if is_snv:
                ref = row["ref"]
                variation = f"{hugo_symbol} {ref}{pos}{alt.split(':')[0]}"
            else:
                ref = None
                variation = f"{hugo_symbol} {alt.split(':')[0]}"

            try:
                variation_norm_resp = await variation_normalizer.normalize_handler.normalize(variation)  # noqa: E501
            except Exception as e:
                logger.error("variation-normalizer unable to normalize %s: %s", variation, str(e))  # noqa: E501
            else:
                if variation_norm_resp and variation_norm_resp.variation:
                    vrs_id = variation_norm_resp.variation.id
                    if vrs_id in self.transformed_data:
                        logger.debug(
                            "duplicate vrs_id (%s) for variation (%s)",
                            vrs_id,
                            variation
                        )

                    mutation, observations = alt.split(":")

                    if is_snv:
                        codon = f"{ref}{pos}"
                        mutation = f"{codon}{mutation}"
                    else:
                        codon = pos

                    self.transformed_data[vrs_id] = {
                        "variation": variation,
                        "codon": codon,
                        "mutation": mutation,
                        "q_value": float(row["qvalue"]),
                        "observations": int(observations),
                        "total_observations": int(row["Mutation_Count"])
                    }
                else:
                    logger.warning("variation-normalizer unable to normalize: %s", variation)  # noqa: E501
