"""Module for getting cancer hotspots data.

Downloads Site: https://www.cancerhotspots.org/#/download
Data URL: https://www.cancerhotspots.org/files/hotspots_v2.xls
"""
from os import remove
import shutil
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
import boto3
from botocore.config import Config

from evidence import DATA_DIR_PATH, logger
from evidence.schemas import SourceMeta, Response, Sources


class CancerHotspots:
    """Class for Cancer Hotspots Data Access."""

    def __init__(
        self, data_url: str = "https://www.cancerhotspots.org/files/hotspots_v2.xls",
        src_dir_path: Path = DATA_DIR_PATH / "cancer_hotspots",
        normalized_data_path: Optional[Path] = None,
        ignore_normalized_data: bool = False
    ) -> None:
        """Initialize Cancer Hotspots class

        :param str data_url: URL to data file
        :param Path src_dir_path: Path to cancer hotspots data directory
        :param Optional[Path] normalized_data_path: Path to normalized cancer
            hotspots file
        :param bool ignore_normalized_data: `True` if only bare init is needed. This
            is intended for developers when using the CLI to normalize cancer hotspots
            data. Ignores path set in `normalized_data_path`.
            `False` will load normalized data from s3 and load normalized
            excel sheet data.
        """
        self.data_url = data_url
        fn = self.data_url.split("/")[-1]
        self.src_dir_path = src_dir_path
        self.src_dir_path.mkdir(exist_ok=True, parents=True)
        self.data_path = self.src_dir_path / fn
        self.og_snv_sheet_name = "SNV-hotspots"
        self.og_indel_sheet_name = "INDEL-hotspots"
        self.new_snv_sheet_name = "snv_hotspots"
        self.new_indel_sheet_name = "indel_hotspots"
        self.source_meta = SourceMeta(label=Sources.CANCER_HOTSPOTS, version="2")

        if not ignore_normalized_data:
            self.normalized_data_path = None
            if not normalized_data_path:
                self.get_normalized_data_path()
            else:
                if normalized_data_path.exists():
                    self.normalized_data_path = normalized_data_path
                else:
                    logger.error(f"The supplied path at `normalized_data_path`, "
                                 f"{normalized_data_path}, for Cancer Hotspots does "
                                 f"not exist.")

            if not self.normalized_data_path:
                raise FileNotFoundError(
                    "Unable to retrieve path for normalized Cancer Hotspots data")

            self.snv_hotspots = pd.read_excel(
                self.normalized_data_path, sheet_name=self.og_snv_sheet_name)
            self.indel_hotspots = pd.read_excel(
                self.normalized_data_path, sheet_name=self.og_indel_sheet_name)

    def get_normalized_data_path(self) -> None:
        """Download latest normalized data from public s3 bucket if it does not already
        exist in data dir and set normalized_data_path
        """
        logger.info("Retrieving normalized data from s3 bucket...")
        s3 = boto3.resource("s3", config=Config(region_name="us-east-2"))
        bucket = sorted(list(s3.Bucket("vicc-normalizers").objects.filter(
            Prefix="evidence_normalization/cancer_hotspots/normalized_hotspots_v").all()), key=lambda o: o.key)  # noqa: E501
        if len(bucket) > 0:
            obj = bucket.pop().Object()
            obj_s3_path = obj.key
            zip_fn = obj_s3_path.split("/")[-1]
            fn = zip_fn[:-4]
            normalized_data_path = self.src_dir_path / fn
            if not normalized_data_path.exists():
                zip_path = self.src_dir_path / zip_fn
                with open(zip_path, "wb") as f:
                    obj.download_fileobj(f)
                shutil.unpack_archive(zip_path, self.src_dir_path)
                remove(zip_path)
                logger.info("Successfully downloaded normalized Cancer Hotspots data")
            else:
                logger.info("Latest normalized Cancer Hotspots data already exists")
            self.normalized_data_path = normalized_data_path
        else:
            logger.warning("Could not find normalized Cancer Hotspots"
                           " data in vicc-normalizers s3 bucket")

    def mutation_hotspots(self, so_id: str, vrs_variation_id: str) -> Response:
        """Get cancer hotspot data for a variant

        :param str so_id: The structural type of the variation
        :param str vrs_variation_id: The VRS digest for the variation
        :return: Mutation hotspots data for variation
        """
        # SO:0001606 - Amino acid
        # SO:0001017 - Silent mutation
        if so_id in ["SO:0001606", "SO:0001017"]:
            data = self.query_snv_hotspots(vrs_variation_id)
        else:
            data = self.query_indel_hotspots(vrs_variation_id)

        return Response(
            data=data if data else dict(),
            source_meta_=self.source_meta
        )

    def query_snv_hotspots(self, vrs_variation_id: str) -> Optional[Dict]:
        """Return data for SNV

        :param str vrs_variation_id: VRS digest for variation
        :return: SNV data for vrs_variation_id
        """
        df = self.snv_hotspots.loc[self.snv_hotspots["vrs_identifier"] == vrs_variation_id]  # noqa: E501
        if df.empty:
            return None

        ref = df["ref"].item()
        pos = df["Amino_Acid_Position"].item()
        alt = df["Variant_Amino_Acid"].item()
        mutation, observations = alt.split(":")
        return {
            "codon": f"{ref}{pos}",
            "mutation": f"{ref}{pos}{mutation}",
            "q_value": df["qvalue"].item(),
            "observations": int(observations),
            "total_observations": int(df["Mutation_Count"].item())
        }

    def query_indel_hotspots(self, vrs_variation_id: str) -> Optional[Dict]:
        """Return data for indel

        :param str vrs_variation_id: VRS digest for variation
        :return: INDEL data for vrs_variation_id
        """
        df = self.indel_hotspots.loc[self.indel_hotspots["vrs_identifier"] == vrs_variation_id]  # noqa: E501
        if df.empty:
            return None

        pos = df["Amino_Acid_Position"].item()
        alt = df["Variant_Amino_Acid"].item()
        mutation, observations = alt.split(":")
        return {
            "codon": pos,
            "mutation": mutation,
            "q_value": df["qvalue"].item(),
            "observations": int(observations),
            "total_observations": int(df["Mutation_Count"].item())
        }
