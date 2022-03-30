"""Module for getting cancer hotspots data.

Downloads Site: https://www.cancerhotspots.org/#/download
Data URL: https://www.cancerhotspots.org/files/hotspots_v2.xls
"""
from os import remove
import shutil
from pathlib import Path
from typing import Dict, Optional, List, Tuple
import csv

import boto3
from botocore.config import Config

from evidence import DATA_DIR_PATH, logger
from evidence.data_sources.base import DataSource
from evidence.schemas import SourceMeta, Response, Sources


class CancerHotspots(DataSource):
    """Class for Cancer Hotspots Data Access."""

    def __init__(
        self, data_url: str = "https://www.cancerhotspots.org/files/hotspots_v2.xls",
        src_dir_path: Path = DATA_DIR_PATH / "cancer_hotspots",
        snv_normalized_data_path: Optional[Path] = None,
        indel_normalized_data_path: Optional[Path] = None,
        ignore_normalized_data: bool = False
    ) -> None:
        """Initialize Cancer Hotspots class

        :param str data_url: URL to data file
        :param Path src_dir_path: Path to cancer hotspots data directory
        :param Optional[Path] snv_normalized_data_path: Path to normalized cancer
            hotspots SNV file
        :param Optional[Path] indel_normalized_data_path: Path to normalized cancer
            hotspots INDEL file
        :param bool ignore_normalized_data: `True` if only bare init is needed. This
            is intended for developers when using the CLI to normalize cancer hotspots
            data. Ignores path set in `normalized_data_path`.
            `False` will load normalized data from s3 and load normalized
            excel sheet data.
        """
        self.data_url = data_url
        self.src_dir_path = src_dir_path
        self.src_dir_path.mkdir(exist_ok=True, parents=True)

        self.source_meta = SourceMeta(label=Sources.CANCER_HOTSPOTS, version="2")
        self.snv_normalized_data_path = None
        self.indel_normalized_data_path = None

        if not ignore_normalized_data:
            if snv_normalized_data_path:
                if snv_normalized_data_path.exists():
                    self.snv_normalized_data_path = snv_normalized_data_path
                else:
                    logger.error(f"The supplied path at `snv_normalized_data_path`, "
                                 f"{snv_normalized_data_path}, for Cancer Hotspots "
                                 f"does not exist.")
            else:
                self.get_normalized_data_path(is_snv=True)

            if not self.snv_normalized_data_path:
                raise FileNotFoundError("Unable to retrieve path for normalized "
                                        "Cancer Hotspots SNV data")

            if indel_normalized_data_path:
                if indel_normalized_data_path.exists():
                    self.indel_normalized_data_path = indel_normalized_data_path
                else:
                    logger.error(f"The supplied path at `indel_normalized_data_path`, "
                                 f"{indel_normalized_data_path}, for Cancer Hotspots "
                                 f"does not exist.")
            else:
                self.get_normalized_data_path(is_snv=False)

            if not self.indel_normalized_data_path:
                raise FileNotFoundError("Unable to retrieve path for normalized "
                                        "Cancer Hotspots INDEL data")

    def get_normalized_data_path(self, is_snv: bool = True) -> None:
        """Download Cancer Hotspots SNV and INDEL data from public s3 bucket if it
        does not already exist in data directory and set the corresponding data path

        :param bool is_snv: `True` if getting SNV data. `False` if getting INDEL data
        """
        data_type = "snv" if is_snv else "indel"
        logger.info(f"Retrieving normalized {data_type} data from s3 bucket...")
        s3 = boto3.resource("s3", config=Config(region_name="us-east-2"))
        prefix = \
            f"evidence_normalization/cancer_hotspots/normalized_{data_type}_hotspots_v"
        bucket = sorted(list(s3.Bucket("vicc-normalizers").objects.filter(Prefix=prefix).all()), key=lambda o: o.key)  # noqa: E501
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
                logger.info(f"Successfully downloaded normalized Cancer Hotspots "
                            f"{data_type} data")
            else:
                logger.info(f"Latest normalized Cancer Hotspots {data_type} data "
                            f"already exists")
            if is_snv:
                self.snv_normalized_data_path = normalized_data_path
            else:
                self.indel_normalized_data_path = normalized_data_path
        else:
            logger.warning(f"Could not find normalized Cancer Hotspots {data_type}"
                           f" data in vicc-normalizers s3 bucket")

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
        return self.format_response(
            Response(data=data if data else dict(), source_meta_=self.source_meta)
        )

    @staticmethod
    def get_row(normalized_data_path: Path,
                vrs_identifier: str) -> Tuple[Optional[List], Optional[List]]:
        """Get row from xls sheet if vrs_identifier matches value in last column

        :param Path normalized_data_path: Path to normalized data file
        :param str vrs_identifier: The vrs_identifier to match on
        :return: Row represented as a list if vrs_identifier match was found and
            headers if match was found
        """
        matched_row = None
        headers = None
        with open(normalized_data_path) as f:
            data = csv.reader(f)
            headers = next(data)
            for row in data:
                if row[headers.index("vrs_identifier")] == vrs_identifier:
                    matched_row = row
                    break
        return matched_row, headers

    def query_snv_hotspots(self, vrs_variation_id: str) -> Optional[Dict]:
        """Return data for SNV

        :param str vrs_variation_id: VRS digest for variation
        :return: SNV data for vrs_variation_id
        """
        row, headers = self.get_row(self.snv_normalized_data_path, vrs_variation_id)
        if not row:
            return None

        ref = row[headers.index("ref")]
        pos = row[headers.index("Amino_Acid_Position")]
        alt = row[headers.index("Variant_Amino_Acid")]
        mutation, observations = alt.split(":")
        return {
            "codon": f"{ref}{pos}",
            "mutation": f"{ref}{pos}{mutation}",
            "q_value": float(row[headers.index("qvalue")]),
            "observations": int(observations),
            "total_observations": int(row[headers.index("Mutation_Count")])
        }

    def query_indel_hotspots(self, vrs_variation_id: str) -> Optional[Dict]:
        """Return data for indel

        :param str vrs_variation_id: VRS digest for variation
        :return: INDEL data for vrs_variation_id
        """
        row, headers = self.get_row(self.indel_normalized_data_path, vrs_variation_id)
        if not row:
            return None

        pos = row[headers.index("Amino_Acid_Position")]
        alt = row[headers.index("Variant_Amino_Acid")]
        mutation, observations = alt.split(":")
        return {
            "codon": pos,
            "mutation": mutation,
            "q_value": float(row[headers.index("qvalue")]),
            "observations": int(observations),
            "total_observations": int(row[headers.index("Mutation_Count")])
        }
