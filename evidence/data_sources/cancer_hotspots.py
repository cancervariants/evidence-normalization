"""Module for getting cancer hotspots data.

Downloads Site: https://www.cancerhotspots.org/#/download
Data URL: https://www.cancerhotspots.org/files/hotspots_v2.xls
"""
from typing import Dict, Optional

import pandas as pd
import requests
from variation.query import QueryHandler

from evidence import DATA_ROOT, logger
from evidence.schemas import SourceMeta, Response, Sources


class CancerHotspots:
    """Class for Cancer Hotspots Data Access."""

    def __init__(self,
                 data_url: str = "https://www.cancerhotspots.org/files/hotspots_v2.xls",
                 fn: str = "hotspots_v2.xls") -> None:
        """Initialize Cancer Hotspots class"""
        self.data_url = data_url
        self.data_path = DATA_ROOT / fn
        self.normalized_data_path = DATA_ROOT / f"normalized_{fn}"
        if not self.normalized_data_path.exists():
            # If normalized hotspots file doesn't exist, create it. Will take some time
            self._add_vrs_identifier_to_data()
        self.snv_hotspots = pd.read_excel(self.normalized_data_path,
                                          sheet_name="SNV-hotspots")
        self.indel_hotspots = pd.read_excel(self.normalized_data_path,
                                            sheet_name="INDEL-hotspots")
        self.source_meta = SourceMeta(
            label=Sources.CANCER_HOTSPOTS.value,
            version="2"
        )

    def _add_vrs_identifier_to_data(self) -> None:
        """Normalize variations in cancer hotspots SNV sheet and adds `vrs_identifier`
        column to dataframe. This should be run each time variation-normalizer
        or Cancer Hotspots releases a new version.
        """
        snv_hotspots = pd.read_excel(self.data_path, sheet_name="SNV-hotspots")
        indel_hotspots = pd.read_excel(self.data_path, sheet_name="INDEL-hotspots")
        variation_normalizer = QueryHandler()

        self._get_normalized_data(snv_hotspots, variation_normalizer, "snv_hotspots")
        self._get_normalized_data(indel_hotspots, variation_normalizer,
                                  "indel_hotspots")

        with pd.ExcelWriter(self.normalized_data_path) as writer:
            snv_hotspots.to_excel(writer, sheet_name="SNV-hotspots", index=False)
            indel_hotspots.to_excel(writer, sheet_name="INDEL-hotspots", index=False)

    @staticmethod
    def _get_normalized_data(df: pd.DataFrame,
                             variation_normalizer: QueryHandler, df_name: str) -> None:
        """Normalize variant and add vrs_identifier column to df

        :param pd.DataFrame df: Dataframe to normalize
        :param QueryHandler variation_normalizer: Variation Normalizer handler
        :param str df_name: Name of df.
            Must be either `snv_hotspots` or `indel_hotspots`
        """
        df["vrs_identifier"] = None
        for i, row in df.iterrows():
            if df_name == "snv_hotspots":
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

    def download_data(self) -> None:
        """Download Cancer Hotspots data."""
        if not self.data_path.exists():
            r = requests.get(self.data_url)
            if r.status_code == 200:
                with open(self.data_path, "wb") as f:
                    f.write(r.content)

    def mutation_hotspots(self, so_id: str, vrs_variation_id: str) -> Optional[Dict]:
        """Get cancer hotspot data for a variant

        :param str so_id: The structural type of the variation
        :param str vrs_variation_id: The VRS digest for the variation
        """
        # SO:0001606 - Amino acid
        # SO:0001017 - Silent mutation
        if so_id in ["SO:0001606", "SO:0001017"]:
            data = self.query_snv_hotspots(vrs_variation_id)
        else:
            data = self.query_indel_hotspots(vrs_variation_id)

        return Response(
            data=data if data else None,
            source_meta_=self.source_meta
        )

    def query_snv_hotspots(self, vrs_variation_id: str) -> Optional[Dict]:
        """Return data for SNV"""
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
        """Return data for indel"""
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
