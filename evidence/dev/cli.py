"""Dev CLI"""
from timeit import default_timer as timer
from datetime import datetime

import requests
import click
from variation.query import QueryHandler
import pandas as pd

from evidence import logger
from evidence.data_sources import CancerHotspots


def log_and_echo_msg(msg: str, log_level: str = "info") -> None:
    """Log and echo message

    :param str msg: Message
    :param str log_level: Logging level. Must be either info, warning, or error
    """
    if log_level == "info":
        logger.info(msg)
    elif log_level == "warning":
        logger.warning(msg)
    else:
        logger.error(msg)
    click.echo(msg)


@click.command()
@click.option(
    "--normalize_cancer_hotspots",
    is_flag=True,
    default=False,
    help="Normalize Cancer Hotspots data"
)
def cli(normalize_cancer_hotspots: bool) -> None:
    """Execute CLI methods

    :param bool normalize_cancer_hotspots: Determines whether or not to normalize
        Cancer Hotspots data
    """
    if normalize_cancer_hotspots:
        normalize_cancer_hotspots_data()


def normalize_cancer_hotspots_data() -> None:
    """Normalize Cancer Hotspots data"""

    def download_data(ch: CancerHotspots) -> None:
        """Download Cancer Hotspots data.

        :param CancerHotspots ch: Cancer Hotspots data source
        """
        if not ch.data_path.exists():
            log_and_echo_msg("Downloading Cancer Hotspots data...")
            r = requests.get(ch.data_url)
            if r.status_code == 200:
                with open(ch.data_path, "wb") as f:
                    f.write(r.content)

    def add_vrs_identifier_to_data(ch: CancerHotspots) -> None:
        """Normalize variations in cancer hotspots SNV sheet and adds `vrs_identifier`
        column to dataframe. Run manually each time variation-normalizer
        or Cancer Hotspots releases a new version.

        :param CancerHotspots ch: Cancer Hotspots data source
        """
        download_data(ch)
        if not ch.data_path.exists():
            log_and_echo_msg("Downloading Cancer Hotspots data was unsuccessful",
                             "warning")
            return

        snv_hotspots = pd.read_excel(ch.data_path,
                                     sheet_name=ch.og_snv_sheet_name)
        indel_hotspots = pd.read_excel(
            ch.data_path, sheet_name=ch.og_indel_sheet_name)
        variation_normalizer = QueryHandler()

        log_and_echo_msg("Normalizing Cancer Hotspots data...")
        start = timer()
        get_normalized_data(
            ch, snv_hotspots, variation_normalizer, ch.new_snv_sheet_name)
        get_normalized_data(
            ch, indel_hotspots, variation_normalizer, ch.new_indel_sheet_name)
        end = timer()
        log_and_echo_msg(f"Normalized Cancer Hotspots data in {(end-start):.5f} s")

        today = datetime.strftime(datetime.today(), "%Y%m%d")
        ch.normalized_data_path = \
            ch.src_dir_path / f"normalized_hotspots_v{ch.source_meta.version}_{today}.xlsx"  # noqa: E501
        with pd.ExcelWriter(ch.normalized_data_path) as writer:
            snv_hotspots.to_excel(
                writer, sheet_name=ch.og_snv_sheet_name, index=False)
            indel_hotspots.to_excel(
                writer, sheet_name=ch.og_indel_sheet_name, index=False)
        log_and_echo_msg(f"Successfully normalized Cancer Hotspots data. "
                         f"Normalized data can be found at: {ch.normalized_data_path}")

    def get_normalized_data(ch: CancerHotspots, df: pd.DataFrame,
                            variation_normalizer: QueryHandler, df_name: str) -> None:
        """Normalize variant and add vrs_identifier column to df

        :param CancerHotspots ch: Cancer Hotspots data source
        :param pd.DataFrame df: Dataframe to normalize
        :param QueryHandler variation_normalizer: Variation Normalizer handler
        :param str df_name: Name of df.
            Must be either `snv_hotspots` or `indel_hotspots`
        """
        df["vrs_identifier"] = None
        for i, row in df.iterrows():
            if df_name == ch.new_snv_sheet_name:
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

    cancer_hotspots = CancerHotspots(ignore_normalized_data=True)
    add_vrs_identifier_to_data(cancer_hotspots)


if __name__ == "__main__":
    cli()
