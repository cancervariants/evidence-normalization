"""Module for ETL cbioportal data"""

import logging
import tarfile
import warnings
from pathlib import Path

import pandas as pd
import requests

from evidence import DATA_DIR_PATH
from evidence.data_sources import CBioPortal
from evidence.dev.etl import ETL_DATA_DIR_PATH

warnings.filterwarnings("ignore")
_logger = logging.getLogger(__name__)


class CBioPortalETLError(Exception):
    """Exceptions for CBioPortal ETL"""


class CBioPortalETL(CBioPortal):
    """Class for cBioPortal ETL methods."""

    def __init__(
        self,
        data_url: str = "https://cbioportal-datahub.s3.amazonaws.com/msk_impact_2017.tar.gz",
        src_dir_path: Path = DATA_DIR_PATH / "cbioportal",
        transformed_mutations_data_path: Path | None = None,
        transformed_case_lists_data_path: Path | None = None,
        ignore_transformed_data: bool = True,
    ) -> None:
        """Initialize cbioportal etl class

        :param str data_url: URL to data file
        :param Path src_dir_path: Path to cbioportal data directory
        :param Optional[Path] transformed_mutations_data_path: Path to transformed
            cbioportal mutations file
        :param Optional[Path] transformed_case_lists_data_path: Path to transformed
            cbioportal case_lists file
        :param bool ignore_transformed_data: `True` if only bare init is needed. This
            is intended for developers when using the CLI to transform cbioportal data.
            Ignores paths set in `transformed_mutations_data_path` and
            `transformed_case_lists_data_path`. `False` will load transformed
            data from s3
        """
        super().__init__(
            data_url,
            src_dir_path,
            transformed_mutations_data_path,
            transformed_case_lists_data_path,
            ignore_transformed_data,
        )
        self.src_dir_etl_path = ETL_DATA_DIR_PATH / "cbioportal"
        self.src_dir_etl_path.mkdir(exist_ok=True, parents=True)
        self.msk_impact_2017_dir = None

    def download_data(self) -> None:
        """Download MSK Impact 2017 data."""
        response = requests.get(self.data_url, stream=True, timeout=5)
        if response.status_code == 200:
            file = tarfile.open(fileobj=response.raw, mode="r|gz")
            file.extractall(path=self.src_dir_etl_path)  # noqa: S202
            self.msk_impact_2017_dir = self.src_dir_etl_path / "msk_impact_2017"
        else:
            _logger.error(
                "Unable to download cBioPortal data. Received status code: %i",
                response.status_code,
            )

    def transform_data(self) -> None:
        """Transform cbioportal data and write to csv files"""
        if not self.msk_impact_2017_dir:
            self.download_data()

        if not self.msk_impact_2017_dir:
            err_msg = "Downloading cBioPortal data was unsuccessful"
            raise CBioPortalETLError(err_msg)

        case_lists_df = self.create_case_lists_df()
        mutations_df = self.create_mutations_df()

        case_lists_data_path = self.src_dir_path / "msk_impact_2017_case_lists.csv"
        mutations_data_path = self.src_dir_path / "msk_impact_2017_mutations.csv"
        case_lists_df.to_csv(case_lists_data_path)
        mutations_df.to_csv(mutations_data_path)
        _logger.info("Successfully transformed cBioPortal data.")

    def create_case_lists_df(self) -> pd.DataFrame:
        """Create case lists data frame

        :return: Dataframe containing case list data
        """
        pathlist = Path(f"{self.msk_impact_2017_dir}/case_lists").glob("case_list_*")
        frames = []
        for p in pathlist:
            tmp = pd.read_csv(str(p), sep="delimiter", header=None)
            tmp = tmp[0].str.split(":", n=1, expand=True)
            tmp = tmp.set_index(0).T
            del tmp["cancer_study_identifier"]
            frames.append(tmp)
        return pd.concat(frames, ignore_index=True)

    def create_mutations_df(self) -> pd.DataFrame:
        """Create mutations data frame

        :return: Dataframe containing mutation data
        """
        return pd.read_csv(
            f"{self.msk_impact_2017_dir}/data_mutations.txt", sep="\t", skiprows=1
        )
