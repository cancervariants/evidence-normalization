"""Module for initializing ETL classes"""
from pathlib import Path
from os import environ

ETL_PATH = Path(__file__).resolve().parents[0]
ETL_DATA_DIR_PATH = environ.get("ETL_DATA_DIR_PATH", ETL_PATH / "data")
