"""Module for initializing ETL classes"""
from pathlib import Path

ETL_PATH = Path(__file__).resolve().parents[0]
ETL_DATA_DIR_PATH = ETL_PATH / "data"
