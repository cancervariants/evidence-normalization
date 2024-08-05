"""The VICC library for normalizing evidence"""
from pathlib import Path
from os import environ

APP_ROOT = Path(__file__).resolve().parents[0]
DATA_DIR_PATH = environ.get("DATA_DIR_PATH", APP_ROOT / "data")
