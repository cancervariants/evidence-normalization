"""Provide utilities for tests."""

from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def evidence_data_dir() -> Path:
    """Return evidence data directory"""
    return Path(__file__).resolve().parents[1] / "evidence" / "data"
