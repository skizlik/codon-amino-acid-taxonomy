"""Make notebooks/codon_utils.py importable from the tests, and share the data."""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "notebooks"))

import codon_utils  # noqa: E402


@pytest.fixture(scope="session")
def raw():
    return codon_utils.read_raw()


@pytest.fixture(scope="session")
def clean(raw):
    df, log = codon_utils.filter_classes(raw)
    return codon_utils.add_groups(df), log
