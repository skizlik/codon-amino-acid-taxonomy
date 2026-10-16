"""Every number tagged <!--r:key--> in README.md and ERRATUM.md must equal results/results.json
at the printed precision.

Tag forms:  0.880<!--r:key-->   +0.112<!--r:key-->   12,963<!--r:key-->
            23%<!--r:key-->   (value x 100)        <0.001<!--r:key-->  (value below the bound)
"""
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
TAG = re.compile(r"(<)?([+\-−]?\d[\d,]*(?:\.\d+)?)(%)?<!--r:([A-Za-z0-9_.\-]+)-->")


def tags(files=("README.md", "ERRATUM.md")):
    return [t for f in files for t in TAG.findall((ROOT / f).read_text(encoding="utf-8"))]


def test_readme_quotes_tagged_numbers():
    assert len(tags(["README.md"])) >= 100
    assert len(tags(["ERRATUM.md"])) >= 5


@pytest.mark.parametrize("below,printed,percent,key", tags())
def test_number_matches_results(below, printed, percent, key):
    results = json.loads((ROOT / "results" / "results.json").read_text())
    assert key in results, f"unknown key {key}"
    value = results[key] * (100 if percent else 1)
    text = printed.replace("−", "-").replace(",", "")
    decimals = len(text.split(".")[1]) if "." in text else 0
    if below:
        assert abs(value) < float(text), f"{key}: {value} is not below {text}"
        return
    shown = f"{value:+.{decimals}f}" if text[0] in "+-" else f"{value:.{decimals}f}"
    assert shown == text, f"{key}: results.json has {value}, README prints {printed}"
