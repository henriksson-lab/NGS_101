#!/usr/bin/env python3
"""Source-transcription boundary for VASA-seq's external supplement."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))

from checks import Check

check = Check()
check.skip("Supplementary Table 12 transcription",
           "the source XLSX is external data and is not committed to the repository")
check.report()
