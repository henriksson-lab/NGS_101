#!/usr/bin/env python3
"""Import-time validation for the Seq-Well S3 construction model."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import seqwell_s3  # noqa: F401
from checks import Check

check = Check()
check.report()
