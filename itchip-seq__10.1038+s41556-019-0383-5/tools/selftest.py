#!/usr/bin/env python3
"""Import-time source transcription and construction validation for itChIP-seq."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]

import itchip  # noqa: F401  (import validates source assemblies and construct interlocks)
from checks import Check

check = Check()
check.report()
