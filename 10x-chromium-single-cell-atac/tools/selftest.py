#!/usr/bin/env python3
"""Import-time validation for the 10x Chromium scATAC model."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import chromium_atac  # noqa: F401
from checks import Check

check = Check()
check.report()
