#!/usr/bin/env python3
"""Import-time construction validation for the scRRBS model."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import scrrbs  # noqa: F401  (construction validates adaptor and primer interlocks)
from checks import Check

check = Check()
check.report()
