#!/usr/bin/env python3
"""Import-time source and construct validation for Drop-ChIP."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import dropchip  # noqa: F401
from checks import Check

check = Check()
check.report()
