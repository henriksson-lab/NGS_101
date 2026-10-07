#!/usr/bin/env python3
"""Source-boundary check for LIANTI-specific sequence transcription."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import lianti as L
from checks import Check

check = Check()
check.section("secondary-source LIANTI transcription")
check("LIANTI transposon", L.TRANSPOSON,
      "AGATGTGTATAAGAGACAGGAACAGAATTTAATACGACTCACTATAGGGAGATGTGTATAAGAGACAG")
check.report()
