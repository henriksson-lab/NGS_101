#!/usr/bin/env python3
"""There are no published vendor sequences to transcribe for SureCell WTA 3'."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
from checks import Check

check = Check()
check.skip("vendor oligo transcription",
           "Illumina Reference Guide #1000000021452 v01 prints reagent names and workflow, but no oligo sequences")
check.report()

