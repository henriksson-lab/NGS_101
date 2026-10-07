#!/usr/bin/env python3
"""Focused CG000184 source-transcription check."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_crispr as C
want="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAAGCAGTGGTATCAACGCAGAG"
if C.FEATURE_SI_REVERSE != want: raise SystemExit("FAIL: CG000184 Feature SI reverse primer transcription")
raise SystemExit(0)
