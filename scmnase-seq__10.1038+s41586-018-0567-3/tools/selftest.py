#!/usr/bin/env python3
"""Exercise the complete source-owned graph; no independent oligo transcription claimed."""
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import scmnase_seq as S
from checks import Check
S.protocol_graph()
Check().report()
