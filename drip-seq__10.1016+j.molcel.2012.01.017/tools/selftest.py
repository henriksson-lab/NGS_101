#!/usr/bin/env python3
"""Build the source-owned graph/page; no independent base transcription is claimed."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem
from batch_page import render
from checks import Check
chem.protocol_graph()
render(chem)
Check().report()
