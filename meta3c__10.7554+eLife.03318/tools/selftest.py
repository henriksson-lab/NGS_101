#!/usr/bin/env python3
"""Construction enforces graph provenance; exact adapter transcription uses the source guard."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem
chem.protocol_graph()
print("PASS meta3C graph construction; source transcription checked separately")
