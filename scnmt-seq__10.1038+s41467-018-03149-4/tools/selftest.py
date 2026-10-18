#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import scnmt_seq
from checks import Check
# Native source-specific strands remain unassigned; graph construction enforces
# molecular continuity and no unverified predecessor sequence is used as a test oracle.
Check().report()
