#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import r2c2 as R
from checks import Check
c=Check(); c.section("R2C2 source transcription")
c("TSO index is seven bases", len(R.CDNA.get("7-nt TSO index")), 7)
c("PCR index is eight bases", len(R.CDNA.get("8-nt Nextera A index")), 8)
c("page model shows repeated circular template", len(R.rca_construct()), 3 * len(R.CIRCLE_TEMPLATE))
c.report()
