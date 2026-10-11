#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "lib")]
from checks import Check
from dsb_capture import get

c = Check()
# Construction itself validates the declared sequencing-primer landing sites.  There
# are no method-specific oligo constants to re-transcribe: Digenome-seq deliberately
# uses a commercial WGS library shell.
get("digenome-seq")
c.report()
