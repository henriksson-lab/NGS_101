#!/usr/bin/env python3
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "lib")]
from checks import Check
from spatial_protocols import HDST_PRECURSOR_WRITTEN
c = Check(); c.section("HDST source transcription")
c("published bead precursor", HDST_PRECURSOR_WRITTEN,
  "/AmC6/UUUUUGACTCGTAATACGACTCACTATAGGGACACGACGCTCTTCCGATCT-Spatial_barcode_Pool1-Bridge1")
c.report()
