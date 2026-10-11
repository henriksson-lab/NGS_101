#!/usr/bin/env python3
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "lib")]
from checks import Check
from dsb_capture import BLESS_PROXIMAL, BLESS_DISTAL
c = Check()
c.section("BLESS Supplementary Table 5 transcription")
c("L1 proximal hairpin", BLESS_PROXIMAL, "TACTACCTCGAGAGTTACGCTAGGGATAACAGGGTAATATAGTTTTTTTCTATATTACCCTGTTATCCCTAGCGTAACTCTCGAGGTAGTA")
c("L3 distal hairpin", BLESS_DISTAL, "CGTCGTCTCGAGAGTTACGCTAGGGATAACAGGGTAATATAGTTTTTTTCTATATTACCCTGTTATCCCTAGCGTAACTCTCGAGACGACG")
c.report()
