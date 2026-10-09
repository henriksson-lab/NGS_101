#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import scnmt_seq
from checks import Check
# The defining paper delegates its two branches to published Smart-seq2/Nextera XT and
# scBS-seq chemistries; their source transcription is checked in those protocol suites.
Check().report()
