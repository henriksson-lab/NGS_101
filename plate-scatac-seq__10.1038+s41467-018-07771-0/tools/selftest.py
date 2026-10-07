#!/usr/bin/env python3
from __future__ import annotations
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import plate_scatac as P
from checks import Check
c=Check(); c.section("2018 Supplementary Methods transcription")
c("N701",P.N701,"CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGG")
c("S502",P.S502,"AATGATACGGCGACCACCGAGATCTACACCTCTCTATTCGTCGGCAGCGTC")
c.report()
