#!/usr/bin/env python3
"""Construct the local source-qualified graph and check original transcriptions."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem
from checks import Check, Source, have
from storage import sources_root
chem.protocol_graph()
check = Check()
source = Source(sources_root() / HERE.parent.name / "nature13802_PMC4224642.xml.txt",
                "Pseudo-seq published author manuscript PMC4224642 Library Preparation")
if have(check, source):
    text = source.path.read_text()
    for oligo in chem.oligos():
        name, modified = oligo.split(" — 5′-", 1)
        modified = modified.removesuffix("-3′")
        check(name + " exact modified source string", modified in text)
check.report()
