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
source = Source(sources_root() / HERE.parent.name / "nature12894_PMC3966492.xml.txt",
                "DMS-seq defining original PMC3966492 Library Generation")
if have(check, source):
    raw = "/5rApp/" + chem.DMS_LINKER.top + "/3ddC/"
    check("exact modified linker in the defining Methods", raw in source.path.read_text())
check.report()
