#!/usr/bin/env python3
"""Construct the protocol-owned declaration and its graph invariants."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem
chem.workflow()

chem.protocol_graph()

from checks import Check, Source, have
from storage import sources_root
check = Check()
source = Source(sources_root() / HERE.parent.name / "PMC8763625_PMC8763625.xml.txt",
                "Cited ASAP-seq original PMC8763625 Methods, transposition and barcoding")
if have(check, source):
    text = source.path.read_text()
    for name in ("BOA", "BOB"):
        check(name + " exact blocked bridge in cited original",
              getattr(chem, "PHOSPHO_" + name) + "/3InvdT/" in text)
check.report()
