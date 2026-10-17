#!/usr/bin/env python3
"""Source-transcription checks for HiRES.

GAT5 and RNA-P7 are exact graph states covered by tools/selftest_source_oligos.py.
GAT5-RT and GAT5-7N carry N windows, so they are checked here against the
Bio-protocol materials list (item 43), with its PDF line wrap rejoined.
"""
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from storage import data_root
import chem as P
from checks import Check, Source, have

c = Check()
protocol = Source(data_root() / "sources" / H.parent.name / "bio-protocol4886.pdf.txt",
                  "https://doi.org/10.21769/BioProtoc.4886 (PDF, converted by tools/doctext.py)")
if have(c, protocol, label="HiRES Bio-protocol primer transcription"):
    text = protocol.path.read_text(encoding="utf-8")
    start = text.index("43. Primers")
    block = text[start:text.index("Solutions", start)]
    printed = {}
    for name, seq in re.findall(r"^([A-Z0-9-]+): ([ACGTNVBDHRYKMSW \n]+?)(?=^[A-Z0-9-]+:|\Z)",
                                block, re.M):
        printed[name] = re.sub(r"\s+", "", seq)
    c("GAT5-RT matches Bio-protocol item 43", P.GAT5_RT, printed.get("GAT5-RT"))
    c("GAT5-7N matches Bio-protocol item 43", P.GAT5_7N, printed.get("GAT5-7N"))
c.report()
