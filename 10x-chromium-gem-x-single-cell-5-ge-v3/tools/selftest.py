#!/usr/bin/env python3
"""Check complete source-printed endpoints and catalog-linked primers against originals."""
from pathlib import Path
import os
import re
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import tenx_5p_v3 as V
from checks import Check, Source, have
check = Check()
source_root = Path(os.environ.get("CHEM_DATA", HERE.parents[1] / "_data")) / "sources" / HERE.parent.name

def normalise(s):
    s = s.replace("rG", "G").replace("cDNA_Insert", "CDNAINSERT")
    s = re.sub(r"([NT])(\d+)", lambda m: m[1] * int(m[2]), s)
    return re.sub(r"[^A-Z]", "", s)

def model(parts):
    s = "".join(p.top for p in parts).replace("B", "N").replace("U", "N").replace("J", "N").replace("I", "N")
    return re.sub(r"X+", "CDNAINSERT", s)

original = Source(source_root / "CG000733_RevA.pdf.txt", "official defining 10x guide; reacquire via tools/get_sources.py")
if have(check, original, label="original guide source transcription"):
    lines = original.path.read_text().splitlines()
    bead = next(l for l in lines if l.startswith("5’-") and l.startswith("5’-" + V.PARTIAL_R1))
    endpoint = next(l for l in lines if l.startswith("5’-AATGAT") and "cDNA_Insert" in l)
    check("complete source-printed bead oligo", model(V.bead_tso()), normalise(bead))
    check("complete source-printed final GEX library", model(V.final_library()), normalise(endpoint))
predecessor = Source(source_root / "CG000331_RevE_predecessor.pdf.txt", "official predecessor guide; same cDNA Primers PN-2000089")
if have(check, predecessor, label="catalog-linked cDNA primer originals"):
    text = predecessor.path.read_text()
    for primer in V.cdna_primers():
        check(primer.name + " source transcription", "5’-" + primer.top() + "-3’" in text)
check.report()
