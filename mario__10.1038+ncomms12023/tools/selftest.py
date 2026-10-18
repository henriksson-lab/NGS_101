#!/usr/bin/env python3
"""Whole modified-oligo transcription against the defining published article."""
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import mario_chem as m
from checks import Check, Source, have
from storage import sources_root

c = Check()
source = Source(sources_root() / HERE.parent.name / "ncomms12023_PMC4931010.xml.txt",
                "Nguyen2016 Methods, Oligonucleotide sequences used in MARIO")
if have(c, source, label="published modified MARIO oligos"):
    text = re.sub(r"\s+", "", source.path.read_text()).upper()
    written = (
        ("RNA linker", "".join("r" + b for b in m.RNA_LINKER_LEFT) + "/iBiodT/" + "".join("r" + b for b in m.RNA_LINKER_RIGHT)),
        ("guard", "T*C*G*C*" + m.LINKER_CDNA[4:]),
        ("RT adapter", "/5rApp/" + m.RT_ADAPTER + "/3ddC/"),
        ("ES-1 RT", "/5Phos/NNAGGTNNN" + m.RT_BODY),
        ("ES-2 RT", "/5Phos/NNCGCCNNNN" + m.RT_BODY),
        ("ES-indirect RT", "/5Phos/NNCATTNNNN" + m.RT_BODY),
        ("Cut oligo", m.CUT_OLIGO + "/3InvdT/"),
        ("DP5", m.DP5), ("DP3", m.DP3), ("PE1", m.PE1), ("PE2", m.PE2),
    )
    for name, oligo in written:
        c(name + " whole source transcription", oligo.upper() in text, True)
c.report()
