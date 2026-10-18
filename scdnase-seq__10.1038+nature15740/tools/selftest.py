#!/usr/bin/env python3
"""Check the obsolete vendor oligo transcription against its original."""
import re
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
from checks import Check,Source,have
from storage import sources_root
import scdnase_seq as S
check=Check()
source=Source(sources_root()/"scdnase-seq__10.1038+nature15740"/
    "illumina_1000000002694_22_ilmn-adapter-sequences.pdf.txt",
    "Illumina Adapter Sequences document1000000002694 v22, obsolete Multiplexing kit pp100–101")
if have(check,source,label="obsolete Multiplexing oligos"):
    text=source.path.read_text(encoding="utf-8")
    block=text.split("Oligonucleotide Sequences for the Multiplexing Sample",1)[1]
    adapter=re.search(r"Multiplexing Adapters\s+5′ P-([ACGT]+)\s+5′ ([ACGT]+)",block)
    check("vendor Multiplexing adapter strands",tuple(adapter.groups()) if adapter else (),
        ("".join(s.top for s in S.adapter_top()),"".join(s.top for s in S.adapter_bottom())))
    primer=re.search(r"Multiplexing PCR Primer 1.0\s+5′ ([ACGT]+)",block)
    check("vendor Multiplexing Primer1.0",primer.group(1) if primer else None,
        "".join(s.top for s in S.p5_primer()))
    members=re.findall(r"PCR Primer, Index \d+\s+5′ ([ACGT]+)",block)
    architecture="".join(s.top for s in S.indexed_p7_primer())
    pattern=architecture.replace("N","[ACGT]")
    check("all vendor indexed primers match the declared six-base/13-base-end architecture",
        bool(members) and all(re.fullmatch(pattern,x) is not None for x in members))
S.protocol_graph()
check.report()
