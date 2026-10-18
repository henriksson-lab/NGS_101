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
source = Source(sources_root() / HERE.parent.name / "s41587-022-01487-9.pdf.txt",
                "GLORI defining Liu2023 Methods, Library preparation and sequencing")
if have(check, source):
    text = "".join(source.path.read_text().split())
    check("source modified RNA linker", "5′rAPP-" + chem.LINKER.top + "-3SpC3" in text)
    check("source exact RT primer", "".join(s.top for s in chem.RT_PRIMER) in text)
    check("source modified UMI DNA adapter", "5′Phos-" + chem.UMI.top + chem.DNA_ADAPTER_ARM.top + "-3SpC3" in text)
    # Compare declared acquisition and analysis handling to the independent Methods,
    # rather than treating a processed/deposited read set as the acquired run.
    import re
    run = re.search(r"2×(\d+)basepairs", text)
    used = re.search(r"onlyreads(R[12])wasusedinourstudy", text)
    if run and used:
        plan = next(iter(chem.protocol_graph().sequencing.values()))
        source_reads = {"Read 1": int(run.group(1)), "Read 2": int(run.group(1))}
        check("acquired paired-end lengths from defining Methods", dict(plan.read_lengths), source_reads)
        analysed = "Read " + used.group(1)[1]
        check("source-used RNA-sense read retained for analysis", plan.read_uses[analysed], "extract_then_include")
        for read in set(source_reads) - {analysed}:
            check("other acquired read remains explicitly unused", plan.read_uses[read], "not_used")
    else:
        check("defining run/analysis passages still identifiable", False)
check.report()
