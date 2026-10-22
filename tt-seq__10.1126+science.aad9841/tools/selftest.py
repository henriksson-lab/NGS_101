#!/usr/bin/env python3
"""TT-seq source run-layout check; geometry is enforced by construction."""
import re
from pathlib import Path
import build_page
from checks import Check, Source, have
from storage import sources_root

build_page.render(build_page.M)
check = Check()
source = Source(sources_root() / Path(__file__).resolve().parents[1].name /
    "aad9841-Schwalb-SM.pdf.txt", "Schwalb2016 published supplementary Methods; user supplied identical copy2026-10-22")
if have(check, source, label="published acquisition recipe"):
    text = source.path.read_text()
    match = re.search(r"Paired-end (\d+) base reads with additional (\d+) base reads of\s+barcodes", text)
    observed = ({"Read 1":int(match[1]), "Index 1 (i7)":int(match[2]), "Read 2":int(match[1])}
                if match else None)
    graph = build_page.M.protocol_graph()
    plan = graph.sequencing["library"]
    check("typed run lengths agree with published paired/index acquisition", dict(plan.read_lengths), observed)
check.report()
