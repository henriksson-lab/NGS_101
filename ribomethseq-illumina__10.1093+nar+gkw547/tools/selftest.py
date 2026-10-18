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
source = Source(sources_root() / HERE.parent.name / "nar_gkw547_supp_supp_gkw547_nar-00545-met-k-2016-File008.doc.txt",
                "Marchand2016 original supplementary Table S2 run-type key")
if have(check, source):
    import re
    text = "".join(source.path.read_text().split())
    key = re.search(r"SR(\d+)=singleread(\d+)nt,PE2\*(\d+)=paired-end2times(\d+)nt", text)
    if key:
        plans = chem.protocol_graph().sequencing
        check("SR acquisition from the published Table S2 key", dict(plans["library_sr50"].read_lengths),
              {"Read 1": int(key.group(2))})
        check("paired acquisition from the published Table S2 key", dict(plans["library_pe75"].read_lengths),
              {"Read 1": int(key.group(4)), "Read 2": int(key.group(4))})
    else:
        check("published Table S2 run key remains identifiable", False)
check.report()
