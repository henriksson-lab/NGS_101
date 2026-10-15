#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
SHARED = HERE.parents[1] / "scrb-seq__10.1101+003236" / "tools"
sys.path.insert(0, str(SHARED))
from selftest_common import run
from checks import Check
from chem import protocol_graph
import scrb

if __name__ == "__main__":
    graph = protocol_graph()
    check = Check()
    check("mcSCRB-seq source-owned graph is publishable", graph.publishable, True)
    check("mcSCRB-seq graph names its own defining source",
          graph.authority, "https://doi.org/10.1038/s41467-018-05347-6")
    check("mcSCRB-seq tagmentation waste is accounted for",
          graph.terminals["tagment_waste"].disposition.value, "discard")
    check.report()
    run(scrb.MCSCRB)

