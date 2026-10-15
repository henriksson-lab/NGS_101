#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
SHARED = HERE.parents[1] / "strt-seq__10.1101+gr.110882.110" / "tools"
sys.path.insert(0, str(SHARED))
from selftest_common import run
from checks import Check
from chem import protocol_graph
import strt

if __name__ == "__main__":
    graph = protocol_graph()
    check = Check()
    check("STRT-C1 source-owned graph is publishable", graph.publishable, True)
    check("STRT-C1 names its own defining source",
          graph.authority, "https://doi.org/10.1038/nmeth.2772")
    check("Unpublished C1 selection stays inferred",
          graph.materials["selected_five_prime_fragments"].inferred, True)
    check.report()
    run(strt.C1)

