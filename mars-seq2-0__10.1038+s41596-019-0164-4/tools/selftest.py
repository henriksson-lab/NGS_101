#!/usr/bin/env python3
import runpy
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from checks import Check
from chem import protocol_graph

graph = protocol_graph()
check = Check()
check("MARS-seq2.0 source-owned graph is publishable", graph.publishable, True)
check("MARS-seq2.0 graph names its own defining source",
      graph.authority, "https://doi.org/10.1038/s41596-019-0164-4")
check("Instrument sequencing-primer oligos remain source-limited",
      bool(graph.sequencing["library"].limitation), True)
check.report()

runpy.run_path(str(HERE.parents[1] / "mars-seq__10.1126+science.1247651" /
                   "tools" / "selftest.py"), run_name="__main__")
