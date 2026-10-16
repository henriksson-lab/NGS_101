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
order = tuple(graph.steps)
required = ("reverse_transcription", "remove_free_primers", "pool",
            "second_strand", "ivt")
check("authors' 2.0 Exo I and pooling order",
      tuple(order.index(s) for s in required),
      tuple(sorted(order.index(s) for s in required)))
check("2.0 Exo I fragments have a discard terminal",
      graph.terminals["exo_primer_fragments"].disposition.value, "discard")
check.report()

runpy.run_path(str(HERE.parents[1] / "mars-seq__10.1126+science.1247651" /
                   "tools" / "selftest.py"), run_name="__main__")
