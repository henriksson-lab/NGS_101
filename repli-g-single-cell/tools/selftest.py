#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
from checks import Check
from protocol_graph import Disposition, StepKind
import repli_g_single_cell as model

c=Check()
c.section("REPLI-g Single Cell native graph")
g=model.protocol_graph()
c("finished graph is publishable", g.publishable, True)
c("library stage is an unresolved boundary, not a canonical Nextera reconstruction",
  g.steps["library_preparation_boundary"].kind, StepKind.MAGIC)
c("boundary reason names the unavailable cited revisions",
  all(k in g.steps["library_preparation_boundary"].reason
      for k in ("August 2014", "15027987 Rev B", "Advanced DNA Single Cell")), True)
c("no borrowed sequencing primers or read lengths",
  (g.sequencing["library"].primers, g.sequencing["library"].read_lengths or {}), ((), {}))
def upstream(material):
    seen, stack = set(), [material]
    while stack:
        step = g.producer_of(stack.pop())
        for x in (step.inputs if step else ()):
            if x.material not in seen:
                seen.add(x.material); stack.append(x.material)
    return seen
lib_up = upstream("library")
c("library descends from the cell genome via the MDA aliquot",
  {"cell", "mda_product", "library_input"} <= lib_up, True)
c("no-template primer-dimer DNA and the stored portion never feed the library",
  lib_up & {"primer_dimer_dna", "primer_dimer_seed", "stored_dna", "no_template_control"}, set())
c("stored remainder is a physical stored terminal",
  g.terminals["stored_dna"].disposition, Disposition.STORED)
c.report()
