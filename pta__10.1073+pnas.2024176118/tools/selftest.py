#!/usr/bin/env python3
from pathlib import Path
import sys
H = Path(__file__).resolve().parent
sys.path[:0] = [str(H), str(H.parents[1] / "lib")]
import chem as P
from checks import Check

c = Check()
c.section("PTA native population behavior")
g = P.protocol_graph()
c("PCR geometry uses incoming native ligation state",
  g.steps["library_pcr"].geometry.template is g.materials["ligated"].state, True)
c("no borrowed sequencing primer geometry", g.sequencing["selected"].primers, ())
c("terminated product records alpha-thio-ddN chemistry",
  "alpha-thio-ddN" in P.amplification_product().segments[-1].note, True)
c("independent original and daughter populations survive amplification model",
  len(g.materials["pta_product"].state.strands), 3)
c.report()
