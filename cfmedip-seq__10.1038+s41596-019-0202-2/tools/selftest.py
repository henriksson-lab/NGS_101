#!/usr/bin/env python3
"""Regression for the native single-strand IP PCR boundary."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import chem
from checks import Check
c=Check()
g=chem.protocol_graph()
# The preceding IP material is not a completed duplex library. The unavailable
# geometry must preserve precisely that incoming Scene rather than manufacture it.
c("IP PCR geometry retains actual production aliquot",
  g.steps["library_pcr"].geometry.template is g.materials["production_template"].state,True)
c.report()
