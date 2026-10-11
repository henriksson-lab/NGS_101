#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H.parents[1]/"lib")]
from checks import Check
from rare_cell_and_affinity import CAP_HUMAN_WASH_MM, CAP_MOUSE_WASH_MM
c=Check(); c.section("CAP-seq source transcription")
c("human / mouse NaCl wash concentrations (mM)",
  (CAP_HUMAN_WASH_MM, CAP_MOUSE_WASH_MM), (600, 560))
c.report()
