#!/usr/bin/env python3
"""Guard the reaction order read from Aw et al. JoVE 2017 sections 5–9."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]

from checks import Check
import chem

graph = chem.protocol_graph()
steps = list(graph.steps)
source_order = (
    "repair_three_prime", "phosphorylate_five_prime", "proximity_ligation",
    "photoreverse", "ligate_rna_linker", "library_rt", "circularize_cdna",
    "library_pcr",
)
check = Check()
check("JoVE end prep, chimera, linker, RT and circle precede PCR",
      tuple(step for step in steps if step in source_order), source_order)
check.report()
