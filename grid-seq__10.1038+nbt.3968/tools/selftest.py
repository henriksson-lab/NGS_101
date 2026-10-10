#!/usr/bin/env python3
import build_page
from checks import Check
c=Check()
c("published GRID DNA linker", build_page.M.LINKER_DNA, "GTTGGAGTTCGGTGTGTGGGAGTGAGCTGTGTC")
c("published mixed bridge structure",
  (build_page.M.BRIDGE.rna_fixed_5, build_page.M.BRIDGE.barcode_length,
   build_page.M.BRIDGE.rna_terminal, build_page.M.BRIDGE.biotin_offset,
   build_page.M.BRIDGE.preadenylated_rna_5),
  ("GUUGGAUUC", 3, "G", 8, True))
build_page.render(build_page.M)
c.report()
