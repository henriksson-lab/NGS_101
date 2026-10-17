#!/usr/bin/env python3
import build_page
from checks import Check
c=Check()
c("published GRID DNA linker", build_page.M.LINKER_DNA, "GTTGGAGTTCGGTGTGTGGGAGTGAGCTGTGTC")
c("published mixed bridge structure",
  (build_page.M.BRIDGE.rna_fixed_5, build_page.M.BRIDGE.barcode_length,
   build_page.M.BRIDGE.rna_terminal, build_page.M.BRIDGE.biotin_offset,
   build_page.M.BRIDGE.preadenylated_rna_5),
  ("GUUGGAUUC", 3, "", 8, True))
# 2019 protocol (PMC7721247, Introduction): "45-bp linker and 20-21-nt RNA and DNA at each end".
c("linker length matches the 2019 protocol's 45 bp", build_page.M.LINKER_LENGTH, 45)
# 2017 Methods: the final PCR band recovered from the gel is 194 bp; tags of 20-21 nt.
c("published 194-bp band within the 20-21-nt tag range",
  len(build_page.M.FINAL_LIBRARY) <= 194 <= len(build_page.M.FINAL_LIBRARY) + 2, True)
build_page.render(build_page.M)
c.report()
