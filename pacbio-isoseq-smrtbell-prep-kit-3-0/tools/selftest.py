#!/usr/bin/env python3
import build_page
build_page.render(build_page.M)
assert build_page.M.barcoded_primer(1) == "CTACACGACGCTCTTCCGATCTACTACACGCAATGAAGTCGCAGGGTTGGG"
assert len(build_page.M.ISOSEQ_BARCODES) == 12
print("PASS  Iso-Seq SMRTbell topology and pore geometry")
