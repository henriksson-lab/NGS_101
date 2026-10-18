#!/usr/bin/env python3
"""Import-time construct validation for 10x 5' V(D)J v3."""
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/"10x-chromium-single-cell-5-vdj"/"tools"),str(HERE.parents[1]/"lib")]
import vdj_v3  # noqa: F401
from checks import Check, Source, have
from storage import sources_root
check = Check()
source = Source(sources_root() / "10x-chromium-single-cell-5-vdj" /
                "CG000331_5v2_UserGuide_RevE.pdf.txt",
                "10x CG000331 Rev E Appendix p.72, cDNA Primers PN-2000089")
if have(check, source, label="catalog-linked cDNA primer transcription"):
    text = source.path.read_text(encoding="utf-8")
    block = text.split("cDNA Primers")[-1].split("Protocol Step 3.2")[0]
    for primer in vdj_v3.CDNA_PRIMERS:
        check("CG000331 Rev E PN-2000089 " + primer.name,
              "5’" + "-" + primer.sequence + "-3’" in block)
# Exercise the complete local declaration; import-time constructors alone do not
# cover the source-owned graph assembly.
vdj_v3.protocol_graph()
check.report()
