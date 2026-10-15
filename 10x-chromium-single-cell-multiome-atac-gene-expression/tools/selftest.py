#!/usr/bin/env python3
from pathlib import Path
import os
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import multiome
from multimodal_page import render
from checks import Check, Source, have
from protocol_graph import SourceDiscrepancy

check = Check()
source_dir = Path(os.environ.get("CHEM_DATA", H.parents[1] / "_data")) / "sources" / H.parent.name
rev_b = Source(source_dir / "CG000338_ChromiumNextGEM_Multiome_ATAC_GEX_User_Guide_RevB.pdf.txt",
               "10x CG000338 Rev B, Appendix p. 80")
rev_g = Source(source_dir / "10x_CG000338_ChromiumNextGEM_Multiome_ATAC_GEX_User_Guide_RevG.pdf.pdf.txt",
               "10x CG000338 Rev G, Appendix pp. 80–81")
if have(check, rev_b, rev_g, label="Multiome guide revision discrepancy"):
    old = rev_b.path.read_text(encoding="utf-8")
    new = rev_g.path.read_text(encoding="utf-8")
    short_p5 = "AATGATACGGCGACCACCGAGATCT-N10"
    full_p5 = "AATGATACGGCGACCACCGAGATCTACAC-N10"
    check("Rev B GEX layout omits ACAC before i5", short_p5 in old)
    check("Rev G GEX layout includes ACAC before i5", full_p5 in new)
    check("model follows Rev G P5", multiome.GEX.get("P5").top + "-N10" == full_p5)
graph = multiome.protocol_graph()
check("guide discrepancy is attached to GEX library",
      any(isinstance(note, SourceDiscrepancy) and note.target == "gex_library"
          for note in graph.annotations))
page = render(multiome)
check("source discrepancy is visible on page",
      "CG000338 Rev B omits four P5 bases" in page)
check.report()
