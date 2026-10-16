#!/usr/bin/env python3
from pathlib import Path
import os
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import flex
from checks import Check, Source, have
from multimodal_page import render
render(flex)
graph = flex.protocol_graph()
assert graph._producer["occupied"] == "partition"
assert graph.steps["seal_probe_pair"].inputs[0].material == "occupied"
assert graph.steps["barcode"].inputs[0].material == "reporter"
assert graph.steps["barcode_protein"].inputs[0].material == "protein_joined"
assert graph._producer["protein_joined"] == "attach_ab"
assert len(graph.sequencing) == 2
print("PASS  adjacent-probe and final-primer geometry build")
print("PASS  original-guide order: GEM partition, nick sealing, barcode transfer")
print("PASS  separate antibody capture and protein sequencing branch")
check = Check()
source_dir = Path(os.environ.get("CHEM_DATA", H.parents[1] / "_data")) / "sources" / H.parent.name
labeling = Source(source_dir / "10x_CG000529_DemonstratedProtocol_CellSurfaceProteinLabelingforFixedRNA_RevC.pdf.txt",
                  "10x CG000529 Rev C p. 4")
guide = Source(source_dir / "10x_CG000674_ChromiumFixedRNAProfiling_SingleplexedSamples_ProteinExpression_BarcodeOligoCapture_UserGuide_Rev_C.pdf.txt",
               "10x CG000674 Rev C Appendix p. 85")
if have(check, labeling, guide, label="Flex antibody and protein library source transcription"):
    written = "/5AmMC12/" + "-".join(
        f"N{len(segment.top)}" if segment.placeholder else segment.top
        for segment in flex.ANTIBODY_OLIGO.segments)
    check("source prints the modeled TotalSeq-C antibody oligo pattern",
          written in labeling.path.read_text(encoding="utf-8"))
    text = guide.path.read_text(encoding="utf-8")
    for name in ("protein constant 1", "protein constant 2", "Nextera Read 2 reverse-complement arm"):
        segment = next(s for s in flex.PROTEIN_LIB.segments if s.name == name)
        check(f"source prints modeled {name}", segment.top in text)
check.report()
