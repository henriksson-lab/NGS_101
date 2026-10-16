#!/usr/bin/env python3
"""OAK source and branch-order regressions."""
from pathlib import Path
import sys

H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from storage import data_root
import chem as P
from checks import Check, Source, have

check=Check()
graph=P.protocol_graph()
order=list(graph.steps)
check("Multiome transposition precedes GEM partitioning",
      order.index("tagment") < order.index("multiome_droplet_partition"))
check("RNA-only branch uses untransposed fixed cells",
      graph.steps["rna_droplet_partition"].inputs[0].material == "fixed_rna")

source_dir=data_root()/"sources"
paper=Source(source_dir/H.parent.name/"s41467-024-53227-z_PMC11499997.xml.txt",
             "Wu et al. OAK original Methods")
if have(check,paper,label="OAK defining Methods"):
    text=paper.path.read_text(encoding="utf-8")
    for name,sequence in (("TSO recognition",P.RNA_TSO_PCR),
                          ("aliquot i5",P.RNA_ALIQUOT_I5),
                          ("partial P5",P.RNA_PARTIAL_P5),
                          ("library i7",P.RNA_LIBRARY_I7)):
        check(name+" source transcription",sequence in text)

vendor=Source(source_dir/"10x-chromium-single-cell-3-ge-v3"/
              "scg_CG000183_ChromiumSingleCell3__v3_UG_Rev-A.pdf.txt",
              "10x CG000183 Rev A original v3 oligo appendix")
if have(check,vendor,label="OAK inherited 10x v3 bead layout"):
    text=vendor.path.read_text(encoding="utf-8")
    bead="CTACACGACGCTCTTCCGATCT-"+"N"*16+"-"+"N"*12+"-"+"T"*30
    check("vendor v3 bead has 16-base barcode and 12-base UMI",bead in text)
    umi=next(segment for segment in P.RNA.segments if segment.name=="UMI")
    check("OAK RNA UMI follows vendor v3 bead",len(umi),12)
check.report()
