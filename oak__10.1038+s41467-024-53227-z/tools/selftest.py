#!/usr/bin/env python3
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from checks import Check,Source,have
from storage import sources_root
c=Check()
import chem as M
M.protocol_graph()
source=Source(sources_root()/H.parent.name/"s41467-024-53227-z_PMC11499997.xml.txt","Defining2024 OAK published Methods")
if have(c,source):
    text="".join(source.path.read_text().split())
    for name,value in (("TSOrecognition",M.RNA_TSO_PCR),("aliquotI5",M.RNA_ALIQUOT_I5),("partialP5",M.RNA_PARTIAL_P5),("libraryI7",M.RNA_LIBRARY_I7)):
        c("published Methods "+name,value in text,True)
vendor=Source(sources_root()/"10x-chromium-single-cell-3-ge-v3"/"scg_CG000183_ChromiumSingleCell3__v3_UG_Rev-A.pdf.txt","10x original CG000183RevA Appendix54, inherited only3′v3 RNA-only bead component")
if have(c,vendor):
    text=vendor.path.read_text()
    c("sourcev3 bead barcode/UMI/T30",("CTACACGACGCTCTTCCGATCT-"+"N"*16+"-"+"N"*12+"-"+"T"*30) in text,True)
c.report()
