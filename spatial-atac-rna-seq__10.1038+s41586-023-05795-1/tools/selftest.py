#!/usr/bin/env python3
"""Transcription of the original 2023 modified RT and library oligos.

The ordinary source guard checks exact unmodified whole strands. This suite also
checks the internal biotin nucleotide and source TSO modification notation.
"""
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check, Source, have
from storage import sources_root
c=Check()
bundle=sources_root()/"spatial-atac-rna-seq__10.1038+s41586-023-05795-1"
table=Source(bundle/"s41586-023-05795-1_supp_41586_2023_5795_MOESM1_ESM.pdf.txt", "Published2023 SupplementaryInformation, Table4")
if have(c,table,label="2023 oligo table"):
    text=table.path.read_text()
    block=text[text.index("RT primer"):text.index("Supplementary Table 5.")]
    def value(name,next_name):
        return re.sub(r"\s+","",block.split(name,1)[1].split(next_name,1)[0])
    c("RT physical internal biotin-dT is a nucleotide in the printed modified oligo",
      "/5Phos/"+P.RT_HANDLE+P.RT_UMI+"/iBiodT/"+P.RT_POLYT,
      value("RT primer","Ligation linker 1"))
    c("first ligation linker",P.LINKER1,value("Ligation linker 1","Ligation linker 2"))
    c("second ligation linker",P.LINKER2,value("Ligation linker 2","PCR Primer 1"))
    c("RNA amplification primer1",P.B_PREFIX,value("PCR Primer 1","PCR Primer 2"))
    c("RNA amplification primer2",P.PCR2,value("PCR Primer 2","N501"))
    c("N501 remains a source-printed reference member",P.N501,value("N501","N701"))
    c("N701 remains a source-printed reference member",P.N701,value("N701","N702"))
main=Source(bundle/"s41586-023-05795-1_PMC10076218.xml.txt", "Published2023 Methods, on-bead template switching")
if have(c,main,label="2023 modified TSO"):
    text=main.path.read_text()
    printed=re.search(r"template switch primer \(([^)]+)\)",text)
    c("source terminal RNA/RNA/LNA-G notation",P.TSO_WRITTEN,printed.group(1) if printed else None)
c.report()
