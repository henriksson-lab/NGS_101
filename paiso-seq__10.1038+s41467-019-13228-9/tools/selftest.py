#!/usr/bin/env python3
"""Source-transcription checks for PAIso-seq."""
from pathlib import Path
import sys
import os
import re
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check,Source,have
from dumbbell import smrtbell_hairpin_geometry
c=Check()
c("TSO-derived handle",P.INSERT.segments[0].top,"AAGCAGTGGTATCAACGCAGAGT")
# Supplementary Table 1, "Oligos for PAIso-seq library construction", rows TSO / RT primer
# / PCR primer.  The bases come from smart-seq2__10.1038+nmeth.2639/tools/smartseq2; these checks pin them to the table, so
# the ACAT TSO linker drawn on this page stays the one PAIso-seq itself published.
c("published TSO ordering string",P.TSO_ORDERING,
  "AAGCAGTGGTATCAACGCAGAGTACATrGrG+G")
c("published RT primer",P.RT_PRIMER,"AAGCAGTGGTATCAACGCAGAGTAC")
c("published PCR primer",P.PCR_PRIMER,"AAGCAGTGGTATCAACGCAGAGT")
c("the published extension family has fifteen members",len(P.BARCODES),15)
c("bulk BC3 ordering string",P.extension_oligo(3),
  "AAGCAGTGGTATCAACGCAGAGTACTACTAGAGTAGCACTCdUTTTTTTTTdUTTTTTTTT")
c("unversioned SMRTbell junction stays unresolved",P.SMRTBELL.junction_unavailable,True)
c("ligating end cannot be presented as a known A/T junction",
  smrtbell_hairpin_geometry(P.SMRTBELL,statement="Kit version unresolved").availability,
  "unavailable")
bundle=Path(os.environ.get("CHEM_DATA",H.parents[1]/"_data"))/"sources"/"paiso-seq__10.1038+s41467-019-13228-9"
source=Source(bundle/"s41467-019-13228-9_supp_41467_2019_13228_MOESM1_ESM.pdf.txt",P.SOURCE)
if have(c,source,label="original Supplementary Table1"):
    table=source.path.read_text().split("Supplementary Table 1.",1)[1].split("Supplementary Table 2.",1)[0]
    compact=re.sub(r"\s+","",table)
    for name,seq in (("TSO",P.TSO_ORDERING),("RT primer",P.RT_PRIMER),("PCR primer",P.PCR_PRIMER)):
        c(name+" transcribed from original table",seq in compact)
    c("all fifteen whole modified extension oligos occur in original table",all(P.extension_oligo(i) in compact for i in range(1,16)))
c.report()
