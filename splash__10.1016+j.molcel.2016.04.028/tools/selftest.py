#!/usr/bin/env python3
"""Guard the reaction order read from Aw et al. JoVE 2017 sections 5–9."""
import sys
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from storage import data_root

from checks import Check, Source, have
import chem

graph = chem.protocol_graph()
steps = list(graph.steps)
source_order = (
    "repair_three_prime", "phosphorylate_five_prime", "proximity_ligation",
    "photoreverse", "ligate_rna_linker", "library_rt", "circularize_cdna",
    "library_pcr",
)
check = Check()
check("JoVE end prep, chimera, linker, RT and circle precede PCR",
      tuple(step for step in steps if step in source_order), source_order)
bundle=data_root()/"sources"/"splash__10.1016+j.molcel.2016.04.028"
source=Source(bundle/"1-s2.0-S1097276516301046-mmc1.pdf.txt",chem.DEFINITION.doi)
if have(check,source,label="original 2016 supplemental oligos"):
    text=source.path.read_text()
    check("original Supplemental Experimental Procedures adapter",re.search(r"3’ adapter primer sequence\. 5'-([ACGT]+)-3'",text)[1],chem.RNA_LINKER)
    from splash_oligos import RT_LEFT,RT_MIDDLE,RT_RIGHT
    primer=re.search(r"5'AGATCGGAAGAGC.*?-3'",text,re.S)[0]
    check("original two-spacer RT primer",re.sub(r"\s+","",primer)[2:-3],RT_LEFT+"/iSp18/"+RT_MIDDLE+"/iSp18/"+RT_RIGHT)
check.report()
