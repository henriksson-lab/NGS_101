#!/usr/bin/env python3
"""Transcribe the active 2020 single-target oligos from published Table 4."""
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chil_chem as C
from checks import Check, Source, have
from storage import data_root
c=Check()
source=Source(data_root()/"sources"/H.parent.name/"s41596-020-0375-8_41596_2020_375_MOESM1_ESM.pdf.txt", "2020 published Supplementary Table 4; reacquire via tools/get_sources.py")
if have(c,source,label="published 2020 active oligos"):
    text=source.path.read_text()
    # The source table's sequence column is separate from names/modification columns.
    compact=re.sub(r"\s+", "", text).replace("[8-merindex]", "N"*8)
    for name, sequence in (("ChILv1-Fw",C.CHIL_V1_FW),("ChILv1-Rv",C.CHIL_V1_RV),
                           ("MEDSB-Fw",C.MEDSB_FW),("MEDSB-Rv",C.MEDSB_RV),
                           ("Read2+N11",C.READ2_N11),("Ad1_noMX",C.AD1_NOMX),
                           ("Ad2_index",C.AD2_INDEX)):
        c("2020 Table4 "+name,sequence in compact,True)
c.report()
