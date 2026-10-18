#!/usr/bin/env python3
"""Published oligo transcription, including physical capture LNA notation."""
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check,Source,have
from storage import sources_root
c=Check();bundle=sources_root()/H.parent.name
table=Source(bundle/"s12864-023-09832-1_supp_12864_2023_9832_MOESM2_ESM.xlsx.txt","Published2023 Additionalfile2 SupplementaryTable4")
if have(c,table,label="published transferred oligos"):
    rows={x[0]:x[1].strip() for line in table.path.read_text().splitlines() if len(x:=line.split("\t"))>=2}
    c("N5 representative transferred strand",P.P5_I5_1_TRANSFERRED,rows.get("P5_i5_1_Universal_Connector_A_C15_ME"))
    c("N7 representative transferred strand",P.P7_I7_1_TRANSFERRED,rows.get("P7_i7_1_Universal_Connector_B_D15_ME"))
    c("nontransfer physical5P", "/5Phos/"+P.TN5ME_REV,rows.get("Tn5ME_Rev"))
main=Source(bundle/"s12864-023-09832-1_PMC10696879.xml.txt","Published2023 Methods179 andSequencing249–261")
if have(c,main,label="published modified capture and sequencing primers"):
    text=re.sub(r"\s+","",main.path.read_text())
    capture=re.search(r"TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNXXXXXXXXXX([T+]+)VN",text)
    c("LNA marks modify the following T; nine sites after firstT3",P.CAPTURE_T30_WRITTEN,capture.group(1) if capture else None)
    for name,seq in [("R1",P.CUSTOM_R1),("R2",P.CUSTOM_R2),("I1",P.CUSTOM_I1),("I2",P.CUSTOM_I2)]:
        c("published custom "+name,seq in text,True)
c.report()
