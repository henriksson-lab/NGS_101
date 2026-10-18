#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import dbit_seq as P
from checks import Check,Source,have
from storage import sources_root
c=Check();P.protocol_graph()
src=Source(sources_root()/H.parent.name/"PMC8132129_supp_mmc2.xlsx.txt","2021 STARProtocol original TableS2; not asserted as2020TableS2")
if have(c,src,label="2021 TableS2 exact implementation oligos"):
 rows={line.split("\t")[0]:line.split("\t")[1] for line in src.path.read_text().splitlines() if "\t" in line}
 for name,expected in (("Primer 1",P.B_PCR_HANDLE),("Primer 2",P.PRIMER_2),("Ligation Linker",P.LIGATION_SPLINT),("Template Switch Primer",P.TEMPLATE_SWITCH_PRIMER),("P5 oligo",P.P5_OLIGO),("P7 oligo (701)",P.P7_701),("Barcode A-1","/5Phos/"+"".join(s.top for s in P.a_oligo()))):
  c("2021 source "+name,rows[name],expected)
c.report()
