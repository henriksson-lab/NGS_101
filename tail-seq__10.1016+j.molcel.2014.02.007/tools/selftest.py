#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import tail_seq as M
from checks import Check, Source, have
from storage import data_root
from chemdraw import revcomp
c=Check();c("tail-facing read is long",M.READ_LENGTHS["Read 2"],251)
geo=Source(data_root()/"sources"/H.parent.name/"geo_GSM1242324_self_full.txt",
           "GEO GSM1242324 SOFT text: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1242324&targ=self&form=text&view=full")
if have(c,geo,label="GEO protocol transcription"):
    t=geo.path.read_text(encoding="utf-8",errors="replace")
    lib={s.name:s for s in M.FINAL_LIBRARY}
    c("GEO prints read-2 delimiter GTCAG","delimiter sequence (‘GTCAG’ as in the direction of read 2)" in t,True)
    c("modelled delimiter is GTCAG in read-2 direction",revcomp(lib["3′-adapter delimiter"].top),"GTCAG")
    c("GEO: degenerate 3′-adapter bases at read-2 cycles 1-15","1st to 15th cycle in read 2 (degenerate bases in 3′ adapter)" in t,True)
    c("modelled degenerate window is 15 nt",len(lib["3′-adapter degenerate bases"].top),15)
    c("GEO read 1 is 51 nt","full read 1 sequences which are 51 nt long" in t,True)
    c("modelled read 1 is 51 nt",M.READ_LENGTHS["Read 1"],51)
    c("GEO order: Ribo-zero, 3′ adapter, RNase T1, streptavidin, PNK, 500-1000 nt gel, 5′ adapter, RT, PCR",
      "depleted of rRNA by using Ribo-zero kit (Epicentre, MRZH11124), ligated to 3′ adapter and partially digested by RNase T1" in t and "gel purified in the range of 500 – 1000 nucleotides. The purified RNAs were ligated to 5′ adapter" in t,True)
c.report()
