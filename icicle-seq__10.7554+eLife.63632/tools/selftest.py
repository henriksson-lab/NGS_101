#!/usr/bin/env python3
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from storage import data_root
import illumina as il
import icicle_seq as M
from checks import Check, Source, have
c=Check()
bundle=data_root()/"sources"/"icicle-seq__10.7554+eLife.63632"
table=Source(bundle/"eLife.63632_elife-63632-supp5-v2.xlsx.txt",M.TABLE5)
if have(c,table,label="ICICLE Table 5 oligo transcription"):
    rows={line.split("\t")[0]:line.split("\t")[1] for line in table.path.read_text().splitlines() if "\t" in line}
    def bases(name):
        printed=re.sub(r"/[^/]*/","",rows[name])
        printed=re.sub(r"^5['’]?-|-['’]?3['’]?$","",printed.strip().replace(" ",""))
        return re.sub(r"[^ACGTNV]","",printed)
    c("MOSAIC_Bot",bases("MOSAIC_Bot"),M.MOSAIC_BOT)
    c("MOSAIC_Bot 5' phosphate printed","/5phos/" in rows["MOSAIC_Bot"].lower())
    c("Tn5ME-s7_Top",bases("Tn5ME-s7_Top"),M.TN5ME_S7_TOP)
    c("Poly-A Top-L",bases("Poly-A Top-L"),M.POLY_T_TOP+"VN"+M.nx.ME)
    c("SI-P5-22",bases("SI-P5-22"),M.SI_P5_22)
    c("F BC Primer",bases("F BC Primer"),M.F_BC_PRIMER)
    c("ADT-Rev-AMP",bases("ADT-Rev-AMP"),M.ADT_REV_AMP)
    c("TruSeq R1 Seq Primer",bases("TruSeq R1 Seq Primer"),M.TRUSEQ_R1_SEQ)
    c("ADT-i7 Primer",bases("ADT-i7 Primer"),il.P7+"N"*8+M.ADT_I7_ARM)
    c("ADT-i7 phosphorothioate before final A",rows["ADT-i7 Primer"].replace(" ","").endswith("CC*A-3'"))
c.report()
