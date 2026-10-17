#!/usr/bin/env python3
"""Source-transcription checks for scSNV-seq."""
import re
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from storage import data_root
import chem as P  # constructors validate identifier metadata and primer sites
from checks import Check, Source, have
c=Check()
table=Source(data_root()/"sources"/H.parent.name/
             "s13059-024-03169-y_supp_13059_2024_3169_MOESM3_ESM.xlsx.txt",
             "Cooper et al. 2024 Additional file 3: Table S2; fetch with tools/get_sources.py")
if have(c, table, label="scSNV-seq primer and pKLV2 table"):
    text=table.path.read_text(encoding="utf-8")
    rows={line.split("\t")[0]: line.split("\t")[1].strip()
          for line in text.splitlines() if line.count("\t")>=1}
    # Table S2 prints the PuroR arm in lower case; the model stores one case only.
    c("printed spike-in RT primers",
      (rows.get("oJR160 rt oligo","").upper(),rows.get("RT_Puro_Singlecell","").upper()),
      (P.GUIDE_RT,P.PURO_RT))
    c("both nested forward primers are one printed sequence",
      (rows.get("Kosuke_Guide_i5longFwd"),rows.get("Puro_i5longFwd")),(P.I5_LONG,P.I5_LONG))
    vector=rows.get("pKLV2 sequence added to hg38","").upper()
    for label,seq in (("scaffold loop 1",P.SCAFFOLD_LOOP1),
                      ("scaffold up to the oJR160 site",P.SCAFFOLD_AFTER_IBAR),
                      ("PuroR 5' UTR",P.PURO_UTR_TAIL),
                      ("PuroR start up to the RT-primer site",P.PURO_CDS_START)):
        c(f"{label} is transcribed verbatim from the published pKLV2 region",seq in vector,True)
    runs=[len(r) for r in re.findall(r"N+",vector)]
    c("published guide, iBAR and PuroR barcode windows",sorted(runs),[6,12,19])
c.report()
