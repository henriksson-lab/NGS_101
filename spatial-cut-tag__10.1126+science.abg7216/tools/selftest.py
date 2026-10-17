#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import spatial_cut_tag as P
import seqprimers as sp
from checks import Check, Source, have
from storage import data_root
c=Check(); c.section("Spatial-CUT&Tag source transcription")
# The 2022 Tn5ME-A appears only in the Science Supplementary Materials (DOCX text
# line 912, "Transposome preparation"); the 2024 table prints a different sequence.
supp = Source(data_root() / "sources" / H.parent.name /
              "science.abg7216_supp_EMS145971-supplement-Supplementary_Material.docx.txt",
              "Deng et al. 2022 Supplementary Materials DOCX (PMC bin, manual download)")
if have(c, supp, label="2022 Supplementary Materials"):
    text = " ".join(supp.path.read_text(errors="replace").split())
    c("2022 Tn5ME-A as printed", f"5′-{P.TN5_ME_A_WRITTEN}-3′" in text, True)
c("representative first barcode oligo", P.BARCODE_A_OLIGO_1,
  "/5Phos/AGGCCAGAGCATTCGAACGTGATAGATGTGTATAAGAGACAG")
c("2024 Supp_Table 1 Tn5ME-A cannot pair with ligation linker 1",
  P.TABLE1_TN5_ME_A_PAIRS_LINKER_1, False)
c("authors' raw Read 1 barcode coordinates",
  [(x.segment, x.cycles) for x in sp.feature_spans(
      P.FINAL_LIBRARY, P.SEQ_PRIMERS[:1], {"Read 1": 150})
   if x.feature.role == "spatial_barcode"],
  [("spatial barcode B", "1–8"), ("spatial barcode A", "39–46")])
c.report()
