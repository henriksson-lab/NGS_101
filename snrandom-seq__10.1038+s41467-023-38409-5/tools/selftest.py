#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import snrandom_seq as M
from checks import Check
c=Check(); c.section("snRandom-seq source transcription")
c("block primer",M.BLOCK_PRIMER,"GAGAATGTGAGTGAAGATGTATGGTGANNNNNNN")
c("cell barcode length",sum(len(M.FINAL_LIBRARY.get(f"cell barcode {i}")) for i in (1,2,3)),30)
c("UMI length",len(M.FINAL_LIBRARY.get("UMI")),8)
c("acrydite bead anchor",M.BEAD_ANCHOR_WRITTEN,"/5Acryd/ATTATATATATUGTGAGTGATGGTTGAGGATGTGTGGAGATA")
c("bead oligo published T tract",M.barcode_parts()[-1].top,"TTTTTTTTTTTTTTTTTTTT")
c("PCR primer1",M.PCR1,"GGAGTTGGAGTGAGTGGATGAGTGATG")
c("PCR primer2",M.PCR2,"GTGAGTGATGGTTGAGGATGTGTGGAGATA")
c.report()
