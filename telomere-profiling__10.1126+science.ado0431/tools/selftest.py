#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import telomere_profiling as M
from checks import Check,run_common
c=Check(); c.section("Telomere Profiling source transcription")
c("NB50 complete TeloTag",M.TELOTAG_TOPS["NB50"],"CCCTCCGATAATGGACTTTGGTAACTTCCTGCGTTGCGTACAGCAATCAGGCGCGCCGAATTCCGATCGTTAATTAATGGTGGACTCTATGACCGTTCA")
c("NB50.C1 complete splint",M.splint_family("NB50")[0],"GAGTCCACCATTAATTAACGATCGGAATTCGGCGCGCCTGATTGCTGTACGCAACGCAGGAAGTTACCAAAGTCCATTATCGGAGGGACCCTAACCCTAACCCTA")
c("capture has release step",len(M.CAPTURE.steps()),3)
import re, storage
from checks import Source, have
SM=Source(storage.sources_root()/"telomere-profiling__10.1126+science.ado0431"/"science.ado0431_sm.pdf.txt",
          "Science Supplementary Materials https://www.science.org/doi/suppl/10.1126/science.ado0431/suppl_file/science.ado0431_sm.pdf (run tools/doctext.py)")
if have(c, SM, label="published Table S3 / fig. S1A transcription"):
    text=SM.path.read_text(encoding="utf-8")
    table=dict(re.findall(r"^(NB\d+)\s*\n\s*/5Phos/([ACGT]+)/3BioTEG/", text, re.M))
    c("Table S3 TeloTags equal the model (14 rows)", table, M.TELOTAG_TOPS)
    rows=[l[3:-3][::-1] for l in text.splitlines() if l.startswith("3’-") and l.endswith("-5’")]
    c("fig. S1A NB50 splints (written 3′→5′) equal the derived family", sorted(rows), sorted(M.splint_family("NB50")))
run_common(c); c.report()
