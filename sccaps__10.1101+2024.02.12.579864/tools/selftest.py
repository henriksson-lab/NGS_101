#!/usr/bin/env python3
import re
import build_page
import storage
from checks import Check, Source, have

M = build_page.M
c = Check(); c("published CAPS+ conversion paths",
             tuple(M.CAPS_PATHS[x].states for x in ("5hmC", "5mC", "C")),
             (("5hmC", "5fC", "5caC", "DHU", "T"), ("5mC", "5mC"), ("C", "C")))

# Source transcription of the imported joint-paper tables, against this bundle's copy.
table = Source(storage.sources_root() / "sccaps__10.1101+2024.02.12.579864" /
               "s13059-025-03708-1_supp_13059_2025_3708_MOESM5_ESM.xlsx.txt",
               "Genome Biology 2025 Additional file 5 (PMC12359873 MOESM5)")
if have(c, table, label="Additional file 5 oligo tables"):
    text = table.path.read_text(encoding="utf-8")
    printed = dict(re.findall(r"(P[57]_i[57]_\d+)_Universal\S*, /5Phos/([ACGT]+)", text))
    printed["pMENTS"] = re.search(r"pMENTS,\s*/5Phos/([ACGT]+)", text).group(1)
    c("Table S3 transfer strands and pMENTS", M.shared.table_s3_oligos(), printed)
    s4 = dict(re.findall(r"(Read 1|Read 2|Index 1|Index 2) [Ss]equencing [Pp]rimer:\s*([ACGT]+)", text))
    c("Table S4 custom sequencing primers", M.shared.TABLE_S4, s4)

build_page.render(M); c.report()
