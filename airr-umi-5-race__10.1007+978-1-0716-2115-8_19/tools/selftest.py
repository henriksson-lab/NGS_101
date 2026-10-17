#!/usr/bin/env python3
import build_page
from checks import Check, Source, have
from storage import data_root

M = build_page.M
c = Check()
chapter = Source(data_root() / "sources" / build_page.H.parents[0].name / "NBK586949_ch19.pdf.txt",
                 "NCBI Bookshelf OA package spr9781071621158_NBK586966.tar.gz, ch19.pdf (DOI 10.1007/978-1-0716-2115-8_19)")
if have(c, chapter, label="published AIRR chapter text"):
    text = " ".join(chapter.path.read_text().split())
    c("chapter i7 index list transcribed (section 2.2.2)",
      all(f"Universal Forward {n}, {d}, {s}." in text for n, (d, s) in enumerate(M.I7_INDEXES, 1)), True)
    c("chapter i5 index list transcribed (section 2.2.3)",
      all(f"IgG Reverse {n}, {d}, {s}." in text for n, (d, s) in enumerate(M.I5_INDEXES, 1)), True)
    c("chapter UMI start and linker/template-switch length",
      f"--start {M.UMI_LENGTH} --len {M.SMART_LINKER_LENGTH}" in text, True)
build_page.render(M); c.report()
