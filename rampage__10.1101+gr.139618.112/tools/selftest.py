#!/usr/bin/env python3
import re
import build_page
import storage
from checks import Check, Source, have

M = build_page.M
build_page.render(M)
print("PASS  RAMPAGE construction and sequencing-primer geometry")

check = Check()
check.section("Table 25B.11.1 TSO barcode transcription")
check("40 barcodes, all distinct six-mers",
      (len(M.TSO_BARCODES), len(set(M.TSO_BARCODES)),
       {len(b) for b in M.TSO_BARCODES}), (40, 40, {6}))
# Current Protocols 'Sequence barcodes': GC 20-80% and pairwise Hamming distance >= 3.
check("transcribed set meets the stated minimum Hamming distance 3",
      M.min_pairwise_hamming() >= 3)
check("transcribed set meets the stated GC range 20-80%",
      all(2 <= sum(c in "GC" for c in b) <= 4 for b in M.TSO_BARCODES))
thesis = Source(storage.sources_root() / "rampage__10.1101+gr.139618.112" / "Batut_Thesis.pdf.txt",
                "https://repository.cshl.edu/id/eprint/31200/1/Batut%20Thesis.pdf "
                "(python3 tools/get_sources.py --slugs rampage__10.1101+gr.139618.112)")
if have(check, thesis):
    printed = re.findall(M.TSO_HANDLE + r"([ACGT]{6})rGrGrG", thesis.path.read_text())
    check("thesis reproduction of Table 25B.11.1 prints the same 40 barcodes as a set",
          sorted(set(printed)), sorted(M.TSO_BARCODES))
check.report()
