#!/usr/bin/env python3
from pathlib import Path
import re
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
from checks import Check, Source, have
import storage
import picoplex_dna_seq as model

c=Check(); c("finished graph is publishable", model.protocol_graph().publishable, True)

# Source transcription: the Figure 7 index table in the manual's text twin lists the
# names in order and then the sequences in the same order (i5 block, then i7 block).
MANUAL = Source(storage.sources_root() / "picoplex-dna-seq" /
                "PicoPLEX-DNA-seq-Manual-QAM-106-002-updated_RB3381_.pdf.txt",
                "python3 tools/get_sources.py --slugs picoplex-dna-seq (catalogue/extra_sources.tsv)")
if have(c, MANUAL, label="index table transcription"):
    text = MANUAL.path.read_text(encoding="utf-8")
    for table in (model.I5_INDEXES, model.I7_INDEXES):
        names = r"\s+".join(table)
        seqs = r"\s+".join(table.values())
        c(f"manual prints {next(iter(table))}… names then sequences in order",
          bool(re.search(names + r"\s+" + seqs, text)), True)
c.report()
