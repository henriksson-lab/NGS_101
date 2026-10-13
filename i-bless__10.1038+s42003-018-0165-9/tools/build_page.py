#!/usr/bin/env python3
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem
from batch_page import render
OUT = HERE.parent / "i-bless.html"
folds = "".join((
    '<h2>Hairpin topology</h2>',
    chem.PROXIMAL_TOPOLOGY.folded_panel(
        "The i-BLESS workflow retains the published BLESS L1 sequence: a 42-bp stem "
        "surrounding a seven-base loop, with the proximal tag at its molecular end."),
    chem.DISTAL_TOPOLOGY.folded_panel(
        "The distinct L3 adaptor has the same stem-loop topology and a different terminal "
        "tag, preserving which sequenced end arose during fragmentation."),
))
page = render(chem).replace('<h2>Reaction workflow</h2>',
                            folds + '<h2>Reaction workflow</h2>', 1)
OUT.write_text(page, encoding="utf-8")
print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
