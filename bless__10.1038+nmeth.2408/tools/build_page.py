#!/usr/bin/env python3
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem
from batch_page import render
OUT = HERE.parent / "bless.html"
folds = "".join((
    '<h2>Hairpin topology</h2>',
    chem.PROXIMAL_TOPOLOGY.folded_panel(
        "Published L1 proximal adaptor: a 42-bp stem surrounding the seven-base loop. "
        "The internal biotin-dT lies in that loop and is retained in the selectable adaptor."),
    chem.DISTAL_TOPOLOGY.folded_panel(
        "Published L3 distal adaptor: the same 42-bp stem and seven-base loop with its "
        "distinct terminal barcode. Both drawings retain the complete ordered sequence."),
))
page = render(chem).replace('<h2>Reaction workflow</h2>',
                            folds + '<h2>Reaction workflow</h2>', 1)
OUT.write_text(page, encoding="utf-8")
print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
