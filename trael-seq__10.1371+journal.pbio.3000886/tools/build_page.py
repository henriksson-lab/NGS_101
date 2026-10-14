#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem
from batch_page import render
OUT=H.parent/"trael-seq.html"
if __name__=="__main__":
    fold = '<h2>Hairpin topology</h2>' + chem.ADAPTER1_TOPOLOGY.folded_panel(
        "The defining paper prints adapter 1 base for base. Its 23-bp stem surrounds "
        "the uracil/biotin loop; the 8-nt UMI and ligating 5′ end remain outside the stem.")
    fold += chem.ADAPTER2_TOPOLOGY.folded_panel(
        "The fully published second adaptor folded base by base: a 31-bp stem, a "
        "four-dU loop and the terminal phosphorothioate-protected 3′ T outside the stem.")
    page = render(chem).replace('<h2>Reaction workflow</h2>',
                                fold + '<h2>Reaction workflow</h2>', 1)
    OUT.write_text(page,encoding="utf-8"); print(f"wrote {OUT}")
