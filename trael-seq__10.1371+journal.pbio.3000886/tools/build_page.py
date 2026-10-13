#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem
from batch_page import render
OUT=H.parent/"trael-seq.html"
if __name__=="__main__":
    fold = '<h2>Hairpin topology</h2>' + chem.ADAPTER2_TOPOLOGY.folded_panel(
        "The fully published second adaptor folded base by base: a 31-bp stem, a "
        "four-dU loop and the terminal phosphorothioate-protected 3′ T outside the stem. "
        "Adapter 1 is not folded here because its complete intervening sequence is not disclosed.")
    page = render(chem).replace('<h2>Reaction workflow</h2>',
                                fold + '<h2>Reaction workflow</h2>', 1)
    OUT.write_text(page,encoding="utf-8"); print(f"wrote {OUT}")
