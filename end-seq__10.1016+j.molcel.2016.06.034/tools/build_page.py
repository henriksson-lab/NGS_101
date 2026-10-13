#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import end_seq as M
from batch_page import render
OUT=H.parent/"end-seq.html"
if __name__ == "__main__":
    folds = "".join((
        '<h2>USER-hairpin topology</h2>',
        M.A1_TOPOLOGY.folded_panel(
            "Adaptor 1 base by base: a 31-bp stem, an eight-base U-bearing loop and "
            "one terminal 3′ T ligation overhang. The two biotin-dT positions lie "
            "in the loop and appear as their underlying T bases."),
        M.A2_TOPOLOGY.folded_panel(
            "Adaptor 2 base by base: a 19-bp stem, a nine-base poly(dU) loop and one "
            "terminal 3′ T ligation overhang."),
    ))
    page = render(M).replace('<h2>Reaction workflow</h2>',
                             folds + '<h2>Reaction workflow</h2>', 1)
    OUT.write_text(page,encoding="utf-8"); print(f"wrote {OUT}")
