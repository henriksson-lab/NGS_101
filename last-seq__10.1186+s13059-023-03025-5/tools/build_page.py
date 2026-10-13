#!/usr/bin/env python3
from pathlib import Path
import sys
H = Path(__file__).resolve().parent
sys.path[:0] = [str(H), str(H.parents[1] / "lib")]
import last_seq as P
from batch_page import render
OUT = H.parent / "last-seq.html"
if __name__ == "__main__":
    fold = '<h2>Capture-hairpin topology</h2>' + P.BC2_HAIRPIN_TOPOLOGY.folded_panel(
        "Published BC2 shown as one representative barcode member. The 11-base terminal "
        "arms form the stem; the PCR handle, T7 promoter, eight-base UMI and six-base "
        "cell barcode occupy the loop. The leading 17 bases remain single-stranded for "
        "annealing and nick ligation to the rU-dT module.")
    page = render(P).replace('<h2>Reaction workflow</h2>',
                             fold + '<h2>Reaction workflow</h2>', 1)
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}")
