#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import rnascope as M
from branched_dna_page import render
from page import info
OUT=HERE.parent/"rnascope.html"


def render_page():
    page = render(M)
    heading = '<h2>(3) Assemble the branched-DNA signal tree</h2>'
    fold = ('<h3>Nucleotide-resolved double-Z fold</h3>'
            + M.SEQUENCE_TOPOLOGY.folded_panel(
                "The patent example with every disclosed base shown. Target RNA is "
                "kept straight; five-thymidine spacers are flexible joints, and the "
                "two 14-nt Z tails jointly form the preamplifier docking duplex."))
    if heading not in page:
        raise ValueError("branched-DNA page no longer exposes the topology insertion point")
    return page.replace(heading, fold + heading + "\n" + info(M.HYBRIDIZATION_SERIES), 1)


if __name__=="__main__": OUT.write_text(render_page(),encoding="utf-8"); print(f"wrote {OUT}")
