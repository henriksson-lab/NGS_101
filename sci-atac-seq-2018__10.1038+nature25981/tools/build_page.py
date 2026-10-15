#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHARED = HERE.parents[1] / "sci-atac-seq__10.1126+science.aab1601" / "tools"
sys.path.insert(0, str(SHARED))
from page_parts import render_page
import sciatac
import chem

OUT = HERE.parent / "sci-atac-seq-2018.html"
if __name__ == "__main__":
    render_page(chem.PROTOCOL, OUT, source=chem.SOURCE,
                graph_factory=chem.protocol_graph,
                final_library=chem.FINAL_LIBRARY,
                pcr_n=chem.PCR_INDEX_NT, inferred=False,
                read_rows=chem.READ_ROWS, primers=chem.PRIMERS)

