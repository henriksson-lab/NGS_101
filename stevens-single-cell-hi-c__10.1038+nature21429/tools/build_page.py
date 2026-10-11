#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import stevens as M
from single_cell_hic_page import render
from chemdraw import workflow_panel
OUT=HERE.parent/"stevens-single-cell-hi-c.html"
if __name__ == "__main__":
    page=render(M)
    extra='<h2>Adapter-end geometry</h2>'+workflow_panel(M.adapter_workflow(),cls="long")
    page=page.replace('<h2>Final sequencing molecule</h2>',extra+'<h2>Final sequencing molecule</h2>')
    OUT.write_text(page, encoding="utf-8"); print(f"wrote {OUT}")
