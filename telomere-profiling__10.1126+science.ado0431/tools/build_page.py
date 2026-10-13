#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import telomere_profiling as M
from telomere import nanopore_telomere_graph
from telomere_page import render_page
OUT=HERE.parent/"telomere-profiling.html"
def render():
    source='Karimian et al. 2024, <a href="https://doi.org/10.1126/science.ado0431">doi:10.1126/science.ado0431</a>.'
    title="Telomere Profiling"
    return render_page(title=title,notes="01_telomere-profiling.html",citation_html=source,panels=M.panels(),entry="motor-loaded ONT adapter ---> one enriched telomere strand through nanopore; no sequencing primer",caveat="Adapter bases not verified from the supplementary protocol are shown as structural placeholders, never inferred sequence.",graph=nanopore_telomere_graph("profiling",title,source))
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
