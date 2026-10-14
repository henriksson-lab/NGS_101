#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import dtm as M
from telomere import nanopore_telomere_graph
from telomere_page import render_page
OUT=HERE.parent/"digital-telomere-measurement.html"
def render():
    source='Sanchez et al. 2024, <a href="https://doi.org/10.1038/s41467-024-49007-4">doi:10.1038/s41467-024-49007-4</a>.'
    title="Digital Telomere Measurement (DTM)"
    return render_page(title=title,notes="01_digital-telomere-measurement.html",citation_html=source,panels=M.panels(),entry="motor-loaded ONT adapter ---> one intact captured strand through nanopore; no sequencing primer",caveat="The capture oligo and seqTether are source-exact. The defining study names historical AMII and NA adapters but does not print their complete motor-bearing assemblies; a later current-kit sequence is not substituted for them.",graph=nanopore_telomere_graph("dtm",title,source),oligos=M.oligos())
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
