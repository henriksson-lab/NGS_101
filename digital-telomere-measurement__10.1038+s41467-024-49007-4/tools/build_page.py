#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import dtm as M
from telomere_page import render_page
OUT=HERE.parent/"digital-telomere-measurement.html"
def render(): return render_page(title="Digital Telomere Measurement (DTM)",notes="01_digital-telomere-measurement.html",citation_html='Sanchez et al. 2024, <a href="https://doi.org/10.1038/s41467-024-49007-4">doi:10.1038/s41467-024-49007-4</a>.',panels=M.panels(),entry="motor-loaded ONT adapter ---> one intact captured strand through nanopore; no sequencing primer",caveat="The article establishes the capture/tether structure and reaction order; the ordered oligo bases are left as explicit unknowns until its supplementary table is locally verified.")
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
