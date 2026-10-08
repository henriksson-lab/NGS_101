#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import teloseq as M
from telomere_page import render_page
OUT=HERE.parent/"telo-seq.html"
def render(): return render_page(title="Telo-seq",notes="01_telo-seq.html",citation_html='Schmidt et al. 2024, <a href="https://doi.org/10.1038/s41467-024-48917-7">doi:10.1038/s41467-024-48917-7</a>.',panels=M.panels(),entry="telomere-side motor-loaded ONT adapter ---> C-rich strand read telomere-outward; no sequencing primer")
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
