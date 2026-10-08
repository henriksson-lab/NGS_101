#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import telobait as M
from telomere_page import render_page
OUT=HERE.parent/"pacbio-telobait.html"
def render(): return render_page(title="PacBio telobait HiFi sequencing",notes="01_pacbio-telobait.html",citation_html='Tham et al. 2023, <a href="https://doi.org/10.1038/s41467-023-35823-7">doi:10.1038/s41467-023-35823-7</a>.',panels=M.panels(),entry="PacBio sequencing primer + polymerase ---> hairpin primer site ---> repeated passes around the telobait SMRTbell")
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
