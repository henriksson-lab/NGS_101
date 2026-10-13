#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import telotag as M
from telomere import nanopore_telomere_graph
from telomere_page import render_page
OUT=HERE.parent/"yeast-telotag.html"
def render():
    source='Sholes et al. 2022, <a href="https://doi.org/10.1101/gr.275868.121">doi:10.1101/gr.275868.121</a>.'
    title="Yeast TeloTag nanopore sequencing"
    return render_page(title=title,notes="01_yeast-telotag.html",citation_html=source,panels=M.panels(),entry="motor-loaded ONT adapter ---> one intact tagged DNA strand through nanopore; no sequencing primer",graph=nanopore_telomere_graph("telotag",title,source))
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
