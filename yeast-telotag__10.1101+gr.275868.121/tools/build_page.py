#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import telotag as M
from telomere_page import render_page
OUT=HERE.parent/"yeast-telotag.html"
def render(): return render_page(title="Yeast TeloTag nanopore sequencing",notes="01_yeast-telotag.html",citation_html='Sholes et al. 2022, <a href="https://doi.org/10.1101/gr.275868.121">doi:10.1101/gr.275868.121</a>.',panels=M.panels(),entry="motor-loaded ONT adapter ---> one intact tagged DNA strand through nanopore; no sequencing primer")
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
