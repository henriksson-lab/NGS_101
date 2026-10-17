#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import teloseq as M
from telomere_page import render_page
OUT=HERE.parent/"telo-seq.html"
protocol_graph=M.protocol_graph

def render():
    return render_page(title=M.TITLE,notes="01_telo-seq.html",citation_html=M.SOURCE,panels=M.panels(),
        entry="telomere-side motor-loaded ONT adapter ---> C-rich strand read telomere-outward; no sequencing primer",
        caveat=M.DISCREPANCY+" The AMII motor adapter and its ligation reagents follow the ONT SQK-NBD111 kit instructions, which the article does not reproduce; their terminal bases remain unresolved, as does whether the 6-nt 3′ extension of S1 (beyond the telorette-core complement) is the AMII-pairing region.",
        graph=protocol_graph(),oligos=M.oligos())
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
