#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import telotag as M
from telomere_page import render_page
OUT=HERE.parent/"yeast-telotag.html"
SOURCE=M.SOURCE
TITLE=M.TITLE
protocol_graph=M.protocol_graph

def render():
    source=SOURCE
    title=TITLE
    return render_page(title=title,notes="01_yeast-telotag.html",citation_html=source,panels=M.panels(),entry="motor-loaded ONT adapter ---> one intact tagged DNA strand through nanopore; no sequencing primer",caveat="Source discrepancies: the final Methods (Sholes et al. 2022, TeloTag addition, JATS line 189) point to Supplemental Table S6 for the five TeloTag primers, but S6 lists plasmids (supplement PDF lines 711–729); Supplemental Table S5 (primer rows, lines 679–709) prints the five sequences. The model uses S5 because it contains the named primers; the table-number error is unexplained. The length of the enzymatically added poly(A) tail varies; the article does not print the EXP-NBD104 barcodes, kit end-preparation reagents or SQK-LSK109 motor-adapter oligos, so native barcoding is a warned kit boundary. TdT can tail any accessible 3′ end; telomere and TeloTag reads are selected computationally after sequencing, not separated physically.",graph=protocol_graph(),oligos=M.oligos())
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
