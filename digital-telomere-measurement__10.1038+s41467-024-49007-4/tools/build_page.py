#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import dtm as M
from telomere_page import render_page
OUT=HERE.parent/"digital-telomere-measurement.html"
SOURCE='Sanchez et al. 2024, <a href="https://doi.org/10.1038/s41467-024-49007-4">doi:10.1038/s41467-024-49007-4</a>.'
TITLE="Digital Telomere Measurement (DTM)"


def render():
    source=SOURCE
    title=TITLE
    return render_page(title=title,notes="01_digital-telomere-measurement.html",citation_html=source,panels=M.panels(),entry="motor-loaded ONT adapter ---> one intact captured strand through nanopore; no sequencing primer",caveat="The capture oligo and seqTether are source-exact (Supplementary Data 3); the barcode/arm partition is derived from that table; the terminal register and drawn telomere lengths are schematic, and the sealed nick at the telomere end is not stated. The defining study names historical AMII and NA adapters but does not print their complete motor-bearing assemblies; a later current-kit sequence is not substituted for them. The cached author bench protocol is dated April 2025 and is a later revision; it does not supply the 2024 chemistry. Product connectivity and gap-fill tracts remain unresolved. Possible unreacted or partially adapted molecules are inferred co-resident classes; their abundance and bead-cleanup carryover are unassigned.",graph=M.protocol_graph(),oligos=M.oligos())
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
