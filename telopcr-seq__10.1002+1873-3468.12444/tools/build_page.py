#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import telopcr as M
from telomere_page import render_page
OUT=HERE.parent/"telopcr-seq.html"
def render():
    graph=M.protocol_graph()
    return render_page(title="TeloPCR-seq", notes="01_telopcr-seq.html", citation_html='Bennett et al. 2016, <a href="https://doi.org/10.1002/1873-3468.12444">doi:10.1002/1873-3468.12444</a>.', panels=M.panels(), entry="PacBio sequencing primer + polymerase ---> hairpin primer site ---> repeated passes around the SMRTbell", caveat="The defining Methods print the terminal anchor and both subtelomeric PCR primers. The historical PacBio hairpin and sequencing-primer site sequences are not disclosed. Library construction followed the manufacturer&#8217;s ~250 bp protocol, whose kit revision, end processing and hairpin ligating ends the article does not state; that step is drawn as an unresolved boundary. Values come from the published FEBS Letters article; the superseded bioRxiv preprint is not used.", graph=graph, oligos=M.oligos())
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
