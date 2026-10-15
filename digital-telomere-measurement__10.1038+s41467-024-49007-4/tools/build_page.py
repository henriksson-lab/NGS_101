#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import dtm as M
import telomere as T
from protocol_graph import ProtocolGraph
from telomere_page import render_page
OUT=HERE.parent/"digital-telomere-measurement.html"
SOURCE='Sanchez et al. 2024, <a href="https://doi.org/10.1038/s41467-024-49007-4">doi:10.1038/s41467-024-49007-4</a>.'
TITLE="Digital Telomere Measurement (DTM)"

def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    current = T.start_telomere_graph(graph)
    tag, current = T.anneal_terminal_tag(graph, current,
        label="published barcoded capture/seqTether duplex")
    current = T.ligate_terminal_tag(graph, current, tag)
    current = T.fill_terminal_gap(graph, current)
    return T.attach_nanopore_motor(graph, current, title=TITLE, authority=SOURCE,
        reviewed=("genomic_dna", "terminal_tag", "anneal_tag", "terminal_ligase",
                  "ligate_terminal_tag", "gap_fill_polymerase", "gap_fill"))

def render():
    source=SOURCE
    title=TITLE
    return render_page(title=title,notes="01_digital-telomere-measurement.html",citation_html=source,panels=M.panels(),entry="motor-loaded ONT adapter ---> one intact captured strand through nanopore; no sequencing primer",caveat="The capture oligo and seqTether are source-exact. The defining study names historical AMII and NA adapters but does not print their complete motor-bearing assemblies; a later current-kit sequence is not substituted for them. Its linked author bench protocol (TelometerLibraryPrep.docx) remains unread because the download failed.",graph=protocol_graph(),oligos=M.oligos())
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
