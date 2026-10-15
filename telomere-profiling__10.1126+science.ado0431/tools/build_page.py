#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import telomere_profiling as M
import telomere as T
from protocol_graph import ProtocolGraph
from telomere_page import render_page
OUT=HERE.parent/"telomere-profiling.html"
SOURCE='Karimian et al. 2024, <a href="https://doi.org/10.1126/science.ado0431">doi:10.1126/science.ado0431</a>.'
TITLE="Telomere Profiling"

def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    current = T.start_telomere_graph(graph)
    current = T.digest_internal_telomere(graph, current)
    tag, current = T.anneal_terminal_tag(graph, current,
        label="published phased TeloTag duplex")
    current = T.ligate_terminal_tag(graph, current, tag)
    current = T.capture_and_release_tag(graph, current)
    return T.attach_nanopore_motor(graph, current, title=TITLE, authority=SOURCE,
        reviewed=("genomic_dna", "internal_nuclease", "digest_internal", "terminal_tag",
                  "anneal_tag", "terminal_ligase", "ligate_terminal_tag",
                  "streptavidin", "capture", "release_nuclease", "release"))

def render():
    source=SOURCE
    title=TITLE
    return render_page(title=title,notes="01_telomere-profiling.html",citation_html=source,panels=M.panels(),entry="motor-loaded ONT adapter ---> one enriched telomere strand through nanopore; no sequencing primer",caveat="The six barcoded TeloTags and their six phase-splint families are published and modeled exactly. Only the complete historical motor-bearing ONT assembly remains unavailable.",graph=protocol_graph(),oligos=M.oligos())
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
