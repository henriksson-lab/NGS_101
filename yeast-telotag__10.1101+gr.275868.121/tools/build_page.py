#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import telotag as M
import telomere as T
from protocol_graph import ProtocolGraph
from telomere_page import render_page
OUT=HERE.parent/"yeast-telotag.html"
SOURCE='Sholes et al. 2022, <a href="https://doi.org/10.1101/gr.275868.121">doi:10.1101/gr.275868.121</a>.'
TITLE="Yeast TeloTag nanopore sequencing"

def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    current = T.start_telomere_graph(graph)
    current = T.add_poly_a_tail(graph, current)
    tag, current = T.anneal_terminal_tag(graph, current,
        label="phased telomere-capture or TeloTag oligo")
    current = T.fill_terminal_gap(graph, current)
    current = T.seal_terminal_nick(graph, current)
    return T.attach_nanopore_motor(graph, current, title=TITLE, authority=SOURCE,
        reviewed=("genomic_dna", "tdt", "polya_tail", "terminal_tag", "anneal_tag",
                  "gap_fill_polymerase", "gap_fill", "nick_ligase", "seal_nick"))

def render():
    source=SOURCE
    title=TITLE
    return render_page(title=title,notes="01_yeast-telotag.html",citation_html=source,panels=M.panels(),entry="motor-loaded ONT adapter ---> one intact tagged DNA strand through nanopore; no sequencing primer",caveat="The five TeloTag primers are printed in Supplemental Table S5. The length of the enzymatically added poly(A) tail varies; the ONT kit does not disclose its historical motor-adapter oligos.",graph=protocol_graph(),oligos=M.oligos())
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
