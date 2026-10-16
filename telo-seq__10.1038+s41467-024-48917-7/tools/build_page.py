#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import teloseq as M
import telomere as T
from protocol_graph import ProtocolGraph
from telomere_page import render_page
OUT=HERE.parent/"telo-seq.html"
SOURCE='Schmidt et al. 2024, <a href="https://doi.org/10.1038/s41467-024-48917-7">doi:10.1038/s41467-024-48917-7</a>.'
TITLE="Telo-seq"

def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    current = T.start_telomere_graph(graph)
    tag, current = T.anneal_terminal_tag(graph, current,
        label="phased telomere-capture or TeloTag oligo")
    current = T.ligate_terminal_tag(graph, current, tag)
    current = T.cut_and_tail_distal_end(graph, current)
    current = T.recruit_nanopore_motor(graph, current)
    return T.attach_nanopore_motor(graph, current, title=TITLE, authority=SOURCE,
        reviewed=("genomic_dna", "terminal_tag", "anneal_tag", "terminal_ligase",
                  "ligate_terminal_tag", "ecorv", "digest_ecorv", "klenow",
                  "tail_distal", "s1_splint", "recruit_motor"))

def render():
    source=SOURCE
    title=TITLE
    return render_page(title=title,notes="01_telo-seq.html",citation_html=source,panels=M.panels(),entry="telomere-side motor-loaded ONT adapter ---> C-rich strand read telomere-outward; no sequencing primer",caveat="Source discrepancy: Schmidt et al. Methods (JATS line 169) specify six phased telorette adapters, but Supplementary Table 7 (PDF lines 1780–1810) labels both the third and fourth distinct sequences T3 and has no T4 row. The model retains all six printed sequences in table order, treating the fourth sequence as the fourth phase; the intended fourth label has not been confirmed by the authors. The cited ONT AMII motor-adapter terminal bases remain unresolved.",graph=protocol_graph())
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
