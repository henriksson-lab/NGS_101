#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import dtm as M
import telomere as T
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from telomere_page import render_page
OUT=HERE.parent/"digital-telomere-measurement.html"
SOURCE='Sanchez et al. 2024, <a href="https://doi.org/10.1038/s41467-024-49007-4">doi:10.1038/s41467-024-49007-4</a>.'
TITLE="Digital Telomere Measurement (DTM)"

def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    current = T.start_telomere_graph(graph)
    tag, current = T.anneal_terminal_tag(graph, current,
        label="published barcoded capture/seqTether duplex")
    # Methods (PMC11189511 Telomere capture library preparation): T4 DNA ligase in
    # rCutSmart + ATP overnight at 37 C, then 65 C 10 min heat inactivation.
    ligase = graph.input("terminal_ligase", "T4 DNA ligase (overnight, 37 °C; 65 °C heat-inactivated)",
                         MaterialKind.CATALYST)
    tagged = graph.product("terminal_tagged", "capture-oligo-ligated telomere DNA",
                           MaterialKind.MOLECULAR, state=current.state)
    graph.react("ligate_terminal_tag", "Ligate the capture duplex to native telomere ends",
                StepKind.LIGATION,
                ((current, InputRole.TEMPLATE), (tag, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (tagged,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The Methods and Supplementary Data 3 give the capture-oligo 5-prime phosphate and its CCTAACC telomere arm, but do not state which nick at the annealed telomere end T4 ligase seals."))
    # Methods: same tube, 4 U Sulfolobus DNA polymerase IV in ThermoPol + dNTPs/ATP,
    # 56 C 2 min then 72 C 15 min.
    current = T.fill_terminal_gap(graph, tagged,
        polymerase="Sulfolobus DNA polymerase IV (4 U; 56 °C 2 min, 72 °C 15 min)",
        source="Sanchez et al. 2024 Methods, Telomere capture library preparation")
    beads = graph.input("pronex_beads", "Promega ProNex size-selection beads (1.6×)",
                        MaterialKind.SELECTION)
    library = graph.product("post_capture", "post-capture sequencing library (ONT elution buffer)",
                            MaterialKind.MOLECULAR, state=current.state)
    waste = graph.product("pronex_waste", "bead supernatant and ethanol washes",
                          MaterialKind.WASTE)
    graph.react("pronex_cleanup", "ProNex bead cleanup of the gap-filled ligation",
                StepKind.CLEANUP,
                ((current, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (library, waste))
    graph.discard(waste, "removed with the supernatant and two 80% ethanol washes")
    return T.attach_nanopore_motor(graph, library, title=TITLE, authority=SOURCE,
        reviewed=("genomic_dna", "terminal_tag", "anneal_tag", "terminal_ligase",
                  "ligate_terminal_tag", "gap_fill_polymerase", "gap_fill",
                  "pronex_beads", "pronex_cleanup", "post_capture"))

def render():
    source=SOURCE
    title=TITLE
    return render_page(title=title,notes="01_digital-telomere-measurement.html",citation_html=source,panels=M.panels(),entry="motor-loaded ONT adapter ---> one intact captured strand through nanopore; no sequencing primer",caveat="The capture oligo and seqTether are source-exact (Supplementary Data 3); the barcode/arm partition is derived from that table, and the sealed nick at the telomere end is not stated. The defining study names historical AMII and NA adapters but does not print their complete motor-bearing assemblies; a later current-kit sequence is not substituted for them. Its linked author bench protocol (TelometerLibraryPrep.docx) remains unread because the download failed.",graph=protocol_graph(),oligos=M.oligos())
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
