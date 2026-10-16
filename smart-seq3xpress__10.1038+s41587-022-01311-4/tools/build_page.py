#!/usr/bin/env python3
"""Build the SMART-seq3xpress protocol schematic."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHARED = HERE.parents[1] / "smart-seq__10.1038+nbt.2282" / "tools"
sys.path.insert(0, str(SHARED))

import page_parts as P  # noqa: E402
from protocol_graph import ProtocolGraph  # noqa: E402

OUT = HERE.parent / "smart-seq3xpress.html"

# The original xpress Methods specify 10-base dual indexes. The protocol's
# attached 384set_Nextera_BC_order.xls supplies complete oligos in two sheets.
# A1 on each sheet is an illustrative pair, not a recovered experimental well.
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Segment, Scene, feature, oligo, panel, protocol_graph_panels, revcomp
from page import head, info
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, SequencingPlan, StepKind, AnnotationSlot
from protocol_views import assemble_page

I5_A1 = "TTGGTACGCG"
I7_A1 = "ATAGGCGCTC"
INDEX_FORWARD = il.P5 + I5_A1 + nx.S5
INDEX_REVERSE = il.P7 + I7_A1 + nx.S7
SEQ_PRIMERS = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2", "R2"))
TAG = "ATTGCGCAATG"
UMI = feature("umi", "umi", "random")
I5 = feature("sample_index_i5", "sample_index", "fixed")
I7 = feature("sample_index_i7", "sample_index", "fixed")


def final_library(*, five_prime=True):
    middle = ([Segment("5-prime tag", TAG, "r1"),
               Segment("UMI", "N" * 8, "umi", placeholder=True, feature=UMI),
               Segment("WW spacer", "WW", "r1", placeholder=True),
               Segment("template-switch GGG", "GGG", "tso")]
              if five_prime else [])
    library = Construct([
        Segment("P5", il.P5, "p5"),
        Segment("illustrative i5 A1", I5_A1, "cbc", feature=I5),
        Segment("s5", nx.S5, "s5"), Segment("ME", nx.ME, "me"),
        *middle, Segment("cDNA", "X" * 32, placeholder=True),
        Segment("ME reverse complement", nx.ME_RC, "me"),
        Segment("s7 reverse complement", nx.S7_RC, "s7"),
        Segment("illustrative i7 A1 reverse complement", revcomp(I7_A1), "cbc", feature=I7),
        Segment("P7 reverse complement", il.P7_RC, "p7")],
        name="SMART-seq3xpress " + ("5-prime" if five_prime else "internal") + " library; illustrative A1/A1")
    require_index_sites(library)
    problems = sp.verify(library, SEQ_PRIMERS,
                         required_roles=tuple(p.role for p in SEQ_PRIMERS))
    if problems:
        raise ValueError("; ".join(problems))
    return library


def require_index_sites(library):
    return PCRGeometry.exact(library,
        PCRPrimer("UD_i5_indices A1", INDEX_FORWARD, anneal_nt=len(INDEX_FORWARD)),
        PCRPrimer("UD_i7_indices A1", INDEX_REVERSE, anneal_nt=len(INDEX_REVERSE)),
        statement="Complete original workbook A1 primers match the displayed illustrative 10-base dual-index endpoint.")


def index_pcr(graph, filled):
    forward = graph.input("index_forward", "Custom i5 A1 primer (illustrative)", MaterialKind.OLIGO,
                          state=Construct([Segment("ordered i5 A1", INDEX_FORWARD)]))
    reverse = graph.input("index_reverse", "Custom i7 A1 primer (illustrative)", MaterialKind.OLIGO,
                          state=Construct([Segment("ordered i7 A1", INDEX_REVERSE)]))
    polymerase = graph.input("index_polymerase", "Phusion DNA polymerase", MaterialKind.CATALYST)
    five, internal = final_library(), final_library(five_prime=False)
    library = graph.product("library", "5-prime library, illustrative A1/A1", MaterialKind.LIBRARY, state=five)
    other = graph.product("internal_library", "Internal library, illustrative A1/A1", MaterialKind.LIBRARY, state=internal)
    graph.react("index_pcr", "Index PCR: source 10-base dual indexes; illustrative A1/A1", StepKind.PCR,
                ((filled, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library, other), geometry=require_index_sites(five))
    for material, construct in ((library, five), (other, internal)):
        graph.sequence(material, SequencingPlan(library=construct, primers=SEQ_PRIMERS,
                       acquired_roles=tuple(p.role for p in SEQ_PRIMERS)))


def render_page():
    graph = protocol_graph()
    # The RT/TSO panels remain the previously checked SMART-seq3xpress declarations;
    # endpoint, index inputs, final panels and sequencing are owned here.
    ordered = P.ss3_oligos("SMART-seq3xpress").split("</seq>", 1)[0] + "</seq>"
    ordered += info("The two terminal phosphorothioate linkages are marked *. The original xpress Methods use 10-base dual indexes. Both A1 oligos below come from the author-linked index workbook; their pairing is illustrative and does not assign a barcode to a published cell.")
    ordered += "<seq>" + oligo("Illustrative i5 A1", [Segment("", INDEX_FORWARD)]) + oligo("Illustrative i7 A1", [Segment("", INDEX_REVERSE)]) + "</seq>"
    finals, sequencing = [], []
    for five_prime in (True, False):
        con = final_library(five_prime=five_prime)
        heading = "5-prime fragments" if five_prime else "Internal fragments"
        finals.extend(["<h3>" + heading + "</h3>", panel(Scene.duplex(con.segments).rows(), cls="long")])
        sequencing.append(sp.section(con, SEQ_PRIMERS, heading="Sequencing: " + heading,
                                     show_read_content=True))
    return assemble_page(graph, [head("SMART-seq3xpress library chemistry"),
        P.preamble("SMART-seq3xpress", complete=True),
        info("This page represents the Illumina-compatible dual-index library before any MGI conversion. The defining paper also reports MGI adapter conversion/circularization and custom sequencing primers; that branch is outside this schematic."),
        AnnotationSlot.INTRO, ordered,
        "<h2>Step-by-step library generation</h2>", protocol_graph_panels(graph, cls="long"),
        AnnotationSlot.AFTER_WORKFLOW, "<h2>Final libraries</h2>", *finals,
        "<h2>Library sequencing</h2>", *sequencing,
        info("Index source: author protocol DOI yxmvmk1yng3p/v1, attached 384set_Nextera_BC_order.xls, sheets UD_i5_indices and UD_i7_indices, A1. The workbook defines available oligos; it does not assign this illustrative pair to an experimental well."),
        P.sources("SMART-seq3xpress"), "</div>"])


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: doi:10.1038/s41587-022-01311-4"
    graph = ProtocolGraph("SMART-seq3xpress", authority=source)
    mrna = P.input_mrna(graph)
    first, enzyme, switch_scene, template = P.reverse_transcribe(
        graph, "SMART-seq3xpress", mrna, rt_label="biotinylated anchored oligo(dT) primer")
    switched = P.template_switch(
        graph, first, enzyme, switch_scene, template,
        tso_label="SMART-seq3xpress UMI template-switch oligo",
        step_label="Template-switch at the transcript 5-prime end")
    amplified = P.preamplify(
        graph, switched, template, family="ss3",
        statement="The published asymmetric forward and reverse primers face inward around the full-length template-switched cDNA.")
    tagged = P.tagment(graph, amplified)
    filled = P.fill_tn5_gaps(graph, tagged)
    index_pcr(graph, filled)
    return P.review_smart_graph(graph, source, family="ss3", internal=True)


def main() -> None:
    OUT.write_text(render_page(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
