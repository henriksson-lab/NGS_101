"""Cusanovich 2018 two-round sci-ATAC chemistry, graph, and page."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] /
                       "sci-atac-seq__10.1126+science.aab1601" / "tools"))

import sciatac as S
from sciatac import indexed_transposome, pool_indexed_nuclei, seg, two_tagmented
from chemdraw import Construct
import illumina as il
import nextera as nx
from page import head, info
from end_combinations import (EndCombinationGeometry, EndProductGroup,
                              LoadedEnd, LoadedEndPool)
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (AnnotationSlot, DetailView, Highlight, InputRole, MaterialKind,
                            MolecularRelation, ProtocolGraph, SequencingPlan, StepKind)
from protocol_views import (OligoEntry, OligoSection, ReadLayoutSection, SourceRead,
                            WorkflowSection, assemble_page)

PROTOCOL = S.SCI18
PAGE_FILE = "sci-atac-seq-2018.html"
SOURCE = ("10.1038/nature25981", PROTOCOL)
PCR_INDEX_NT = 10
SOURCE_P5_1 = "AATGATACGGCGACCACCGAGATCTACACCTCCATCGAGTCGTCGGCAGCGTC"
FINAL_LIBRARY = S.two_level_library(PROTOCOL, pcr_n=PCR_INDEX_NT, inferred=False)
PRIMERS = S.TWO_LEVEL_PRIMERS
READ_ROWS = (
    SourceRead("Read 1", "50", "genomic insert", "Read 1", "include"),
    SourceRead("Index 1", "8 + 27 dark + 10", "t7 barcode + PCR i7",
               "Index 1 (i7)", "extract_and_demultiplex"),
    SourceRead("Index 2", "8 + 21 dark + 10", "t5 barcode + PCR i5",
               "Index 2 (i5)", "extract_and_demultiplex"),
    SourceRead("Read 2", "50", "genomic insert", "Read 2", "include"),
)


def protocol_graph() -> ProtocolGraph:
    protocol = PROTOCOL
    doi = SOURCE[0]
    source = "https://doi.org/" + doi
    graph = ProtocolGraph(protocol, authority=source)
    nuclei = graph.input(
        "nuclei", "accessible chromatin in nuclei", MaterialKind.BIOLOGICAL,
        state=Construct([seg("genomic DNA", "X" * 40, placeholder=True)]),
        primary=True)
    reviewed = ["nuclei"]
    p5_sequence = il.P5 + "N" * PCR_INDEX_NT + nx.S5
    p7_sequence = il.P7 + "N" * PCR_INDEX_NT + nx.S7
    p5_primer = PCRPrimer("P5 index primer", p5_sequence, anneal_nt=len(nx.S5))
    p7_primer = PCRPrimer("P7 index primer", p7_sequence, anneal_nt=len(nx.S7))
    transposome = graph.input(
        "indexed_transposomes", "indexed T5 and T7 transposomes",
        MaterialKind.ADAPTER, state=indexed_transposome("T5"),
        chemical_state=LoadedEndPool(
            (LoadedEnd("T5", "".join(seg.top for seg in S.t5())),
             LoadedEnd("T7", "".join(seg.top for seg in S.t7()))), nx.ME))
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagmented = graph.product(
        "tagmented_fragments", "heterologous T5/T7 tagmentation products",
        MaterialKind.MOLECULAR, state=two_tagmented())
    t5_t5 = graph.product("t5_t5", "T5/T5 tagmentation product",
                          MaterialKind.WASTE, state=S.tagmented_pair("T5", "T5"))
    t7_t7 = graph.product("t7_t7", "T7/T7 tagmentation product",
                          MaterialKind.WASTE, state=S.tagmented_pair("T7", "T7"))
    graph.react("tagment", "Tagment nuclei in the first indexing plate",
                StepKind.TAGMENTATION,
                ((nuclei, InputRole.TEMPLATE),
                 (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagmented, t5_t5, t7_t7),
                geometry=EndCombinationGeometry(
                    (p5_primer, p7_primer),
                    (EndProductGroup("tagmented_fragments", (("T5", "T7"),)),
                     EndProductGroup("t5_t5", (("T5", "T5"),)),
                     EndProductGroup("t7_t7", (("T7", "T7"),))),
                    "index_pcr"),
                relations=(MolecularRelation.END_COMBINATIONS,))
    graph.discard(t5_t5, "lacks the T7 PCR entry site")
    graph.discard(t7_t7, "lacks the T5 PCR entry site")
    pooled, pool_reviewed = pool_indexed_nuclei(
        graph, tagmented, prefix="round1", state=two_tagmented(),
        label="first-round indexed tagmented nuclei from other wells")
    gap_polymerase = graph.input("gap_polymerase", "gap-fill DNA polymerase",
                                 MaterialKind.CATALYST)
    final = FINAL_LIBRARY
    filled = graph.product(
        "gap_filled_fragments", "gap-filled redistributed fragments",
        MaterialKind.MOLECULAR, state=final)
    graph.react("gap_fill", "Fill the 9-nt tagmentation gaps",
                StepKind.EXTENSION,
                ((pooled, InputRole.TEMPLATE),
                 (gap_polymerase, InputRole.CATALYST)), (filled,))
    pcr_n = PCR_INDEX_NT
    p5 = graph.input(
        "p5_index_primer", "second-round P5 index primer", MaterialKind.OLIGO,
        state=Construct([seg("P5 index primer", p5_sequence, placeholder=True)]))
    p7 = graph.input(
        "p7_index_primer", "second-round P7 index primer", MaterialKind.OLIGO,
        state=Construct([seg("P7 index primer", p7_sequence, placeholder=True)]))
    polymerase = graph.input("pcr_polymerase", "PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(
        "library", f"{protocol} sequencing library", MaterialKind.LIBRARY,
        state=final)
    geometry = PCRGeometry.exact(
        final, p5_primer, p7_primer,
        statement="The 2018 oligo table prints both indexed PCR primers.")
    graph.react("index_pcr", "Add second-round P5/P7 indices by PCR",
                StepKind.PCR,
                ((filled, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=geometry)
    run_primers = PRIMERS
    graph.sequence(library, SequencingPlan(
        library=final, primers=run_primers,
        acquired_roles=tuple(p.role for p in run_primers)))
    reviewed.extend([
        "indexed_transposomes", "tn5", "tagment", *pool_reviewed,
        "gap_polymerase", "gap_fill", "p5_index_primer", "p7_index_primer",
        "pcr_polymerase", "index_pcr", "library",
    ])
    graph.declare_reviewed_scope(reviewed, source)
    graph.annotate(Highlight(
        "tagment", "Which Tn5 products enter the library",
        "The T5/T7 end pair supplies both second-round PCR primer sites. "
        "T5/T5 and T7/T7 products are separate discarded graph outputs drawn below; "
        "the declared end-combination check tests those primer-site assignments."))
    graph.annotate(DetailView(
        "tagmented_fragments", "Retained first-round fragment",
        "The T5/T7 product carries first-round indexing at both ends before gap fill."))
    return graph.require_publishable()


def oligo_section() -> OligoSection:
    """Source-specific oligos, assembled from local choices and canonical adapters."""
    rows = [
        OligoEntry("T5 indexed transposon", tuple(S.t5())),
        OligoEntry("T7 indexed transposon", tuple(S.t7())),
        OligoEntry("pMENTS", (seg("ME'", nx.ME_RC, "me"),), mods="5Phos"),
        OligoEntry("Read 1", (seg("extension A", S.R1_EXT, "r1"), seg("ME", nx.ME, "me"))),
        OligoEntry("Read 2", (seg("extension B", S.R2_EXT, "r1"), seg("ME", nx.ME, "me"))),
        OligoEntry("P5 PCR primer", (seg("P5", il.P5, "p5"),
                                     seg("i5", "N" * PCR_INDEX_NT, "cbc", placeholder=True),
                                     seg("s5", nx.S5, "s5"))),
        OligoEntry("P7 PCR primer", (seg("P7", il.P7, "p7"),
                                     seg("i7", "N" * PCR_INDEX_NT, "cbc", placeholder=True),
                                     seg("s7", nx.S7, "s7"))),
    ]
    return OligoSection("Key oligos", tuple(rows))


def render_page() -> str:
    graph = protocol_graph()
    preamble = (f'<div class="wrap" data-protocol-complete="true"><h1>{PROTOCOL} '
                f'&mdash; combinatorial single-cell chromatin accessibility</h1>\n'
                + info(f'Defining source: <a href="https://doi.org/{SOURCE[0]}">'
                       f'doi:{SOURCE[0]}</a>.'))
    return assemble_page(graph, [head(f"{PROTOCOL} library chemistry"),
                                 preamble, AnnotationSlot.INTRO,
                                 oligo_section(), WorkflowSection("Library generation"),
                                 AnnotationSlot.AFTER_WORKFLOW,
                                 ReadLayoutSection(
                                     "library", READ_ROWS,
                                     "Read roles and contents are source-declared; "
                                     "compound dark-cycle recipes remain in the layout "
                                     "table above.", '</div>')])


def selftest() -> None:
    from checks import Check
    check = Check()
    # Cusanovich et al. 2018 Supplementary Table S12, P5_1_PCR_Primer.
    check("P5_1_PCR_Primer", il.P5 + "CTCCATCGAG" + nx.S5, SOURCE_P5_1)
    check("the published ten-base PCR indices are established in the construct",
          [(s.name, s.inferred) for s in FINAL_LIBRARY
           if s.name in ("PCR i5", "PCR i7'")],
          [("PCR i5", False), ("PCR i7'", False)])
    page = render_page()
    check("both discarded Tn5 branches have molecular panels on the public page",
          [f'<section data-graph-target="{target}" data-graph-disposition="discard">'
           in page and
           '<svg' in page.split(f'<section data-graph-target="{target}"', 1)[1].split(
               '</section>', 1)[0]
           for target in ("t5_t5", "t7_t7")],
          [True, True])
    check.report()
