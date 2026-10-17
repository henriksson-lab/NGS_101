#!/usr/bin/env python3
"""Shared renderer for individual sci-ATAC-seq protocols."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import sciatac as S
from sciatac import indexed_transposome, pool_indexed_nuclei, two_tagmented
from chemdraw import (Construct, Segment, annotation_rows, oligo, panel,
                      protocol_graph_panels, strand_row)
import illumina as il
import nextera as nx
from page import head, info, table
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind, AnnotationSlot)
from protocol_views import assemble_page
import seqprimers as sp

PAPERS = {
    S.SCI15: ("10.1126/science.aab1601", "sci-ATAC-seq"),
    S.SCI3: ("10.1126/science.aba7612", "sci-ATAC-seq3"),
}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def preamble(protocol: str, *, complete: bool = False,
             source: tuple[str, str] | None = None, caveat: str | None = None) -> str:
    doi, title = source or PAPERS[protocol]
    if caveat is None and protocol == S.SCI15:
        caveat = ('<div class="caveat"><b>Inferred PCR indices.</b> The 2015 Supplementary Materials (Materials and Methods, “Sample Processing”) amplify 20 cycles “while incorporating standard Nextera library barcodes following the protocol described in” Adey 2014 and Amini 2014, and state that after the first three experiments the authors “switched to a set of 10bp library barcodes that were otherwise identical to the standard Nextera library adapters”. No PCR-primer bases are printed. Dotted 8-nt PCR indices draw the standard Nextera-kit variant of the first three experiments (8-bp length and four-part barcode order also in the authors’ GEO record GSM1647121); the later 10-nt variant is not drawn. Tn5 release, gap fill and PCR mix follow the cited Amini 2014 protocol by reference. '
                  '<b>Inferred Index 2 priming.</b> The supplement sequences on a MiSeq (v2 300-cycle kit, paired 51-bp reads, custom recipe) with custom primers “as described in” Adey 2014 and Amini 2014; Amini prints no Index 2 primer, so the i5 read is drawn from the flow-cell P5 oligo, reading PCR i5 before the Tn5 i5. Index-read cycle numbers are not stated. '
                  '<b>Published predecessor correction:</b> <a href="https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fng.3119/MediaObjects/41588_2014_BFng3119_MOESM2_ESM.xlsx">Amini 2014 Supplementary Table 4</a>, Correction Notice and corrected rows 1–5, reports that its original online table omitted the first five indexed Tn5 oligos; they were restored on 20 May 2015. Their corrected 8-base i5 values are TATAGCCT, ATAGAGGC, CCTATCCT, GGCTCTGA and AGGCGAAG. The model uses Amini’s shared Tn5 architecture with an unspecified first-well index; the gated 2015 supplement prevents assignment of these Amini members to this protocol.</div>')
    caveat = caveat or ""
    complete_attr = ' data-protocol-complete="true"' if complete else ""
    return f'''<div class="wrap"{complete_attr}><h1>{title} &mdash; combinatorial single-cell chromatin accessibility</h1>
{info(f'Defining source: <a href="https://doi.org/{doi}">doi:{doi}</a>.')}{caveat}'''


def two_oligos(protocol: str, *, pcr_n: int | None = None,
               inferred: bool | None = None) -> str:
    pcr_n = 8 if pcr_n is None else pcr_n
    inf = protocol == S.SCI15 if inferred is None else inferred
    rows = [
        oligo("T5 indexed transposon", S.t5()), oligo("T7 indexed transposon", S.t7()),
        oligo("pMENTS", [seg("ME'", nx.ME_RC, "me")], mods="5Phos"),
        oligo("Read 1", [seg("extension A", S.R1_EXT, "r1"), seg("ME", nx.ME, "me")]),
        oligo("Read 2", [seg("extension B", S.R2_EXT, "r1"), seg("ME", nx.ME, "me")]),
        oligo("P5 PCR primer", [seg("P5", il.P5, "p5"), seg("i5", "N" * pcr_n, "cbc", placeholder=True, inferred=inf), seg("s5", nx.S5, "s5")]),
        oligo("P7 PCR primer", [seg("P7", il.P7, "p7"), seg("i7", "N" * pcr_n, "cbc", placeholder=True, inferred=inf), seg("s7", nx.S7, "s7")]),
    ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def final_panel(protocol: str) -> str:
    lib = S.library(protocol)
    caption = "Final four-part cell-barcode library."
    if protocol == S.SCI15: caption = "INFERRED — final library using the compatible 8-nt PCR-index layout; dotted indices are unavailable in the defining supplement."
    if protocol == S.SCI3: caption = "Final three-level library assembled from Supplementary Table S7."
    return panel([strand_row(lib), strand_row(lib, "bottom"), *annotation_rows(lib)], caption=caption)


def protocol_graph(protocol: str, *, doi: str | None = None,
                   final_library: Construct | None = None,
                   pcr_n: int | None = None,
                   inferred: bool | None = None,
                   primers: tuple | None = None) -> ProtocolGraph:
    if protocol not in PAPERS and doi is None:
        raise ValueError(f"unknown sci-ATAC protocol: {protocol!r}")
    source = "https://doi.org/" + (doi or PAPERS[protocol][0])
    graph = ProtocolGraph(protocol, authority=source)
    nuclei = graph.input(
        "nuclei", "accessible chromatin in nuclei", MaterialKind.BIOLOGICAL,
        state=Construct([seg("genomic DNA", "X" * 40, placeholder=True)]),
        primary=True)
    reviewed = ["nuclei"]

    if protocol != S.SCI3:
        inferred = protocol == S.SCI15 if inferred is None else inferred
        transposome = graph.input(
            "indexed_transposomes", "indexed T5 and T7 transposomes",
            MaterialKind.ADAPTER, state=indexed_transposome("T5"),
            inferred=inferred,
            source=(source + " (2015 PCR-index supplement unavailable)"
                    if inferred else None))
        tn5 = graph.input("tn5", ("indexed Tn5 2.5 µM in Nextera TD buffer, 2,500 cells per "
                                  "well, 37 °C 30 min; stop 40 mM EDTA, 37 °C 15 min (2015 SM)"
                                  if protocol == S.SCI15 else "Tn5 transposase"),
                          MaterialKind.CATALYST)
        tagmented = graph.product(
            "tagmented_fragments", "heterologous T5/T7 tagmentation products",
            MaterialKind.MOLECULAR, state=two_tagmented(), inferred=inferred,
            source=(source + " (2015 PCR-index layout inferred from compatible kit)"
                    if inferred else None))
        tagment_waste = graph.product("tagment_waste", "same-end tagmentation products",
                                      MaterialKind.WASTE)
        graph.react("tagment", "Tagment nuclei in the first indexing plate",
                    StepKind.TAGMENTATION,
                    ((nuclei, InputRole.TEMPLATE),
                     (transposome, InputRole.STRUCTURAL),
                     (tn5, InputRole.CATALYST)), (tagmented, tagment_waste))
        graph.discard(tagment_waste,
                      "same-end T5/T5 and T7/T7 products lack one terminal PCR site")
        pooled, pool_reviewed = pool_indexed_nuclei(
            graph, tagmented, prefix="round1", state=two_tagmented(),
            label="first-round indexed tagmented nuclei from other wells")
        # Science 2015 main text (PMC4836442 JATS line 45): pool, dilute,
        # redistribute 15-25 nuclei per second-plate well by cell sorter, lyse, PCR.
        sorter = graph.input(
            "second_sort", ("FACSAria II: 15 nuclei (FSC/SSC) or 25 DAPI-positive nuclei into "
                            "20 µl EB per well (2015 SM)" if protocol == S.SCI15 else
                            "cell sorter: 15–25 nuclei into each of 96 second-plate wells"),
            MaterialKind.SELECTION)
        redistributed = graph.product(
            "redistributed_nuclei", "15–25 pooled nuclei per second-round well",
            MaterialKind.MOLECULAR, state=two_tagmented())
        unsorted = graph.product(
            "unsorted_nuclei", "pooled nuclei not sorted into a second-round well",
            MaterialKind.WASTE)
        graph.react("redistribute", "Dilute and re-sort pooled nuclei into second-round wells",
                    StepKind.SAMPLING,
                    ((pooled, InputRole.TEMPLATE), (sorter, InputRole.SELECTION)),
                    (redistributed, unsorted), changes_molecule=False)
        graph.discard(unsorted, "not sorted into an assayed second-round well")
        lysis = graph.input(
            "lysis", ("Tn5 release per cited Amini 2014 protocol (0.1% SDS, 55 °C 15 min); "
                      "2015 SM names no conditions" if protocol == S.SCI15 else
                      "nuclear lysis (conditions not modelled)"),
            MaterialKind.CONDITION)
        released = graph.product(
            "released_fragments", "tagmented DNA released in each second-round well",
            MaterialKind.MOLECULAR, state=two_tagmented())
        graph.react("lyse", "Lyse sorted nuclei", StepKind.EXTRACTION,
                    ((redistributed, InputRole.TEMPLATE), (lysis, InputRole.CONDITION)),
                    (released,))
        gap_polymerase = graph.input("gap_polymerase", ("NPM 72 °C 3 min gap fill (cited Amini 2014 protocol)"
                                                        if protocol == S.SCI15 else "gap-fill DNA polymerase"),
                                     MaterialKind.CATALYST)
        final = final_library or S.library(protocol)
        filled = graph.product(
            "gap_filled_fragments", "gap-filled redistributed fragments",
            MaterialKind.MOLECULAR, state=S.gap_filled_fragment(), inferred=inferred,
            source=(source + " (2015 PCR-index layout inferred from compatible kit)"
                    if inferred else None))
        graph.react("gap_fill", "Fill the 9-nt tagmentation gaps",
                    StepKind.EXTENSION,
                    ((released, InputRole.TEMPLATE),
                     (gap_polymerase, InputRole.CATALYST)), (filled,))
        pcr_n = 8 if pcr_n is None else pcr_n
        p5_sequence = il.P5 + "N" * pcr_n + nx.S5
        p7_sequence = il.P7 + "N" * pcr_n + nx.S7
        p5 = graph.input(
            "p5_index_primer", "second-round P5 index primer", MaterialKind.OLIGO,
            state=Construct([seg("P5 index primer", p5_sequence, placeholder=True)]),
            inferred=inferred,
            source=(source + " (2015 primer sequence unavailable)" if inferred else None))
        p7 = graph.input(
            "p7_index_primer", "second-round P7 index primer", MaterialKind.OLIGO,
            state=Construct([seg("P7 index primer", p7_sequence, placeholder=True)]),
            inferred=inferred,
            source=(source + " (2015 primer sequence unavailable)" if inferred else None))
        polymerase = graph.input("pcr_polymerase", ("20 cycles (2015 SM); NPM per cited Amini 2014 protocol"
                                                    if protocol == S.SCI15 else "PCR polymerase"),
                                 MaterialKind.CATALYST)
        library = graph.product(
            "library", f"{protocol} sequencing library", MaterialKind.LIBRARY,
            state=final, inferred=inferred,
            source=(source + " (2015 PCR-index layout inferred from compatible kit)"
                    if inferred else None))
        geometry = (PCRGeometry.unavailable(
            final, reason="The gated 2015 supplement prevents exact PCR-primer transcription.")
                    if inferred else
                    PCRGeometry.exact(
                        final, PCRPrimer("P5 index primer", p5_sequence),
                        PCRPrimer("P7 index primer", p7_sequence),
                        statement="The 2018 oligo table prints both indexed PCR primers."))
        graph.react("index_pcr", "Add second-round P5/P7 indices by PCR",
                    StepKind.PCR,
                    ((filled, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                     (p7, InputRole.STRUCTURAL),
                     (polymerase, InputRole.CATALYST)), (library,),
                    geometry=geometry)
        run_primers = primers or (S.SCI15_PRIMERS if protocol == S.SCI15
                                  else S.TWO_LEVEL_PRIMERS)
        graph.sequence(library, SequencingPlan(
            library=final, primers=run_primers,
            acquired_roles=tuple(p.role for p in run_primers)))
        reviewed.extend([
            "indexed_transposomes", "tn5", "tagment", *pool_reviewed,
            "second_sort", "redistribute", "lysis", "lyse", "gap_polymerase", "gap_fill", "p5_index_primer", "p7_index_primer",
            "pcr_polymerase", "index_pcr", "library",
        ])

    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def sequencing(protocol: str, *, final_library: Construct | None = None,
               read_rows: list[tuple[str, str, str]] | None = None,
               primers: tuple | None = None) -> str:
    lib = final_library or S.library(protocol)
    primers = primers or (S.SCI3_PRIMERS if protocol == S.SCI3 else
                          S.SCI15_PRIMERS if protocol == S.SCI15 else S.TWO_LEVEL_PRIMERS)
    if sp.verify(lib, primers): raise ValueError("read layout does not match final library")
    if read_rows is not None:
        rows = read_rows
    elif protocol == S.SCI15:
        rows = [("Read 1", "51", "genomic insert"), ("Index 1", "8 + 27 constant + PCR i7 (cycles not stated)", "first- and second-round barcodes"), ("Index 2", "PCR i5 + 21 constant + 8 (cycles not stated)", "second- and first-round barcodes"), ("Read 2", "51", "genomic insert")]
    else:
        rows = [("Read 1", "51", "genomic insert"), ("Index 1", "10 + 15 dark + 10", "N7 + i7"), ("Index 2", "10 + 15 dark + 10", "N5 + i5"), ("Read 2", "51", "genomic insert")]
    specs = []
    for role, cycles, content in rows:
        role = {"Index 1": "Index 1 (i7)", "Index 2": "Index 2 (i5)"}.get(role, role)
        fixed = int(cycles) if str(cycles).isdigit() else None
        use = "include" if role.startswith("Read") else "extract_and_demultiplex"
        specs.append(sp.read_spec(role, fixed, content.replace(" + ", "; "), use))
    return ('<h2>Read layout</h2>' + sp.diagram(lib, primers)
            + table(("Read", "Cycles", "Content"), rows)
            + sp.read_spec_table(specs, "Read roles and contents are source-declared; compound dark-cycle recipes remain in the layout table above.")
            + '</div>')


def render_page(protocol: str, out: Path | None = None, *,
                source: tuple[str, str] | None = None,
                graph_factory=None,
                final_library: Construct | None = None,
                pcr_n: int | None = None,
                inferred: bool | None = None,
                read_rows: list[tuple[str, str, str]] | None = None,
                primers: tuple | None = None,
                oligo_html: str | None = None,
                caveat: str | None = None) -> str:
    """Render one page; a protocol whose oligos are not the two-level T5/T7 set
    supplies its own Key-oligos markup (and any preamble caveat)."""
    graph = graph_factory() if graph_factory else protocol_graph(protocol)
    body = "<h2>Library generation</h2>\n" + protocol_graph_panels(graph)
    if oligo_html is None and protocol == S.SCI3:
        raise ValueError("sci-ATAC-seq3 supplies its Table S7 oligo list")
    olig = (oligo_html if oligo_html is not None
            else two_oligos(protocol, pcr_n=pcr_n, inferred=inferred))
    rendered = assemble_page(graph, [
        head(f"{protocol} library chemistry"),
        preamble(protocol, complete=True, source=source, caveat=caveat), AnnotationSlot.INTRO,
        olig, body, AnnotationSlot.AFTER_WORKFLOW,
        sequencing(protocol, final_library=final_library,
                   read_rows=read_rows, primers=primers)])
    if out is not None:
        out.write_text(rendered, encoding="utf-8")
        print(f"wrote {out}  ({out.stat().st_size:,} bytes)")
    return rendered
