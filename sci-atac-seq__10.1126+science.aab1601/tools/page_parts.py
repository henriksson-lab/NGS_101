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
from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Construct, Scene, Segment, annotation_rows, complement_segments,
                      oligo, panel, protocol_graph_panels, revcomp, strand_row)
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
             source: tuple[str, str] | None = None) -> str:
    doi, title = source or PAPERS[protocol]
    if protocol == S.SCI15:
        caveat = ('<div class="caveat"><b>Inferred PCR indices.</b> The defining paper’s gated supplement contains the PCR-primer list. Dotted 8-nt PCR indices use the compatible Amini 2014 Nextera-kit layout; their exact 2015 length and bases remain unavailable. '
                  '<b>Published predecessor correction:</b> <a href="https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fng.3119/MediaObjects/41588_2014_BFng3119_MOESM2_ESM.xlsx">Amini 2014 Supplementary Table 4</a>, Correction Notice and corrected rows 1–5, reports that its original online table omitted the first five indexed Tn5 oligos; they were restored on 20 May 2015. Their corrected 8-base i5 values are TATAGCCT, ATAGAGGC, CCTATCCT, GGCTCTGA and AGGCGAAG. The model uses Amini’s shared Tn5 architecture with an unspecified first-well index; the gated 2015 supplement prevents assignment of these Amini members to this protocol.</div>')
    elif protocol == S.SCI3:
        caveat = ""
    else:
        caveat = ""
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


def sci3_oligos() -> str:
    inf = False
    rows = [
        oligo("N5 ligation oligo", [seg("head", S.N5_HEAD, "r2", inferred=inf), seg("N5 barcode", "A" * 10, "cbc", placeholder=True, inferred=inf), seg("tail", S.N5_TAIL, "r2", inferred=inf)]),
        oligo("N7 ligation oligo", [seg("head", S.N7_HEAD, "r3", inferred=inf), seg("N7 barcode", "B" * 10, "cbc", placeholder=True, inferred=inf), seg("tail", S.N7_TAIL, "r3", inferred=inf)]),
        oligo("P5 PCR primer", [seg("P5", il.P5, "p5"), seg("i5", "I" * 10, "cbc", placeholder=True), seg("N5 head", S.N5_HEAD, "r2", inferred=inf)]),
        oligo("P7 PCR primer", [seg("P7", il.P7, "p7"), seg("i7", "I" * 10, "cbc", placeholder=True), seg("N7 head", S.N7_HEAD, "r3", inferred=inf)]),
    ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def final_panel(protocol: str) -> str:
    lib = S.library(protocol)
    caption = "Final four-part cell-barcode library."
    if protocol == S.SCI15: caption = "INFERRED — final library using the compatible 8-nt PCR-index layout; dotted indices are unavailable in the defining supplement."
    if protocol == S.SCI3: caption = "Final three-level library assembled from Supplementary Table S7."
    return panel([strand_row(lib), strand_row(lib, "bottom"), *annotation_rows(lib)], caption=caption)


def two_steps(protocol: str) -> str:
    return f'''<h2>Library generation</h2>
<h3>(1) Assemble indexed T5 and T7 transposomes</h3>{panel(indexed_transposome("T5").rows(), cls="small")}{panel(indexed_transposome("T7").rows(), cls="small")}
<h3>(2) Tagment nuclei in the first indexing plate</h3>{panel(two_tagmented().rows(), caption="Representative amplifiable T5...T7 product before gap fill; each end records the first-round well.")}
<h3>(3) Pool nuclei, redistribute, fill gaps and PCR-index</h3>{info('A P5/P7 primer pair records the second-round well. Same-end T5...T5 and T7...T7 fragments lack one primer site.')}
<h3>(4) Final library</h3>{final_panel(protocol)}'''


def ligation_scene(which: str) -> Scene:
    if which == "N5":
        joined = [seg("head", S.N5_HEAD, "r2"), seg("barcode", "A" * 10, "cbc", placeholder=True), seg("tail", S.N5_TAIL, "r2"), seg("s5 first 8", nx.S5[:8], "s5"), seg("s5 rest", nx.S5[8:], "s5"), seg("ME", nx.ME, "me")]
        bridge = [seg("splint", revcomp(S.N5_TAIL + nx.S5[:8]), "r2")]
        pair = ("splint", "tail")
    else:
        joined = [seg("head", S.N7_HEAD, "r3"), seg("barcode", "B" * 10, "cbc", placeholder=True), seg("tail", S.N7_TAIL, "r3"), seg("s7 first 8", nx.S7[:8], "s7"), seg("s7 rest", nx.S7[8:], "s7"), seg("ME", nx.ME, "me")]
        bridge = [seg("splint", revcomp(S.N7_TAIL + nx.S7[:8]), "r3")]
        pair = ("splint", "tail")
    sc = Scene(); sc.strand("joined strand", joined); sc.anneal("splint", bridge, to="joined strand", pair=pair, unpaired=())
    sc.mark("joined strand", "tail", "ligation nick", through=("s5 first 8" if which == "N5" else "s7 first 8"))
    return sc


def ligation_reagent_scene(which: str) -> Scene:
    head, tail = ((S.N5_HEAD, S.N5_TAIL) if which == "N5"
                  else (S.N7_HEAD, S.N7_TAIL))
    adapter = [seg(f"{which} head", head),
               seg(f"{which} barcode", "N" * 10, "cbc", placeholder=True),
               seg(f"{which} tail", tail)]
    sc = Scene(); sc.strand("adapter", adapter, label=f"{which} ligation oligo")
    sc.ligation_end("adapter", f"{which} tail", kind="splinted",
                    polarity="3-prime", owner="adapter")
    sc.note("adapter", "end-specific splint aligns this 3′ end to the phosphorylated Tn5 end")
    return sc


def sci3_steps() -> str:
    return f'''<h2>Library generation</h2>
<h3>(1) Tagment nuclei with plain Nextera Tn5, then phosphorylate transferred 5' ends</h3>{info('T4 PNK creates the 5\'-phosphate required for ligation; no cell barcode is added by Tn5.')}
<h3>(2) First split-pool ligation: N5 barcode onto s5 ends</h3>{panel(ligation_scene("N5").rows(), cls="small", caption="An 8+8 splint bridges the N5 tail and s5.")}
<h3>(3) Second split-pool ligation: N7 barcode onto s7 ends</h3>{panel(ligation_scene("N7").rows(), cls="small", caption="The end-specific N7 splint analogously joins only s7 ends.")}
<h3>(4) Redistribute nuclei, fill gaps and add i5/i7 by PCR</h3>{info('Only an N5...N7 fragment carries both 15-nt PCR landing sites.')}
<h3>(5) Final library</h3>{final_panel(S.SCI3)}'''


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
        tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
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
        gap_polymerase = graph.input("gap_polymerase", "gap-fill DNA polymerase",
                                     MaterialKind.CATALYST)
        final = final_library or S.library(protocol)
        filled = graph.product(
            "gap_filled_fragments", "gap-filled redistributed fragments",
            MaterialKind.MOLECULAR, state=final, inferred=inferred,
            source=(source + " (2015 PCR-index layout inferred from compatible kit)"
                    if inferred else None))
        graph.react("gap_fill", "Fill the 9-nt tagmentation gaps",
                    StepKind.EXTENSION,
                    ((pooled, InputRole.TEMPLATE),
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
        polymerase = graph.input("pcr_polymerase", "PCR polymerase",
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
        run_primers = primers or S.TWO_LEVEL_PRIMERS
        graph.sequence(library, SequencingPlan(
            library=final, primers=run_primers,
            acquired_roles=tuple(p.role for p in run_primers)))
        reviewed.extend([
            "indexed_transposomes", "tn5", "tagment", *pool_reviewed,
            "gap_polymerase", "gap_fill", "p5_index_primer", "p7_index_primer",
            "pcr_polymerase", "index_pcr", "library",
        ])

    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def sequencing(protocol: str, *, final_library: Construct | None = None,
               read_rows: list[tuple[str, str, str]] | None = None,
               primers: tuple | None = None) -> str:
    lib = final_library or S.library(protocol)
    primers = primers or (S.SCI3_PRIMERS if protocol == S.SCI3 else S.TWO_LEVEL_PRIMERS)
    if sp.verify(lib, primers): raise ValueError("read layout does not match final library")
    if read_rows is not None:
        rows = read_rows
    elif protocol == S.SCI15:
        rows = [("Read 1", "not available", "genomic insert"), ("Index 1", "8 + 27 constant + PCR i7", "first- and second-round barcodes"), ("Index 2", "8 + 21 constant + PCR i5", "first- and second-round barcodes"), ("Read 2", "not available", "genomic insert")]
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
                primers: tuple | None = None) -> str:
    graph = graph_factory() if graph_factory else protocol_graph(protocol)
    body = "<h2>Library generation</h2>\n" + protocol_graph_panels(graph)
    olig = (two_oligos(protocol, pcr_n=pcr_n, inferred=inferred)
            if protocol != S.SCI3 else sci3_oligos())
    rendered = assemble_page(graph, [
        head(f"{protocol} library chemistry"),
        preamble(protocol, complete=True, source=source), AnnotationSlot.INTRO,
        olig, body, AnnotationSlot.AFTER_WORKFLOW,
        sequencing(protocol, final_library=final_library,
                   read_rows=read_rows, primers=primers)])
    if out is not None:
        out.write_text(rendered, encoding="utf-8")
        print(f"wrote {out}  ({out.stat().st_size:,} bytes)")
    return rendered
