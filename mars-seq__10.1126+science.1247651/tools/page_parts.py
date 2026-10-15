"""Sparse individual page renderers for MARS-seq and MARS-seq2.0."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import marsseq as M
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Segment, oligo, panel, protocol_graph_panels, strand_row
from page import head, info, table
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from rna_amplification import barcoded_ivt_front


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


PAPERS = {
    "MARS-seq": ("10.1126/science.1247651", "Science", "2014"),
    "MARS-seq2.0": ("10.1038/s41596-019-0164-4", "Nature Protocols", "2019"),
}


def protocol_graph(protocol: str) -> ProtocolGraph:
    if protocol not in PAPERS:
        raise ValueError(f"unknown MARS-seq protocol: {protocol!r}")
    source = f"https://doi.org/{PAPERS[protocol][0]}"
    graph = ProtocolGraph(protocol, authority=source)
    primary = graph.input("mrna", "polyadenylated RNA from one well",
                          MaterialKind.BIOLOGICAL, state=M.rt_scene(protocol),
                          primary=True)
    arna, reviewed = barcoded_ivt_front(
        graph, primary, rt_primer_label=f"{protocol} RT1",
        rt_state=M.rt_scene(protocol), first_state=M.rt_scene(protocol),
        duplex_state=M.ds_cdna(protocol), arna_state=M.amplified_rna(protocol))
    reviewed.insert(0, "mrna")

    fragmentation = graph.input("fragmentation_condition",
                                "RNA fragmentation condition",
                                MaterialKind.CONDITION)
    fragmented = graph.product("fragmented_arna", "fragmented amplified RNA",
                               MaterialKind.MOLECULAR,
                               state=M.amplified_rna(protocol))
    graph.react("fragment_arna", "Fragment amplified RNA",
                StepKind.FRAGMENTATION,
                ((arna, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CONDITION)), (fragmented,))

    adapter_state = (M.mars2_ligation_adapter() if protocol == "MARS-seq2.0"
                     else Construct([
                         seg("ordered pool barcode", M.POOL_ORDERED, "cbc"),
                         seg("random diversity", "N" * M.MARS1_POOL_RANDOM_NT,
                             placeholder=True),
                         seg("ligation/RT2 handle", M.LIG_CONSTANT, "r1"),
                     ], name="MARS-seq ligation adapter"))
    adapter = graph.input("pool_adapter", "pool-barcoded RNA ligation adapter",
                          MaterialKind.ADAPTER, state=adapter_state)
    ligase = graph.input("rna_ligase", "RNA ligase", MaterialKind.CATALYST)
    ligated = graph.product("ligated_arna", "pool-adapter-ligated amplified RNA",
                            MaterialKind.MOLECULAR,
                            state=M.ligated_arna(protocol))
    geometry = AdapterLigationGeometry.unavailable(
        "The source prints the adapter sequence and terminal modifications, but the "
        "compact model does not establish both ligation-facing RNA end geometries.")
    graph.react("ligate_pool_adapter", "Ligate the pool-barcoded RNA adapter",
                StepKind.LIGATION,
                ((fragmented, InputRole.TEMPLATE),
                 (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,), geometry=geometry)

    if protocol == "MARS-seq":
        rt2 = graph.input("rt2", "RT2 primer", MaterialKind.OLIGO,
                          state=Construct([seg("RT2", M.MARS1_RT2_PRINTED)]))
    else:
        rt2 = graph.input(
            "rt2", "RT2 primer (sequence not disclosed)", MaterialKind.OLIGO,
            inferred=True, source=source + " (RT2 sequence not disclosed)")
    reverse_transcriptase = graph.input("rt2_reverse_transcriptase",
                                       "RT2 reverse transcriptase",
                                       MaterialKind.CATALYST)
    final = M.final_library(protocol)
    inferred_boundary = protocol == "MARS-seq2.0"
    pcr_template = graph.product(
        "pcr_template", "second-strand library-PCR template",
        MaterialKind.MOLECULAR, state=final, inferred=inferred_boundary,
        source=(source + " (RT2 and PCR-primer sequences not disclosed)"
                if inferred_boundary else source))
    graph.react("second_reverse_transcription",
                "Reverse-transcribe the adapter-ligated amplified RNA with RT2",
                StepKind.REVERSE_TRANSCRIPTION,
                ((ligated, InputRole.TEMPLATE), (rt2, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (pcr_template,))

    if protocol == "MARS-seq":
        f_primer = PCRPrimer("P5_Rd1 PCR", M.MARS1_P5_PCR)
        r_primer = PCRPrimer("P7_Rd2 PCR", M.MARS1_P7_PCR)
        forward = graph.input("pcr_forward", f_primer.name, MaterialKind.OLIGO,
                              state=Construct([seg(f_primer.name, f_primer.sequence)]))
        reverse = graph.input("pcr_reverse", r_primer.name, MaterialKind.OLIGO,
                              state=Construct([seg(r_primer.name, r_primer.sequence)]))
        pcr_geometry = PCRGeometry.exact(
            final, f_primer, r_primer,
            statement="Supplementary Tables S7–S8 print both library-PCR primers.")
    else:
        primer_source = source + " (library-PCR primer sequences not disclosed)"
        forward = graph.input("pcr_forward", "library-PCR forward primer",
                              MaterialKind.OLIGO, inferred=True, source=primer_source)
        reverse = graph.input("pcr_reverse", "library-PCR reverse primer",
                              MaterialKind.OLIGO, inferred=True, source=primer_source)
        pcr_geometry = PCRGeometry.unavailable(
            final, reason="The defining source does not disclose both PCR primer sequences.")
    polymerase = graph.input("pcr_polymerase", "PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(
        "library", f"{protocol} sequencing library", MaterialKind.LIBRARY,
        state=final, inferred=inferred_boundary,
        source=(source + " (outer PCR arms not disclosed)"
                if inferred_boundary else source))
    graph.react("library_pcr", "Amplify the final sequencing library",
                StepKind.PCR,
                ((pcr_template, InputRole.TEMPLATE),
                 (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=pcr_geometry)
    if protocol == "MARS-seq":
        graph.sequence(library, SequencingPlan(
            library=final, primers=M.SEQ_PRIMERS,
            acquired_roles=("Read 1", "Read 2")))
    else:
        graph.sequence(library, SequencingPlan(
            limitation=("The source declares Read 1 and Read 2 contents, but does not "
                        "disclose exact sequencing-primer sites.")))
    reviewed.extend((
        "fragmentation_condition", "fragment_arna", "pool_adapter", "rna_ligase",
        "ligate_pool_adapter", "rt2", "rt2_reverse_transcriptase",
        "second_reverse_transcription", "pcr_forward", "pcr_reverse",
        "pcr_polymerase", "library_pcr", "library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def render_page(protocol: str) -> str:
    graph = protocol_graph(protocol)
    doi, journal, year = PAPERS[protocol]
    author = "Jaitin et al." if protocol == "MARS-seq" else "Keren-Shaul et al."
    exact = protocol == "MARS-seq2.0"
    caveat = ("The 2014 author supplement publishes RT1, the ligation adapters and both "
              "library-PCR primers; the final library below is assembled from those oligos."
              if not exact else
              "The supplements publish RT1 and the pool-barcoded ligation adapters, but not "
              "RT2 or the PCR-primer sequences; dotted outer arms stop at that boundary.")
    oligos = [oligo("RT1", M.rt1_segments(protocol))]
    if exact:
        oligos.append(oligo("representative pool ligation adapter",
                            list(M.mars2_ligation_adapter()), mods="/5Phos/ … /3SpC3/"))
    else:
        oligos.extend([
            oligo("representative pool ligation adapter",
                  [seg("pool barcode", M.POOL_ORDERED, "cbc"),
                   seg("random diversity", "N" * M.MARS1_POOL_RANDOM_NT,
                       placeholder=True),
                   seg("ligation/RT2 handle", M.LIG_CONSTANT, "r1")],
                  mods="/5Phos/ … /3SpC3/"),
            oligo("P5_Rd1 PCR", [seg("P5 + Read 1", M.MARS1_P5_PCR, "r1")]),
            oligo("P7_Rd2 PCR", [seg("P7 + Read 2", M.MARS1_P7_PCR, "r2")]),
        ])
    final = M.final_library(protocol)
    if exact:
        r1 = Construct([seg("random diversity", "N" * 5, placeholder=True),
                        seg("pool barcode", M.POOL_READ, "cbc"),
                        seg("cDNA", "XXXXXXXX...", placeholder=True)])
        r2 = Construct([seg("well barcode", "N" * 7, "cbc", placeholder=True),
                        seg("UMI", "N" * 8, "umi", placeholder=True),
                        seg("poly(T), ignored", "TTTTT")])
        # The Supplementary Information's read-design fields print the worked example
        # R1_design = 3I.4P.50M with R2_design = 7W.8R.5I.  Read 1's 3I matches the
        # MARS-seq1 N3 adapter, so 5I is what the 2.0 N5 adapter implies -- derived here,
        # never published, and labelled as such.
        read_caption = ("Published MARS-seq2.0 read design: 3I.4P.50M and 7W.8R.5I. "
                        "DERIVED — with this version's N5 ligation adapter Read 1 is "
                        "5I.4P instead; the 5I is derived from the adapter, not published.")
        read_rows = [("Read 1", "N5, 4-nt pool barcode, sense cDNA"),
                     ("Read 2", "7-nt well barcode, 8-nt UMI, poly(T)")]
    else:
        r1 = Construct([seg("random diversity", "N" * 3, placeholder=True),
                        seg("pool barcode", M.POOL_READ, "cbc"),
                        seg("cDNA", "XXXXXXXX...", placeholder=True)])
        r2 = Construct([seg("well barcode", "N" * 6, "cbc", placeholder=True),
                        seg("UMI", "N" * 4, "umi", placeholder=True),
                        seg("poly(T), ignored", "TTTTT")])
        read_caption = "Published MARS-seq read design: N3, 4-nt pool barcode, then cDNA; Read 2 reports cell barcode and UMI."
        read_rows = [("Read 1", "N3, 4-nt pool barcode, sense cDNA"),
                     ("Read 2", "6-nt well barcode, 4-nt UMI, poly(T)")]
    return "\n".join([
        head(f"{protocol} library chemistry"),
        '<div class="wrap" data-protocol-complete="true">',
        f'<h1>{protocol} &mdash; three-level barcoded RNA-seq with T7 amplification</h1>',
        info(f'Defining source: <a href="https://doi.org/{doi}">{author}, '
             f'<i>{journal}</i> ({year})</a>.'),
        f'<div class="caveat"><b>Sequence boundary.</b> {caveat}</div>',
        f'<h2>Key oligos</h2><seq>{"".join(oligos)}</seq>',
        '<h2>Library construction</h2>', protocol_graph_panels(graph),
        (sp.unavailable_diagram(final, "the defining sources do not print the library-primer sequences")
         if exact else sp.diagram(final, M.SEQ_PRIMERS)),
        '<h2>Read layout</h2>',
        panel([strand_row(r1, prefix="Read 1  ", suffix=""),
               strand_row(r2, prefix="Read 2  ", suffix="")], cls="small",
              caption=read_caption),
        table(("Read", "Reports"), read_rows),
        sp.read_spec_table([
            sp.read_spec("Read 1", None, read_rows[0][1], "extract_then_include"),
            sp.read_spec("Read 2", None, read_rows[1][1], "extract"),
        ], "Read contents are source-declared; the defining sources do not fix a complete run recipe."),
        '</div>',
    ])
