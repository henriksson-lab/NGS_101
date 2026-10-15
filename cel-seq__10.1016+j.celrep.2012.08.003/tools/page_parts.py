"""Sparse, source-specific renderers for CEL-Seq and CEL-Seq2."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import celseq as C
import seqprimers as sp
from chemdraw import Construct, Segment, oligo, panel, protocol_graph_panels, strand_row
from page import head, info, table
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from rna_amplification import barcoded_ivt_front


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


PAPERS = {
    "CEL-Seq": ("10.1016/j.celrep.2012.08.003", "Cell Reports", "2012"),
}


def preamble(protocol: str, *, complete: bool = False,
             source: tuple[str, str, str] | None = None) -> str:
    doi, journal, year = source or PAPERS[protocol]
    supporting = ""
    if protocol == "CEL-Seq":
        supporting = (' Supporting primer source: <a href="https://doi.org/'
                      '10.1186/s13059-016-0938-8">CEL-Seq2 Table S2</a>.')
    complete_attr = ' data-protocol-complete="true"' if complete else ""
    return f"""<div class="wrap"{complete_attr}>
<h1>{protocol} &mdash; barcoded 3' RNA-seq by T7 linear amplification</h1>
{info(f'Defining source: <a href="https://doi.org/{doi}">Hashimshony et al., '
      f'<i>{journal}</i> ({year})</a>.{supporting}')}
<div class="caveat"><b>Small-RNA library boundary.</b> The accessible primary sources do
not print the Illumina library-primer sequences. Dotted outer arms are unpublished
placeholders, not the secondary-source TruSeq Small RNA reconstruction.</div>"""


def oligos(protocol: str, design: C.PrimerDesign | None = None) -> str:
    rows = [oligo(f"{protocol} RT primer", C.rt_primer_segments(protocol, design))]
    if protocol == "CEL-Seq2":
        if design is None:
            raise ValueError("CEL-Seq2 requires its source-specific primer design")
        rows.append(oligo("randomhexRT primer", [
            seg("library-RT tail", design.randomhex_tail, "r2"),
            seg("random hexamer", "N" * design.randomhex_nt, placeholder=True),
        ]))
    return f"<h2>Key oligos</h2><seq>{''.join(rows)}</seq>"


def cel1_readout() -> str:
    final = C.final_library("CEL-Seq")
    reads = [
        ("Read 1", "cell barcode, then poly(T); exact primer sequence unavailable"),
        ("Index", "library index; exact primer and index design unavailable"),
        ("Read 2", "transcript sequence from the 3' tag"),
    ]
    return f"""<h2>Read layout</h2>
{sp.unavailable_diagram(final, "exact library-primer sequences are unavailable")}
{table(("Read", "Reports"), reads)}
{sp.read_spec_table([
    sp.read_spec("Read 1", None, "cell barcode, then poly(T)", "extract"),
    sp.read_spec("Index 1 (i7)", None, "library index", "demultiplex"),
    sp.read_spec("Read 2", None, "transcript sequence from the 3' tag", "include"),
], "Read roles are source-declared; exact primer sequences and cycle counts are unavailable.")}"""


def protocol_graph(protocol: str, *, design: C.PrimerDesign | None = None,
                   doi: str | None = None) -> ProtocolGraph:
    if protocol not in PAPERS and doi is None:
        raise ValueError(f"unknown CEL-Seq protocol: {protocol!r}")
    doi = doi or PAPERS[protocol][0]
    source = f"https://doi.org/{doi}"
    graph = ProtocolGraph(protocol, authority=source)

    primary = graph.input("mrna", "polyadenylated RNA from one well",
                          MaterialKind.BIOLOGICAL, state=C.rt_scene(protocol, design),
                          primary=True)
    arna, reviewed = barcoded_ivt_front(
        graph, primary, rt_primer_label=f"{protocol} RT primer",
        rt_state=C.rt_scene(protocol, design), first_state=C.rt_scene(protocol, design),
        duplex_state=C.double_stranded_cdna(protocol, design),
        arna_state=C.arna(protocol, design))
    reviewed.insert(0, "mrna")

    shear = graph.input("fragmentation_condition", "RNA fragmentation condition",
                        MaterialKind.CONDITION)
    fragmented = graph.product("fragmented_arna", "fragmented amplified RNA",
                               MaterialKind.MOLECULAR, state=C.arna(protocol, design))
    graph.react("fragment_arna", "Fragment amplified RNA", StepKind.FRAGMENTATION,
                ((arna, InputRole.TEMPLATE), (shear, InputRole.CONDITION)),
                (fragmented,))
    reviewed.extend(("fragmentation_condition", "fragment_arna"))

    final_state = C.final_library(protocol, design)
    if protocol == "CEL-Seq":
        adapter = graph.input(
            "rna_adapter", "RNA 3′ adapter (sequence not disclosed)",
            MaterialKind.ADAPTER, state=C.cel1_adapter_geometry(), inferred=True,
            source=source + " (adapter bases and terminal geometry not disclosed)")
        ligase = graph.input("rna_ligase", "RNA ligase", MaterialKind.CATALYST)
        ligated = graph.product("ligated_arna", "3′-adapter-ligated aRNA",
                                MaterialKind.MOLECULAR,
                                state=C.cel1_ligated_arna(), inferred=True,
                                source=source + " (adapter sequence not disclosed)")
        graph.react("ligate_adapter", "Ligate the RNA 3′ adapter",
                    StepKind.LIGATION,
                    ((fragmented, InputRole.TEMPLATE),
                     (adapter, InputRole.STRUCTURAL),
                     (ligase, InputRole.CATALYST)), (ligated,),
                    geometry=C.cel1_adapter_geometry())
        library_rt_primer = graph.input(
            "library_rt_primer", "library reverse-transcription primer",
            MaterialKind.OLIGO, inferred=True,
            source=source + " (exact library RT primer not disclosed)")
        library_rt = graph.input("library_reverse_transcriptase",
                                 "library reverse transcriptase",
                                 MaterialKind.CATALYST)
        pcr_template = graph.product(
            "pcr_template", "library-PCR template", MaterialKind.MOLECULAR,
            state=final_state, inferred=True,
            source=source + " (outer library arms not disclosed)")
        graph.react("library_reverse_transcription",
                    "Reverse-transcribe the adapter-ligated aRNA",
                    StepKind.REVERSE_TRANSCRIPTION,
                    ((ligated, InputRole.TEMPLATE),
                     (library_rt_primer, InputRole.STRUCTURAL),
                     (library_rt, InputRole.CATALYST)), (pcr_template,))
        reviewed.extend(("rna_adapter", "rna_ligase", "ligate_adapter",
                         "library_rt_primer", "library_reverse_transcriptase",
                         "library_reverse_transcription"))
    else:
        randomhex = graph.input("randomhex_rt", "randomhexRT primer",
                                MaterialKind.OLIGO,
                                state=C.random_rt_scene(design))
        library_rt = graph.input("library_reverse_transcriptase",
                                 "library reverse transcriptase",
                                 MaterialKind.CATALYST)
        pcr_template = graph.product(
            "pcr_template", "random-primed library-PCR template",
            MaterialKind.MOLECULAR, state=final_state, inferred=True,
            source=source + " (outer library arms not disclosed)")
        graph.react("library_reverse_transcription",
                    "Reverse-transcribe fragmented aRNA with randomhexRT",
                    StepKind.REVERSE_TRANSCRIPTION,
                    ((fragmented, InputRole.TEMPLATE),
                     (randomhex, InputRole.STRUCTURAL),
                     (library_rt, InputRole.CATALYST)), (pcr_template,))
        reviewed.extend(("randomhex_rt", "library_reverse_transcriptase",
                         "library_reverse_transcription"))

    forward = graph.input(
        "pcr_forward", "library-PCR forward primer (sequence not disclosed)",
        MaterialKind.OLIGO, inferred=True,
        source=source + " (library-PCR primer sequences not disclosed)")
    reverse = graph.input(
        "pcr_reverse", "indexed library-PCR reverse primer (sequence not disclosed)",
        MaterialKind.OLIGO, inferred=True,
        source=source + " (library-PCR primer sequences not disclosed)")
    polymerase = graph.input("pcr_polymerase", "PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(
        "library", f"{protocol} sequencing library", MaterialKind.LIBRARY,
        state=final_state, inferred=True,
        source=source + " (outer library arms not disclosed)")
    graph.react("library_pcr", "Amplify and index the sequencing library",
                StepKind.PCR,
                ((pcr_template, InputRole.TEMPLATE),
                 (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    final_state,
                    reason="The defining sources do not disclose both library-PCR primers."))
    graph.sequence(library, SequencingPlan(
        limitation=("The published read roles are retained below, but the defining "
                    "sources do not disclose exact sequencing-primer sites.")))
    reviewed.extend(("pcr_forward", "pcr_reverse", "pcr_polymerase",
                     "library_pcr", "library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def render_page(protocol: str, *, design: C.PrimerDesign | None = None,
                source: tuple[str, str, str] | None = None,
                graph_factory=None, readout=None) -> str:
    graph = graph_factory() if graph_factory else protocol_graph(protocol)
    readout_html = readout() if readout else cel1_readout()
    workflow = "<h2>Library construction</h2>\n" + protocol_graph_panels(graph)
    return "\n".join([head(f"{protocol} library chemistry"),
                      preamble(protocol, complete=True, source=source),
                      oligos(protocol, design), workflow, readout_html, "</div>"])
