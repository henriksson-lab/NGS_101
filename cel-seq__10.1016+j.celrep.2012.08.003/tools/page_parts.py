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
        supporting = (' Primer designs and bench steps: the 2012 Supplemental Information '
                      '&ldquo;CEL-Seq Protocol&rdquo; (all 40 RT primers); the same design is '
                      'restated in <a href="https://doi.org/10.1186/s13059-016-0938-8">CEL-Seq2 '
                      'Table S2</a>.')
    complete_attr = ' data-protocol-complete="true"' if complete else ""
    return f"""<div class="wrap"{complete_attr}>
<h1>{protocol} &mdash; barcoded 3' RNA-seq by T7 linear amplification</h1>
{info(f'Defining source: <a href="https://doi.org/{doi}">Hashimshony et al., '
      f'<i>{journal}</i> ({year})</a>.{supporting}')}
<div class="caveat"><b>Small-RNA library boundary.</b> The 2012 CEL-Seq protocol names
the Illumina TruSeq Small RNA kit reagents (RA3 3&prime; adapter, RTP, RP1, indexed RPIX) but
does not print their sequences (&ldquo;available from Illumina&rdquo;), and the kit revision
used in 2012 has not been read. Dotted outer arms are therefore placeholders, not a
TruSeq Small RNA reconstruction.</div>"""


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
    sp.read_spec("Read 1", "at least 15", "8-nt cell barcode, then poly(T)", "extract"),
    sp.read_spec("Index 1 (i7)", None, "library index (RPIX), read when needed", "demultiplex"),
    sp.read_spec("Read 2", 50, "transcript sequence from the 3' tag", "include"),
], "Read roles and lengths are source-declared (HiSeq 2000 paired-end: at least 15 bases "
   "for Read 1, 50 for Read 2, plus the Illumina index when needed); exact primer "
   "sequences are unavailable.")}"""


def graph_for(protocol: str, *, design: C.PrimerDesign | None = None,
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
        rt_state=C.rt_primer(protocol, design), first_state=C.rt_scene(protocol, design),
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
        protocol_src = source + " (Supplemental Information, CEL-Seq Protocol)"
        phosphatase = graph.input("phosphatase", "Antarctic phosphatase (NEB M0289), 37 °C 30 min",
                                  MaterialKind.CATALYST, source=protocol_src)
        dephosphorylated = graph.product("dephosphorylated_arna",
                                         "phosphatase-treated fragmented aRNA",
                                         MaterialKind.MOLECULAR, state=C.arna(protocol, design))
        graph.react("dephosphorylate_arna", "Dephosphorylate fragmented aRNA ends",
                    StepKind.DEPHOSPHORYLATION,
                    ((fragmented, InputRole.TEMPLATE), (phosphatase, InputRole.CATALYST)),
                    (dephosphorylated,))
        pnk = graph.input("pnk", "T4 polynucleotide kinase (NEB M0201) with ATP, 37 °C 60 min",
                          MaterialKind.CATALYST, source=protocol_src)
        kinased = graph.product("kinased_arna", "PNK-treated fragmented aRNA (RNeasy cleanup)",
                                MaterialKind.MOLECULAR, state=C.arna(protocol, design))
        graph.react("phosphorylate_arna", "Phosphorylate fragmented aRNA 5′ ends",
                    StepKind.PHOSPHORYLATION,
                    ((dephosphorylated, InputRole.TEMPLATE), (pnk, InputRole.CATALYST)),
                    (kinased,))
        reviewed.extend(("phosphatase", "dephosphorylate_arna", "pnk", "phosphorylate_arna"))
        adapter = graph.input(
            "rna_adapter", "Illumina TruSeq Small RNA kit RA3 3′ adapter, diluted 5-fold (sequence not printed)",
            MaterialKind.ADAPTER, state=C.cel1_adapter_geometry(), inferred=True,
            source=protocol_src + "; RA3 bases and terminal geometry not printed")
        ligase = graph.input("rna_ligase", "T4 RNA ligase 2, truncated (NEB M0242), 28 °C 1 h",
                             MaterialKind.CATALYST, source=protocol_src)
        ligated = graph.product("ligated_arna", "3′-adapter-ligated aRNA",
                                MaterialKind.MOLECULAR,
                                state=C.cel1_ligated_arna(), inferred=True,
                                source=protocol_src + "; adapter sequence not printed")
        graph.react("ligate_adapter", "Ligate the RNA 3′ adapter",
                    StepKind.LIGATION,
                    ((kinased, InputRole.TEMPLATE),
                     (adapter, InputRole.STRUCTURAL),
                     (ligase, InputRole.CATALYST)), (ligated,),
                    geometry=C.cel1_adapter_geometry())
        library_rt_primer = graph.input(
            "library_rt_primer", "Illumina RNA RT Primer (RTP; sequence not printed)",
            MaterialKind.OLIGO, inferred=True,
            source=protocol_src + "; RTP bases not printed")
        library_rt = graph.input("library_reverse_transcriptase",
                                 "SuperScript II reverse transcriptase, 50 °C 1 h",
                                 MaterialKind.CATALYST, source=protocol_src)
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

    if protocol == "CEL-Seq":
        fwd_label = "Illumina RNA PCR Primer RP1 (sequence not printed)"
        rev_label = "Illumina indexed RNA PCR Primer RPIX (sequence not printed)"
        pol_label = "Illumina PCR mix (PML), 12 cycles (11–15), 30-s extension"
        pol_source = source + " (Supplemental Information, CEL-Seq Protocol; Experimental Procedures)"
    else:
        fwd_label = "library-PCR forward primer (sequence not disclosed)"
        rev_label = "indexed library-PCR reverse primer (sequence not disclosed)"
        pol_label = "PCR polymerase"
        pol_source = None
    forward = graph.input(
        "pcr_forward", fwd_label,
        MaterialKind.OLIGO, inferred=True,
        source=source + " (library-PCR primer sequences not disclosed)")
    reverse = graph.input(
        "pcr_reverse", rev_label,
        MaterialKind.OLIGO, inferred=True,
        source=source + " (library-PCR primer sequences not disclosed)")
    polymerase = graph.input("pcr_polymerase", pol_label,
                             MaterialKind.CATALYST,
                             **({"source": pol_source} if pol_source else {}))
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


def protocol_graph() -> ProtocolGraph:
    """The defining CEL-Seq chemistry for source and molecular audits."""
    return graph_for("CEL-Seq")


def render_page(protocol: str, *, design: C.PrimerDesign | None = None,
                source: tuple[str, str, str] | None = None,
                graph_factory=None, readout=None) -> str:
    graph = graph_factory() if graph_factory else graph_for(protocol, design=design)
    readout_html = readout() if readout else cel1_readout()
    workflow = "<h2>Library construction</h2>\n" + protocol_graph_panels(graph)
    return "\n".join([head(f"{protocol} library chemistry"),
                      preamble(protocol, complete=True, source=source),
                      oligos(protocol, design), workflow, readout_html, "</div>"])
