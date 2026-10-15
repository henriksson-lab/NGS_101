"""TeloPCR-seq molecular workflow (Bennett et al. 2016)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from adapter_ligation import reported_adapter_end
from chemdraw import Construct, Row, Segment, oligo, revcomp
from dumbbell import (Dumbbell, dumbbell_rows, smrtbell_circuit,
                      smrtbell_hairpin_geometry, smrtbell_sequencing_plan)
from endprep import dA_tailed_scene, repair_and_dA_tail
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, StepKind)
from telomere import TelomereEnd

END = TelomereEnd(repeat="GGTTACA", duplex_copies=3, overhang_copies=2)
# Bennett et al., Materials and Methods, "Telomere PCR ... and high-throughput sequencing".
ANCHOR = "AGGGTTAATAAGCGGCCGCGTCG"  # 5′ phosphate; 3′ AmC3 blocker
# The Methods name this primer as "the reverse complement of the anchor sequence"
# and print the same bases, so derive it instead of retyping it.
ANCHOR_PCR_PRIMER = revcomp(ANCHOR)
SUBTELOMERIC_PRIMER_1 = "CGAGGCTGCGGGTTACAAGGTTACG"
SUBTELOMERIC_PRIMER_2 = "GGGAATTTAGGAAGTGCGGTAAGTTG"
AMPLICON = Construct([Segment("subtelomere primer 1 flank", SUBTELOMERIC_PRIMER_1),
                      Segment("telomere", END.repeat * 4, "r1"),
                      Segment("terminal anchor", ANCHOR, "r3")], name="TeloPCR amplicon")
SMRTBELL = Dumbbell(AMPLICON, "N" * 16, "N" * 16, name="PacBio circular library")
SOURCE = "Bennett et al., FEBS Letters (2016), doi:10.1002/1873-3468.12444"

def oligos():
    return (
        oligo("5′-phosphorylated terminal anchor, 3′ AmC3 blocked",
              (Segment("anchor", ANCHOR),), five="5′-[Phos]", three="[AmC3]-3′"),
        oligo("Subtelomeric PCR primer 1", (Segment("primer 1", SUBTELOMERIC_PRIMER_1),)),
        oligo("Subtelomeric PCR primer 2", (Segment("primer 2", SUBTELOMERIC_PRIMER_2),)),
        oligo("Anchor PCR primer", (Segment("anchor PCR primer", ANCHOR_PCR_PRIMER),)),
    )

def panels():
    return [("Attach a single-stranded anchor to the native chromosome terminus", END.scene(),
             "CircLigase joins the 5′-phosphorylated anchor to the native 3′ telomere end; the anchor bears a 3′ AmC3 blocker."),
            ("Amplify from a subtelomere-specific primer to the anchor primer",
             [Row(chunks=[("subtelomere primer  --->  [subtelomere][telomere repeats][terminal anchor]  <--- anchor primer", None, False)])],
             "The paper prints two subtelomeric primers and the shared reverse-complement anchor primer."),
            ("Convert the amplicon into a circular PacBio template", dumbbell_rows(SMRTBELL),
             "Hairpin ligation closes both amplicon ends; ** marks the sealed insert–hairpin junctions.")]


def protocol_graph():
    graph = ProtocolGraph("TeloPCR-seq", authority=SOURCE)
    chromosome = graph.input(
        "chromosome_end", "S. pombe chromosome end with native telomere terminus",
        MaterialKind.BIOLOGICAL, state=END.scene(), primary=True)
    anchor = graph.input(
        "anchor", "5′-phosphorylated, 3′-AmC3 single-stranded terminal anchor", MaterialKind.ADAPTER,
        state=ANCHOR, source=SOURCE + " (Materials and Methods)")
    circligase = graph.input("circligase", "CircLigase", MaterialKind.CATALYST)
    anchored = graph.product(
        "anchored_end", "anchor-tagged native chromosome terminus",
        MaterialKind.MOLECULAR, state=END.scene(), inferred=True,
        source=SOURCE + " (anchor junction not explicitly sequenced)")
    graph.react(
        "attach_anchor", "Attach the single-stranded anchor to the native terminus",
        StepKind.LIGATION,
        ((chromosome, InputRole.TEMPLATE), (anchor, InputRole.STRUCTURAL),
         (circligase, InputRole.CATALYST)), (anchored,),
        geometry=reported_adapter_end(
            "TeloPCR anchor", kind="single-strand", polarity="5-prime",
            terminal_sequence=ANCHOR,
            statement="The Methods give a 5′ phosphate and a 3′ AmC3 blocker; CircLigase joins the anchor 5′ end to the native telomere 3′ end."))
    subtelomere_primer = graph.input(
        "subtelomere_primer", "chromosome-end-specific subtelomere PCR primer",
        MaterialKind.OLIGO, state=(SUBTELOMERIC_PRIMER_1, SUBTELOMERIC_PRIMER_2),
        source=SOURCE + " (Materials and Methods)")
    anchor_primer = graph.input(
        "anchor_primer", "shared terminal-anchor PCR primer", MaterialKind.OLIGO,
        state=ANCHOR_PCR_PRIMER, source=SOURCE + " (Materials and Methods)")
    polymerase = graph.input("pcr_polymerase", "telomere PCR polymerase",
                             MaterialKind.CATALYST)
    amplicon = graph.product(
        "amplicon", AMPLICON.name, MaterialKind.MOLECULAR, state=AMPLICON,
        inferred=True, source=SOURCE + " (representative subtelomere and repeat count shown)")
    graph.react(
        "telomere_pcr", "Amplify from the subtelomere to the terminal anchor",
        StepKind.PCR,
        ((anchored, InputRole.TEMPLATE), (subtelomere_primer, InputRole.STRUCTURAL),
         (anchor_primer, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
        (amplicon,), geometry=PCRGeometry.unavailable(
            AMPLICON, reason="The paper prints both PCR primers, but does not give a complete source-defined genomic template for exact binding coordinates."))
    end_mix = graph.input(
        "end_mix", "PacBio end-repair and dA-tailing mix", MaterialKind.CATALYST)
    end_product = repair_and_dA_tail(AMPLICON, inferred=True)
    tailed = graph.product(
        "tailed", "end-repaired dA-tailed TeloPCR amplicon",
        MaterialKind.MOLECULAR, state=dA_tailed_scene(end_product))
    graph.react("end_prep", "Repair and dA-tail the telomere amplicon",
                StepKind.END_REPAIR,
                ((amplicon, InputRole.TEMPLATE), (end_mix, InputRole.CATALYST)),
                (tailed,))
    hairpins = graph.input(
        "hairpins", "T-overhang PacBio hairpin adapters", MaterialKind.ADAPTER,
        inferred=True, source=SOURCE + " (hairpin-body bases proprietary)")
    smrt_ligase = graph.input("smrt_ligase", "SMRTbell ligase", MaterialKind.CATALYST)
    circuit = smrtbell_circuit(SMRTBELL)
    library = graph.product(
        "library", SMRTBELL.name, MaterialKind.LIBRARY, state=circuit,
        inferred=True, source=SOURCE + " (PacBio hairpin-body bases unavailable)")
    graph.react(
        "ligate_hairpins", "Ligate hairpins to both amplicon ends", StepKind.LIGATION,
        ((tailed, InputRole.TEMPLATE), (hairpins, InputRole.STRUCTURAL),
         (smrt_ligase, InputRole.CATALYST)), (library,),
        geometry=smrtbell_hairpin_geometry(
            SMRTBELL,
            statement="PacBio library preparation ligates T-overhang hairpins to dA-tailed amplicons."))
    graph.sequence(library, smrtbell_sequencing_plan(
        SMRTBELL, circuit,
        limitation="The paper specifies PacBio circular sequencing but does not print the proprietary hairpin primer-site bases."))
    graph.declare_reviewed_scope((
        "chromosome_end", "anchor", "circligase", "attach_anchor",
        "subtelomere_primer", "anchor_primer", "pcr_polymerase", "telomere_pcr",
        "end_mix", "end_prep", "hairpins", "smrt_ligase", "ligate_hairpins",
        "library"), SOURCE)
    return graph.require_publishable()
