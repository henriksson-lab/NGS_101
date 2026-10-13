"""TeloPCR-seq molecular workflow (Bennett et al. 2016)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Row, Segment
from dumbbell import (Dumbbell, dumbbell_rows, smrtbell_circuit,
                      smrtbell_hairpin_geometry, smrtbell_sequencing_plan)
from endprep import dA_tailed_scene, repair_and_dA_tail
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, StepKind)
from telomere import TelomereEnd

END = TelomereEnd(repeat="GGTTACA", duplex_copies=3, overhang_copies=2)
ANCHOR = "N" * 24  # structure published; exact ordered bases are in the inaccessible supplement
AMPLICON = Construct([Segment("subtelomere primer flank", "X" * 18, placeholder=True),
                      Segment("telomere", END.repeat * 4, "r1"),
                      Segment("terminal anchor", ANCHOR, "r3", placeholder=True)], name="TeloPCR amplicon")
SMRTBELL = Dumbbell(AMPLICON, "N" * 16, "N" * 16, name="PacBio circular library")
SOURCE = "Bennett et al., FEBS Letters (2016), doi:10.1002/1873-3468.12444"

def panels():
    return [("Attach a single-stranded anchor to the native chromosome terminus", END.scene(),
             "CircLigase joins the anchor to the native telomere end; the exact supplement oligo is not fabricated here."),
            ("Amplify from a subtelomere-specific primer to the anchor primer",
             [Row(chunks=[("subtelomere primer  --->  [subtelomere][telomere repeats][terminal anchor]  <--- anchor primer", None, False)])],
             "Each chromosome end is amplified with its subtelomeric primer and the shared anchor primer."),
            ("Convert the amplicon into a circular PacBio template", dumbbell_rows(SMRTBELL),
             "Hairpin ligation closes both amplicon ends; ** marks the sealed insert–hairpin junctions.")]


def protocol_graph():
    graph = ProtocolGraph("TeloPCR-seq", authority=SOURCE)
    chromosome = graph.input(
        "chromosome_end", "S. pombe chromosome end with native telomere terminus",
        MaterialKind.BIOLOGICAL, state=END.scene(), primary=True)
    anchor = graph.input(
        "anchor", "single-stranded terminal anchor", MaterialKind.ADAPTER,
        inferred=True, source=SOURCE + " (ordered anchor bases unavailable)")
    circligase = graph.input("circligase", "CircLigase", MaterialKind.CATALYST)
    anchored = graph.product(
        "anchored_end", "anchor-tagged native chromosome terminus",
        MaterialKind.MOLECULAR, state=END.scene(), inferred=True,
        source=SOURCE + " (ordered anchor bases unavailable)")
    graph.react(
        "attach_anchor", "Attach the single-stranded anchor to the native terminus",
        StepKind.LIGATION,
        ((chromosome, InputRole.TEMPLATE), (anchor, InputRole.STRUCTURAL),
         (circligase, InputRole.CATALYST)), (anchored,),
        geometry=AdapterLigationGeometry.unavailable(
            "The source establishes CircLigase-mediated terminal anchor attachment, "
            "but the inaccessible supplement contains the ordered anchor bases."))
    subtelomere_primer = graph.input(
        "subtelomere_primer", "chromosome-end-specific subtelomere PCR primer",
        MaterialKind.OLIGO, inferred=True,
        source=SOURCE + " (ordered primer bases unavailable)")
    anchor_primer = graph.input(
        "anchor_primer", "shared terminal-anchor PCR primer", MaterialKind.OLIGO,
        inferred=True, source=SOURCE + " (ordered primer bases unavailable)")
    polymerase = graph.input("pcr_polymerase", "telomere PCR polymerase",
                             MaterialKind.CATALYST)
    amplicon = graph.product(
        "amplicon", AMPLICON.name, MaterialKind.MOLECULAR, state=AMPLICON,
        inferred=True, source=SOURCE + " (ordered primer and anchor bases unavailable)")
    graph.react(
        "telomere_pcr", "Amplify from the subtelomere to the terminal anchor",
        StepKind.PCR,
        ((anchored, InputRole.TEMPLATE), (subtelomere_primer, InputRole.STRUCTURAL),
         (anchor_primer, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
        (amplicon,), geometry=PCRGeometry.unavailable(
            AMPLICON, reason="The source-defined primer roles are known, but their order-sheet bases are unavailable."))
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
        inferred=True, source=SOURCE + " (anchor and hairpin-body bases unavailable)")
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
