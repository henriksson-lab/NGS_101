"""TeloPCR-seq molecular workflow (Bennett et al. 2016)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from adapter_ligation import reported_adapter_end
from chemdraw import Construct, Row, Segment, oligo, revcomp
from dumbbell import (Dumbbell, dumbbell_rows, smrtbell_circuit,
                      smrtbell_hairpin_geometry, smrtbell_sequencing_plan)
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SourceClaim,
                            SourceDiscrepancy, StepKind)
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
# The article delegates library construction to the manufacturer's ~250 bp protocol
# (Methods, "Initial pipeline"); it does not state the end processing or hairpin
# ligating ends, so the closed topology is drawn without a junction overhang.
SMRTBELL = Dumbbell(AMPLICON, "N" * 16, "N" * 16, name="PacBio circular library",
                    insert_overhang="", adapter_overhang="", junction_unavailable=True)
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
             "Hairpin ligation closes both amplicon ends; the kit end processing and ligating ends are not stated in the article.")]


def protocol_graph():
    graph = ProtocolGraph("TeloPCR-seq", authority=SOURCE)
    chromosome = graph.input(
        "chromosome_end", "S. pombe chromosome end with native telomere terminus",
        MaterialKind.BIOLOGICAL, state=END.scene(), primary=True)
    anchor = graph.input(
        "anchor", "5′-phosphorylated, 3′-AmC3 single-stranded terminal anchor", MaterialKind.ADAPTER,
        state=ANCHOR, source=SOURCE + " (Materials and Methods)")
    circligase = graph.input(
        "circligase", "CircLigase II ssDNA ligase (Epicentre) with kit buffer and betaine; 60 °C 2 h, then 80 °C 10 min",
        MaterialKind.CATALYST)
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
    polymerase = graph.input(
        "pcr_polymerase",
        "Taq (Roche, 5% DMSO; primer 1), iProof GC buffer (primer 2) or LongAmp Taq 2X; 45 cycles, 2 µL ligation per 50 µL",
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
    pool_reactions = graph.input(
        "pcr_replicates", "five further 50 µL TeloPCR reactions of the same sample",
        MaterialKind.MOLECULAR, state=AMPLICON, inferred=True,
        source=SOURCE + " (Methods: 6 tubes of 50 µL PCR products pooled)")
    pooled = graph.product("pooled", "pooled TeloPCR products (6 × 50 µL)",
                           MaterialKind.MOLECULAR, state=AMPLICON, inferred=True)
    graph.react("pool_pcr", "Pool six 50 µL PCR reactions", StepKind.POOL,
                ((amplicon, InputRole.TEMPLATE), (pool_reactions, InputRole.TEMPLATE)),
                (pooled,), changes_molecule=False)
    filters = graph.input(
        "primer_removal", "Amicon Ultra 0.5 mL centrifugal filter, then Qiagen PCR purification kit",
        MaterialKind.SELECTION)
    purified = graph.product("purified", "purified telomere PCR products",
                             MaterialKind.MOLECULAR, state=AMPLICON, inferred=True)
    free_primers = graph.product("free_primers", "free primers, filtrate and column washes",
                                 MaterialKind.WASTE)
    graph.react("remove_primers", "Concentrate, buffer-exchange and column-purify to remove free primers",
                StepKind.CLEANUP,
                ((pooled, InputRole.TEMPLATE), (filters, InputRole.SELECTION)),
                (purified, free_primers), changes_molecule=False)
    graph.discard(free_primers, "free primers removed by centrifugal filtration and column purification")
    kit = graph.input(
        "pacbio_250bp_protocol", "Pacific Biosciences ~250 bp library protocol (Yale Center for Genome Analysis)",
        MaterialKind.CONDITION)
    prepared = graph.product(
        "end_prepared", "amplicon with kit-processed ends (end chemistry not stated)",
        MaterialKind.MOLECULAR, state=AMPLICON, inferred=True)
    graph.react("end_prep", "Process amplicon ends per the manufacturer's ~250 bp protocol",
                StepKind.MAGIC, ((purified, InputRole.TEMPLATE), (kit, InputRole.CONDITION)),
                (prepared,),
                reason="The article says only that libraries were generated per the manufacturer's ~250 bp protocol; the historical kit revision, end repair and any tailing are not stated.")
    hairpins = graph.input(
        "hairpins", "PacBio SMRTbell hairpin adapters (kit-supplied)", MaterialKind.ADAPTER,
        inferred=True, source=SOURCE + " (hairpin-body bases proprietary)")
    smrt_ligase = graph.input("smrt_ligase", "SMRTbell ligase", MaterialKind.CATALYST)
    circuit = smrtbell_circuit(SMRTBELL)
    library = graph.product(
        "library", SMRTBELL.name, MaterialKind.LIBRARY, state=circuit,
        inferred=True, source=SOURCE + " (PacBio hairpin-body bases unavailable)")
    graph.react(
        "ligate_hairpins", "Ligate hairpins to both amplicon ends", StepKind.LIGATION,
        ((prepared, InputRole.TEMPLATE), (hairpins, InputRole.STRUCTURAL),
         (smrt_ligase, InputRole.CATALYST)), (library,),
        geometry=smrtbell_hairpin_geometry(
            SMRTBELL,
            statement="The article names PacBio circular libraries from the manufacturer's ~250 bp protocol but not the hairpin ligating ends or end processing."))
    graph.sequence(library, smrtbell_sequencing_plan(
        SMRTBELL, circuit,
        limitation="The paper specifies PacBio RS circular consensus sequencing but does not print the proprietary hairpin primer-site bases; reads lacking both amplifying primers (≤2 mismatches) or any GGTTAC were excluded computationally."))
    graph.annotate(SourceDiscrepancy(
        "remove_primers", "Free-primer removal device named differently",
        (SourceClaim("Materials and Methods, Telo-PCR (PMC5561429 JATS line 48)",
                     "Amicon Ultra 0.5 mL centrifugal filter (Millipore), then Qiagen PCR purification kit"),
         SourceClaim("Figure 1A legend (PMC5561429 JATS line 374)",
                     "purified to remove free primers using a Microcon concentrator")),
        "The model follows the Methods: Amicon Ultra 0.5 mL filtration followed by column purification.",
        "Both statements are in the published FEBS Letters article (author manuscript, PMC5561429).",
        "Both are Millipore centrifugal concentrators; the article does not state the filter cut-off, so which device or whether both were used is unresolved."))
    graph.declare_reviewed_scope((
        "chromosome_end", "anchor", "circligase", "attach_anchor",
        "subtelomere_primer", "anchor_primer", "pcr_polymerase", "telomere_pcr",
        "pcr_replicates", "pool_pcr", "primer_removal", "remove_primers",
        "pacbio_250bp_protocol", "end_prep", "hairpins", "smrt_ligase", "ligate_hairpins",
        "library"), SOURCE)
    return graph.require_publishable()
