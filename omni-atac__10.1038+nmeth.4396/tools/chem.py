"""Protocol-owned ATAC variant declaration and reaction graph."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from protocol_families.atac import AtacVariant, library, workflow as build_workflow, oligos
from pcr import PCRGeometry
from protocol_graph import Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal
from protocol_model import molecule, role

VARIANT = AtacVariant(
        key="omni-atac",
        title="Omni-ATAC — low-background chromatin accessibility sequencing",
        source='<a href="https://doi.org/10.1038/nmeth.4396">Corces et al., Nature Methods (2017)</a> and Supplementary Protocol 1.',
        summary="NP-40, digitonin and Tween permeabilize cells, a Tween wash removes mitochondria, and PBS is retained during Tn5 transposition; the resulting accessible-DNA fragments enter the original single-index ATAC-seq amplification.",
        caveat="This page follows the defining 2017 Omni-ATAC implementation with Illumina TDE1 and the original Ad1/Ad2 single-index primers. Later protocols often substitute home-made Tn5 or dual-index primers; those substitutions are not silently folded into this construct.",
        first_step="Lyse 50,000 cells in multi-detergent ATAC-RSB, then wash nuclei",
        first_state="washed, permeabilized nuclei depleted of free mitochondria",
        first_note="The defining recipe uses 0.1% NP-40, 0.1% Tween-20 and 0.01% digitonin for 3 minutes on ice, followed by a 0.1% Tween-20 wash before pelleting.",
        transposition_action="Transpose accessible chromatin in the PBS-containing Omni mix",
        transposition_note="The 50-microlitre reaction contains 25 microlitres 2x TD buffer, 2.5 microlitres TDE1, 16.5 microlitres PBS, 0.01% digitonin and 0.1% Tween-20; incubate 30 minutes at 37 degrees C with 1,000-rpm mixing.",
        cleanup_note="A silica cleanup stops transposition. Five initial PCR cycles fill the nine-base Tn5 gaps and install flow-cell/index sequence; qPCR on an aliquot determines the minimum additional cycles.",
        sequencing_intro="The defining method uses paired-end Illumina sequencing so the two reads enter accessible genomic DNA from opposite mosaic ends. It does not prescribe one immutable genomic-read length; the single eight-base Index 1 read reports the Ad2 sample index.",
        read_lengths={"Index 1 (i7)": 8},
    )
TITLE = VARIANT.title
NOTES = "01_omni_atac.html"
SOURCE = VARIANT.source
SUMMARY = VARIANT.summary
CAVEAT = VARIANT.caveat
workflow = lambda: build_workflow(VARIANT)
FINAL_LIBRARY, SEQ_PRIMERS = library(TITLE.split(" — ", 1)[0])
FINAL_CAPTION = "The original Ad1_noMX and Ad2 architecture is shown: no i5 index is present. Tn5 supplies the complete 19-bp mosaic ends, while PCR appends P5 on the S5 side and P7 plus one eight-base i7 sample index on the S7 side."
READ_LENGTHS = VARIANT.read_lengths
SEQUENCING_INTRO = VARIANT.sequencing_intro
SEQUENCING_ENDING = ""
SEQUENCING_UNAVAILABLE = ""

def protocol_graph():
    variant = VARIANT
    source = SOURCE
    graph = ProtocolGraph(TITLE, authority=source)
    cells = graph.input("cells", "fresh input cells", MaterialKind.BIOLOGICAL,
                        state=molecule("fresh cells", (role("nuclear and mitochondrial DNA"),)),
                        primary=True)
    permeabilization = graph.input("permeabilization_mix", variant.first_step,
                                   MaterialKind.CONDITION)
    nuclei = graph.product("nuclei", variant.first_state, MaterialKind.MOLECULAR,
                           state=molecule(variant.first_state, (role("accessible chromatin"),)))
    graph.react("permeabilize", variant.first_step, StepKind.PERMEABILIZATION,
                ((cells, InputRole.TEMPLATE),
                 (permeabilization, InputRole.CONDITION)), (nuclei,))
    tde1 = graph.input("tde1_adapter", "TDE1 S5/S7 mosaic-end adapter mixture",
                       MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "TDE1 Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "heterologously tagmented accessible-DNA fragments",
                           MaterialKind.MOLECULAR, state=nuclei.state)
    non_library = graph.product("non_library", "S5/S5, S7/S7 and other unusable products",
                                MaterialKind.WASTE)
    graph.react("tagment", variant.transposition_action, StepKind.TAGMENTATION,
                ((nuclei, InputRole.TEMPLATE), (tde1, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, non_library))
    graph.discard(non_library, "lacks the heterologous S5/S7 end pair required for amplification")
    cleanup = graph.input("cleanup_beads", "post-transposition cleanup matrix",
                          MaterialKind.SELECTION)
    purified = graph.product("purified", "purified tagmented accessible-DNA fragments",
                             MaterialKind.MOLECULAR, state=tagged.state)
    cleanup_waste = graph.product("cleanup_waste", "removed transposition reagents and short products",
                                  MaterialKind.WASTE)
    graph.react("cleanup", "Stop transposition and purify tagged DNA", StepKind.CLEANUP,
                ((tagged, InputRole.TEMPLATE), (cleanup, InputRole.SELECTION)),
                (purified, cleanup_waste))
    graph.discard(cleanup_waste, "is removed during the post-transposition cleanup")
    ad1 = graph.input("ad1_primer", "Ad1_noMX amplification primer", MaterialKind.OLIGO)
    ad2 = graph.input("ad2_primer", "indexed Ad2 amplification primer", MaterialKind.OLIGO)
    polymerase = graph.input("polymerase", "ATAC gap-fill/PCR polymerase",
                            MaterialKind.CATALYST)
    preamp = graph.product("preamp", "five-cycle pre-amplified ATAC library",
                           MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    graph.react("pre_amplify", "Gap-fill and pre-amplify for five cycles", StepKind.PCR,
                ((purified, InputRole.TEMPLATE), (ad1, InputRole.STRUCTURAL),
                 (ad2, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (preamp,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,
                    reason="The exact final primer sites are checked on the completed library; the staggered pre-gap-fill substrate is not represented as a Construct."))
    qPCR = graph.input("qpcr_rule", "aliquot-qPCR cycle-selection rule",
                       MaterialKind.CONDITION)
    final = graph.product("library", "complexity-preserving ATAC sequencing library",
                          MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    qc = graph.product("qpcr_aliquot", "consumed qPCR aliquot", MaterialKind.QC)
    graph.react("choose_cycles", "Use aliquot qPCR to choose the additional cycles",
                StepKind.SAMPLING,
                ((preamp, InputRole.TEMPLATE), (qPCR, InputRole.CONDITION)),
                (final, qc), changes_molecule=False)
    graph.mark_terminal(Terminal(qc.id, Disposition.QC))
    graph.sequence(final, SequencingPlan(library=FINAL_LIBRARY,
                                         primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope([
        "cells", "permeabilization_mix", "permeabilize", "tde1_adapter", "tn5",
        "tagment", "cleanup_beads", "cleanup", "ad1_primer", "ad2_primer",
        "polymerase", "pre_amplify", "qpcr_rule", "choose_cycles", "library"], source)
    return graph.require_publishable()
