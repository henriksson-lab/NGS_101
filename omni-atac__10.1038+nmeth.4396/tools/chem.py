"""Protocol-owned ATAC variant declaration and reaction graph."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nextera as nx
from protocol_families.atac import (AtacVariant, ATAC_AD1_NOMX, ATAC_AD2_1, ATAC_INDEX_1,
                                  library, workflow as build_workflow)
from pcr import PCRGeometry
from protocol_graph import Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal
from protocol_model import molecule, role

VARIANT = AtacVariant(
        key="omni-atac",
        title="Omni-ATAC — low-background chromatin accessibility sequencing",
        source='<a href="https://doi.org/10.1038/nmeth.4396">Corces et al., Nature Methods (2017)</a> and Supplementary Protocol 1.',
        summary="NP-40, digitonin and Tween permeabilize cells, a Tween wash removes mitochondria, and PBS is retained during Tn5 transposition; the resulting accessible-DNA fragments enter the original single-index ATAC-seq amplification.",
        caveat="This page follows the defining 2017 Omni-ATAC implementation: an in-house Tn5 transposase (Picelli et al. 2014, produced by an author for all experiments) in a home-made DMF-containing TD buffer, and the original Ad1/Ad2 single-index primers. The article does not print the Tn5 loading oligos; the S5/S7 mosaic-end adapter geometry is inferred from the Ad1/Ad2 primer 3' ends. Later commercial-TDE1 or dual-index implementations are not folded into this construct.",
        first_step="Lyse 50,000 cells in multi-detergent ATAC-RSB, then wash nuclei",
        first_state="washed, permeabilized nuclei depleted of free mitochondria",
        first_note="The defining recipe uses 0.1% NP-40, 0.1% Tween-20 and 0.01% digitonin for 3 minutes on ice, followed by a 0.1% Tween-20 wash before pelleting.",
        transposition_action="Transpose accessible chromatin in the PBS-containing Omni mix",
        transposition_note="The 50-microlitre reaction contains 25 microlitres home-made 2x TD buffer, 2.5 microlitres in-house Tn5 (100 nM final), 16.5 microlitres PBS, 0.01% digitonin and 0.1% Tween-20; incubate 30 minutes at 37 degrees C with 1,000-rpm mixing.",
        cleanup_note="A Zymo DNA Clean & Concentrator-5 column stops transposition. NEBNext 2x Master Mix with Ad1/Ad2 fills the nine-base Tn5 gaps (72 degrees C, 5 min) and pre-amplifies for five cycles; qPCR on a 10% aliquot determines the minimum additional cycles.",
        sequencing_intro="All deep sequencing in the defining article used 2 x 75-bp paired-end reads on a HiSeq4000 (low-depth runs: 2 x 75 bp on a MiSeq), so the two reads enter accessible genomic DNA from opposite mosaic ends. The article does not state the Index 1 cycle count; the eight-base value shown is the length of the Ad2 sample index.",
        read_lengths={"Read 1": 75, "Index 1 (i7)": 8, "Read 2": 75},
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


def oligos():
    """Displayed oligos; the Tn5 loading oligos are not printed by the defining article."""
    return (
        f"Ad1_noMX — 5′-{ATAC_AD1_NOMX}-3′",
        f"Ad2.1_TAAGGCGA (example Ad2 member) — 5′-{ATAC_AD2_1}-3′",
        f"Ad2.1 sample-sheet Index 1 value — 5′-{ATAC_INDEX_1}-3′",
        f"19-bp Tn5 mosaic end at each genomic junction (in-house Tn5 loading oligos not printed) — 5′-{nx.ME}-3′",
    )

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
    # Corces 2017 Online Methods (PMC5623106 JATS lines 307-313) cite Picelli 2014
    # (ref 26) for the transposase; Author Contributions (line 898) state that S.W.C.
    # produced the Tn5 used in all experiments.  2x TD buffer is the home-made
    # Tris/MgCl2/20% DMF recipe of Supplementary Protocol 1.
    tde1 = graph.input("tde1_adapter",
                       "S5/S7 mosaic-end adapters loaded on in-house Tn5 (loading oligos not printed; inferred from Ad1/Ad2 3' ends)",
                       MaterialKind.ADAPTER, inferred=True)
    tn5 = graph.input("tn5", "in-house Tn5 transposase (Picelli et al. 2014), 100 nM in home-made TD buffer with PBS",
                      MaterialKind.CATALYST)
    tagged = graph.product("tagged", "heterologously tagmented accessible-DNA fragments",
                           MaterialKind.MOLECULAR, state=nuclei.state)
    non_library = graph.product("non_library", "S5/S5, S7/S7 and other unusable products",
                                MaterialKind.WASTE)
    graph.react("tagment", variant.transposition_action, StepKind.TAGMENTATION,
                ((nuclei, InputRole.TEMPLATE), (tde1, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, non_library))
    graph.discard(non_library, "lacks the heterologous S5/S7 end pair required for amplification")
    cleanup = graph.input("cleanup_beads", "Zymo DNA Clean & Concentrator-5 column (21 µl elution)",
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
    polymerase = graph.input("polymerase", "NEBNext 2x Master Mix (72 °C 5 min gap fill, then 98/63/72 °C cycles)",
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
    remainder = graph.product("preamp_remainder", "remaining 45 µl pre-amplified reaction",
                              MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    qc = graph.product("qpcr_aliquot", "consumed 5-µl SYBR qPCR aliquot", MaterialKind.QC)
    graph.react("choose_cycles", "Use aliquot qPCR to choose the additional cycles",
                StepKind.SAMPLING,
                ((preamp, InputRole.TEMPLATE), (qPCR, InputRole.CONDITION)),
                (remainder, qc), changes_molecule=False)
    graph.mark_terminal(Terminal(qc.id, Disposition.QC))
    # Supplementary Protocol 1, Final amplification and cleanup: the remaining 45 ul
    # returns to the cycler with no added reagents, then a second Zymo DCC-5 column.
    extra_cycles = graph.input("extra_cycles", "qPCR-chosen additional cycles, no added reagents",
                               MaterialKind.CONDITION)
    amplified = graph.product("amplified", "final-amplified ATAC library",
                              MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    graph.react("final_amplify", "Run the qPCR-chosen additional cycles", StepKind.PCR,
                ((remainder, InputRole.TEMPLATE), (ad1, InputRole.STRUCTURAL),
                 (ad2, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST),
                 (extra_cycles, InputRole.CONDITION)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,
                    reason="Continues the pre-amplification reaction with the same Ad1/Ad2 primers; primer sites are checked on the completed library."))
    final_column = graph.input("final_column", "second Zymo DCC-5 column (20 µl water elution)",
                               MaterialKind.SELECTION)
    final = graph.product("library", "column-purified Omni-ATAC library",
                          MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    final_waste = graph.product("final_column_waste", "PCR reagents removed by the column",
                                MaterialKind.WASTE)
    graph.react("final_cleanup", "Purify the final PCR on a Zymo column", StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (final_column, InputRole.SELECTION)),
                (final, final_waste))
    graph.discard(final_waste, "is removed in the column flow-through and washes")
    # Online Methods, Sequencing (JATS lines 267-272): pooled libraries for the
    # HiSeq4000 were PAGE size-selected (>100 bp) to remove excess primers.
    page_gel = graph.input("page_gel", "PAGE gel size selection of pooled libraries (>100 bp)",
                           MaterialKind.SELECTION)
    selected = graph.product("size_selected", "PAGE-selected pooled Omni-ATAC libraries",
                             MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    primer_waste = graph.product("excess_primers", "excess primers (<100 bp)",
                                 MaterialKind.WASTE)
    graph.react("page_select", "Size-select pooled libraries on PAGE before HiSeq4000 sequencing",
                StepKind.SELECTION,
                ((final, InputRole.TEMPLATE), (page_gel, InputRole.SELECTION)),
                (selected, primer_waste))
    graph.discard(primer_waste, "is left in the gel outside the excised >100-bp region")
    graph.sequence(selected, SequencingPlan(library=FINAL_LIBRARY,
                                            primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope([
        "cells", "permeabilization_mix", "permeabilize", "tde1_adapter", "tn5",
        "tagment", "cleanup_beads", "cleanup", "ad1_primer", "ad2_primer",
        "polymerase", "pre_amplify", "qpcr_rule", "choose_cycles", "library",
        "preamp_remainder", "extra_cycles", "amplified", "final_amplify",
        "final_column", "final_column_waste", "final_cleanup", "page_gel",
        "size_selected", "excess_primers", "page_select"], source)
    return graph.require_publishable()
