"""Protocol-owned ATAC variant declaration and reaction graph."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from protocol_families.atac import AtacVariant, workflow as build_workflow
from batch_ngs import nextera_library, seg
import illumina as il
import nextera as nx
from pcr import PCRGeometry
from protocol_graph import Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal
from protocol_model import molecule, role

VARIANT = AtacVariant(
        key="fast-atac",
        title="Fast-ATAC — one-step accessibility sequencing for primary blood cells",
        source='<a href="https://doi.org/10.1038/ng.3646">Corces et al., Nature Genetics (2016)</a> and its Online Methods.',
        summary="Digitonin permeabilizes 5,000 sorted primary blood cells inside the Tn5 reaction itself. The defining study then used dual-index ATAC primers and a NextSeq 76×8×8×76 run.",
        caveat="The article cites the 2015 single-cell ATAC method for modified primers; its supplement prints both i5 and i7 primer families. The particular pair used for each Fast-ATAC sample is not identified here, so both index sequences are placeholders.",
        first_step="Sort 5,000 primary blood cells and pellet them cold",
        first_state="small intact-cell pellet",
        first_note="The defining workflow sorts purified blood populations, pellets them at 500 r.c.f. for 5 minutes at 4 degrees C, and removes the supernatant without a separate detergent lysis.",
        transposition_action="Permeabilize and transpose together with digitonin",
        transposition_note="The 50-microlitre reaction contains 25 microlitres 2x TD buffer, 2.5 microlitres TDE1, 0.5 microlitre 1% digitonin and 22 microlitres water; incubate 30 minutes at 37 degrees C with 300-rpm mixing.",
        cleanup_note="MinElute cleanup stops the one-tube reaction. Five initial PCR cycles fill the Tn5 gaps and install the complete library arms; an aliquot qPCR chooses the remaining cycles.",
        sequencing_intro="The defining article's Online Methods specify a NextSeq paired-end dual-index 76×8×8×76 run. Read 1 and Read 2 enter accessible DNA from opposite mosaic ends; the index reads identify the selected Ad2 i7 and Ad1 i5 primers.",
        read_lengths={"Read 1": 76, "Index 1 (i7)": 8, "Index 2 (i5)": 8, "Read 2": 76},
        dual_index=True,
    )
TITLE = VARIANT.title
NOTES = "01_fast_atac.html"
SOURCE = VARIANT.source
SUMMARY = VARIANT.summary
CAVEAT = VARIANT.caveat
workflow = lambda: build_workflow(VARIANT)
FINAL_LIBRARY, SEQ_PRIMERS = nextera_library(
    [seg("accessible genomic DNA", "X" * 44, placeholder=True)],
    "Fast-ATAC dual-index library", dual_index=True, i5_nt=8, i7_nt=8)
FINAL_CAPTION = "The 2015 dual-index Ad1/Ad2 family cited by the defining article is shown. Tn5 supplies complete 19-bp mosaic ends; PCR appends P5 and an eight-base i5 on the S5 side, and P7 plus an eight-base i7 on the S7 side. Exact sample index members are not assigned."
READ_LENGTHS = VARIANT.read_lengths
SEQUENCING_INTRO = VARIANT.sequencing_intro
SEQUENCING_ENDING = ""
SEQUENCING_UNAVAILABLE = ""

def oligos():
    ad1_1 = il.P5 + "TAGATCGC" + nx.S5 + nx.ME[:10]
    ad2_1 = il.P7 + "TCGCCTTA" + nx.S7 + nx.ME[:7]
    return (
        f"v2_Ad1.1_TAGATCGC — 5′-{ad1_1}-3′ (illustrative i5 member; 2015 supplement)",
        f"v2_Ad2.1_TAAGGCGA — 5′-{ad2_1}-3′ (illustrative i7 member; 2015 supplement)",
        f"TDE1 mosaic end at each junction — 5′-{nx.ME}-3′",
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
    ad1 = graph.input("ad1_primer", "indexed v2_Ad1 amplification primer family", MaterialKind.OLIGO)
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
