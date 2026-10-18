"""Protocol-owned ATAC variant declaration and reaction graph."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nextera as nx
from types import SimpleNamespace
from dataclasses import replace
from chemdraw import Construct, Segment, feature
import seqprimers as sp
import illumina as il
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal
from protocol_model import molecule, role

VARIANT = SimpleNamespace(
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
        sequencing_intro="All deep sequencing in the defining article used 2 x 75-bp paired-end reads on a HiSeq4000 (low-depth runs: 2 x 75 bp on a MiSeq), so the two reads enter accessible genomic DNA from opposite mosaic ends. The eight-base Ad2 sample index is molecular content; its acquired Index1 cycle count remains unassigned.",
        read_lengths={"Read 1": 75, "Read 2": 75},
    )
TITLE = VARIANT.title
NOTES = "01_omni_atac.html"
SOURCE = VARIANT.source
SUMMARY = VARIANT.summary
CAVEAT = VARIANT.caveat
# Original Ad1_noMX/Ad2 geometry is an explicitly scoped predecessor
# reference for the SupplementaryProtocol1 Ad1/Ad2 reagents.
ATAC_AD1_NOMX=il.P5+nx.S5+nx.ME[:7]
ATAC_AD2_1=il.P7+"TCGCCTTA"+nx.S7+nx.ME[:6]
ATAC_INDEX_1="TAAGGCGA"


def final_library():
    lib=Construct([Segment("P5",il.P5,"p5"),Segment("S5",nx.S5,"s5"),
        Segment("left mosaic end",nx.ME,"me",inferred=True),role("accessible genomic DNA"),
        Segment("right mosaic end reverse complement",nx.ME_RC,"me",inferred=True),Segment("S7 reverse complement",nx.S7_RC,"s7"),
        Segment("i7 reverse complement","I"*8,"cbc",placeholder=True,feature=feature("sample_index_i7","sample_index","unknown")),
        Segment("P7 reverse complement",il.P7_RC,"p7")],name="Omni-ATAC original-Ad1/Ad2 reference")
    primers=(sp.NEXTERA["R1"],sp.NEXTERA["I1"],sp.NEXTERA["R2"])
    problems=sp.verify(lib,primers,required_roles=tuple(p.role for p in primers))
    if problems:raise ValueError("Omni-ATAC library: "+"; ".join(problems))
    return lib,primers


def primer_states():
    return (Construct([Segment("Ad1_noMX",ATAC_AD1_NOMX,inferred=True)]),
            Construct([Segment("P7",il.P7,inferred=True),Segment("i7 oligo","N"*8,placeholder=True,inferred=True,feature=feature("sample_index_i7","sample_index","unknown")),
                       Segment("S7 and6ME bases",nx.S7+nx.ME[:6],inferred=True)]))


def gap_filled_template():
    con=Construct([Segment("S5",nx.S5,"s5"),Segment("left mosaic end",nx.ME,"me"),
        role("accessible genomic DNA"),Segment("right mosaic end reverse complement",nx.ME_RC,"me"),
        Segment("S7 reverse complement",nx.S7_RC,"s7")],name="gap-filled heterologous transposition template")
    return Construct([replace(seg,inferred=True) for seg in con],name=con.name)

FINAL_LIBRARY,SEQ_PRIMERS=final_library()
FINAL_CAPTION = "INFERRED — complete Tn5 loading-end assignment; The original Ad1_noMX and Ad2 architecture is shown: no i5 index is present. Tn5 supplies the complete 19-bp mosaic ends, while PCR appends P5 on the S5 side and P7 plus one eight-base i7 sample index on the S7 side."
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
    permeabilization = graph.input("permeabilization_mix", "0.1%NP-40/0.1%Tween20/0.01%digitonin in ATAC-RSB;3min ice, then0.1%Tween wash",
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
    mix=graph.input("transposition_mix","50µl mix:25µl home-made2×TD,16.5µl PBS,0.01%digitonin,0.1%Tween20;37°C30min1000rpm",MaterialKind.CONDITION)
    tagged = graph.product("tagged", "INFERRED — heterologously tagmented accessible-DNA end reference",
                           MaterialKind.MOLECULAR, state=molecule("staggered transposition product",(role("inferred S5-bearing end"),role("gapped genomic DNA"),role("inferred S7-bearing end"))))
    non_library = graph.product("non_library", "S5/S5, S7/S7 and other unusable products",
                                MaterialKind.MOLECULAR)
    graph.react("tagment", variant.transposition_action, StepKind.TAGMENTATION,
                ((nuclei, InputRole.TEMPLATE), (tde1, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST),(mix,InputRole.CONDITION)), (tagged, non_library))
    graph.mark_terminal(Terminal(non_library.id,Disposition.UNRESOLVED,
        "Missing intended S5/S7 geometry does not establish failure to amplify or physical removal by the column. Class-specific recovery remains unresolved."))
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
    ad1 = graph.input("ad1_primer", "Ad1_noMX amplification primer reference", MaterialKind.OLIGO,state=primer_states()[0],inferred=True)
    ad2 = graph.input("ad2_primer", "indexed Ad2 amplification primer reference", MaterialKind.OLIGO,state=primer_states()[1],inferred=True)
    polymerase = graph.input("polymerase", "NEBNext 2x Master Mix (72 °C 5 min gap fill, then 98/63/72 °C cycles)",
                            MaterialKind.CATALYST)
    extended=graph.product("gap_filled","gap-filled accessible-DNA template",MaterialKind.MOLECULAR,state=gap_filled_template(),inferred=True)
    graph.react("gap_fill","INFERRED —72°C5min extension of transposition-end reference",StepKind.EXTENSION,
        ((purified,InputRole.TEMPLATE),(polymerase,InputRole.CATALYST)),(extended,))
    forward,reverse=("".join(x.top for x in con) for con in primer_states())
    preamp = graph.product("preamp", "five-cycle pre-amplified ATAC library",
                           MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    graph.react("pre_amplify", "INFERRED — five-cycle preamplification of the end-reference template", StepKind.PCR,
                ((extended, InputRole.TEMPLATE), (ad1, InputRole.STRUCTURAL),
                 (ad2, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (preamp,), geometry=PCRGeometry.exact(extended.state,PCRPrimer("Ad1",forward),PCRPrimer("Ad2",reverse),
                    statement="Source-qualified primer-family footprints are checked on the gap-filled precursor; this relation does not validate transposition or gap filling."))
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
                (amplified,), geometry=PCRGeometry.exact(remainder.state,PCRPrimer("Ad1",forward),PCRPrimer("Ad2",reverse),
                    statement="Continued PCR uses the actual preamplified remainder with the same primer-family architectures; physical index members remain unassigned."))
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
                                            primers=tuple(SEQ_PRIMERS),acquired_roles=("Read 1","Read 2"),read_lengths=READ_LENGTHS,
        limitation="The defining article specifies2×75 reads; the molecular Ad2 index length does not establish acquired index cycles."))
    graph.declare_reviewed_scope([
        "cells", "permeabilization_mix", "permeabilize", "tde1_adapter", "tn5",
        "transposition_mix", "tagment", "cleanup_beads", "cleanup", "ad1_primer", "ad2_primer",
        "polymerase", "gap_fill", "pre_amplify", "qpcr_rule", "choose_cycles", "library",
        "preamp_remainder", "extra_cycles", "amplified", "final_amplify",
        "final_column", "final_column_waste", "final_cleanup", "page_gel",
        "size_selected", "excess_primers", "page_select"], source)
    return graph.require_publishable()
