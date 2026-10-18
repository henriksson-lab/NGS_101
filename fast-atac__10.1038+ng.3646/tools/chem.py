"""Protocol-owned ATAC variant declaration and reaction graph."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from types import SimpleNamespace
from chemdraw import Construct, Segment, feature
import seqprimers as sp
import illumina as il
import nextera as nx
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal
from protocol_model import molecule, role

SOURCE_OLIGO_EXCEPTIONS={
    "ad1_primer":"2015 Nature14590 SupplementaryTable1 MOESM36 predecessor original lives in its own bundle; local suite checks every published v2_Ad1 family member.",
    "ad2_primer":"2015 Nature14590 SupplementaryTable1 MOESM36 predecessor original lives in its own bundle; local suite checks every published v2_Ad2 family member."
}

VARIANT = SimpleNamespace(
        key="fast-atac",
        title="Fast-ATAC — one-step accessibility sequencing for primary blood cells",
        source='<a href="https://doi.org/10.1038/ng.3646">Corces et al., Nature Genetics (2016)</a> and its Online Methods.',
        summary="Digitonin permeabilizes 5,000 sorted primary blood cells inside the Tn5 reaction itself. The defining study then used dual-index ATAC primers and a NextSeq 76×8×8×76 run.",
        caveat="The article states only the sort, pellet, digitonin/Tn5 transposition and MinElute steps itself; amplification is cited to the 2015 Current Protocols ATAC-seq workflow (gap-fill, five cycles, qPCR-chosen extra cycles, NEBNext master mix), which this page shows as an inherited procedure. Modified primers are cited to the 2015 single-cell ATAC study, whose supplement prints both i5 and i7 primer families. The particular pair used for each Fast-ATAC sample is not identified here, so both index sequences are placeholders.",
        first_step="Sort 5,000 primary blood cells and pellet them cold",
        first_state="small intact-cell pellet",
        first_note="The defining workflow sorts purified blood populations, pellets them at 500 r.c.f. for 5 minutes at 4 degrees C, and removes the supernatant without a separate detergent lysis.",
        transposition_action="Permeabilize and transpose together with digitonin",
        transposition_note="The 50-microlitre reaction contains 25 microlitres 2x TD buffer, 2.5 microlitres TDE1, 0.5 microlitre 1% digitonin and 22 microlitres water; incubate 30 minutes at 37 degrees C with 300-rpm mixing.",
        cleanup_note="MinElute Reaction Cleanup stops the one-tube reaction. The article then amplifies \"as described previously\" (Buenrostro et al. 2015, Curr Protoc Mol Biol): a 72 degrees C gap-fill, five initial PCR cycles, and an aliquot qPCR that chooses the remaining cycles.",
        sequencing_intro="The defining article's Online Methods specify a NextSeq paired-end dual-index 76×8×8×76 run. Read 1 and Read 2 enter accessible DNA from opposite mosaic ends; the index reads identify the selected Ad2 i7 and Ad1 i5 primers.",
        read_lengths={"Read 1": 76, "Index 1 (i7)": 8, "Index 2 (i5)": 8, "Read 2": 76},
        dual_index=True,
    )
TITLE = VARIANT.title
NOTES = "01_fast_atac.html"
SOURCE = VARIANT.source
SUMMARY = VARIANT.summary
CAVEAT = VARIANT.caveat
def final_library():
    # Fast-ATAC cites the2015 single-cell Ad1/Ad2 dual-index primer family;
    # both8-base member assignments remain unknown for the chosen sample.
    lib=Construct([Segment("P5",il.P5,"p5"),
        Segment("i5","I"*8,"cbc",placeholder=True,feature=feature("sample_index_i5","sample_index","unknown")),
        Segment("S5",nx.S5,"s5"),Segment("left mosaic end",nx.ME,"me"),
        role("accessible genomic DNA"),Segment("right mosaic end reverse complement",nx.ME_RC,"me"),
        Segment("S7 reverse complement",nx.S7_RC,"s7"),
        Segment("i7 reverse complement","I"*8,"cbc",placeholder=True,feature=feature("sample_index_i7","sample_index","unknown")),
        Segment("P7 reverse complement",il.P7_RC,"p7")],name="Fast-ATAC2015-family dual-index library")
    primers=(sp.NEXTERA["R1"],sp.NEXTERA["I1"],sp.NEXTERA["I2"],sp.NEXTERA["R2"])
    problems=sp.verify(lib,primers,required_roles=tuple(p.role for p in primers))
    if problems:raise ValueError("Fast-ATAC library: "+"; ".join(problems))
    return lib,primers


def primer_states():
    return (Construct([Segment("P5",il.P5),Segment("i5","N"*8,placeholder=True,feature=feature("sample_index_i5","sample_index","unknown")),
                       Segment("S5 and10ME bases",nx.S5+nx.ME[:10])]),
            Construct([Segment("P7",il.P7),Segment("i7 oligo","N"*8,placeholder=True,feature=feature("sample_index_i7","sample_index","unknown")),
                       Segment("S7 and7ME bases",nx.S7+nx.ME[:7])]))


def gap_filled_template():
    return Construct([Segment("S5",nx.S5,"s5"),Segment("left mosaic end",nx.ME,"me"),
        role("accessible genomic DNA"),Segment("right mosaic end reverse complement",nx.ME_RC,"me"),
        Segment("S7 reverse complement",nx.S7_RC,"s7")],name="gap-filled heterologous transposition template")

FINAL_LIBRARY,SEQ_PRIMERS=final_library()
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
    cells = graph.input("cells", "FACS-sorted primary blood cells in FACS buffer",
                        MaterialKind.BIOLOGICAL,
                        state=molecule("sorted cells", (role("nuclear and mitochondrial DNA"),)),
                        primary=True)
    pelleting = graph.input("pellet_condition", "500 r.c.f., 5 min, 4 °C fixed-angle spin",
                            MaterialKind.CONDITION)
    pellet = graph.product("pellet", variant.first_state, MaterialKind.MOLECULAR,
                           state=molecule(variant.first_state, (role("intact cells with chromatin"),)))
    supernatant = graph.product("supernatant", "removed FACS-buffer supernatant",
                                MaterialKind.WASTE)
    # Corces 2016 Online Methods: 5,000 sorted cells are pelleted and all supernatant
    # is removed; there is no separate lysis step before the transposition mix.
    graph.react("pellet_cells", variant.first_step, StepKind.CLEANUP,
                ((cells, InputRole.TEMPLATE),
                 (pelleting, InputRole.CONDITION)), (pellet, supernatant))
    graph.discard(supernatant, "is removed in two pipetting steps before transposition")
    tde1 = graph.input("tde1_adapter", "TDE1 S5/S7 mosaic-end adapter mixture",
                       MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "TDE1 Tn5 transposase (Illumina FC-121-1030)", MaterialKind.CATALYST)
    digitonin = graph.input("digitonin_td", "2x TD buffer plus 0.5 µl 1% digitonin (Promega G9441) in 50 µl",
                            MaterialKind.CONDITION)
    tagged = graph.product("tagged", "heterologously tagmented accessible-DNA fragments",
                           MaterialKind.MOLECULAR, state=molecule("staggered transposition product",(role("S5-bearing end"),role("gapped genomic DNA"),role("S7-bearing end"))))
    non_library = graph.product("non_library", "S5/S5, S7/S7 and other unusable products",
                                MaterialKind.MOLECULAR)
    graph.react("tagment", variant.transposition_action, StepKind.TAGMENTATION,
                ((pellet, InputRole.TEMPLATE), (tde1, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST), (digitonin, InputRole.CONDITION)),
                (tagged, non_library))
    graph.mark_terminal(Terminal(non_library.id,Disposition.UNRESOLVED,
        "Missing intended S5/S7 geometry does not establish failure to amplify or physical removal by the column. Class-specific recovery remains unresolved."))
    cleanup = graph.input("cleanup_column", "QIAGEN MinElute Reaction Cleanup column (Cat# 28204)",
                          MaterialKind.SELECTION)
    purified = graph.product("purified", "purified tagmented DNA in 10 µl elution buffer",
                             MaterialKind.MOLECULAR, state=tagged.state)
    cleanup_waste = graph.product("cleanup_waste", "removed transposition reagents, protein and detergent",
                                  MaterialKind.WASTE)
    graph.react("cleanup", "Stop transposition and purify tagged DNA", StepKind.CLEANUP,
                ((tagged, InputRole.TEMPLATE), (cleanup, InputRole.SELECTION)),
                (purified, cleanup_waste))
    graph.discard(cleanup_waste, "is removed during the post-transposition column cleanup")
    # Amplification below is inherited: the defining article says only "amplified and
    # purified as described previously" (ref. 48, Buenrostro 2015 Curr Protoc) "with
    # modified primers" (ref. 23, Buenrostro 2015 Nature MOESM36 v2_Ad1/v2_Ad2).
    ad1 = graph.input("ad1_primer", "indexed v2_Ad1 amplification primer family", MaterialKind.OLIGO,state=primer_states()[0])
    ad2 = graph.input("ad2_primer", "indexed v2_Ad2 amplification primer family", MaterialKind.OLIGO,state=primer_states()[1])
    polymerase = graph.input("polymerase", "NEBNext High-Fidelity 2x PCR Master Mix (cited 2015 workflow)",
                            MaterialKind.CATALYST)
    extended=graph.product("gap_filled","gap-filled accessible-DNA template",MaterialKind.MOLECULAR,state=gap_filled_template())
    graph.react("gap_fill","72°C5min extension of transposition ends",StepKind.EXTENSION,
        ((purified,InputRole.TEMPLATE),(polymerase,InputRole.CATALYST)),(extended,))
    forward,reverse=("".join(x.top for x in con) for con in primer_states())
    preamp = graph.product("preamp", "five-cycle pre-amplified ATAC library",
                           MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    graph.react("pre_amplify", "Five pre-amplification cycles after gap filling", StepKind.PCR,
                ((extended, InputRole.TEMPLATE), (ad1, InputRole.STRUCTURAL),
                 (ad2, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (preamp,), geometry=PCRGeometry.exact(extended.state,PCRPrimer("Ad1",forward),PCRPrimer("Ad2",reverse),
                    statement="Source-qualified primer-family footprints are checked on the gap-filled precursor; this relation does not validate transposition or gap filling."))
    qPCR = graph.input("qpcr_rule", "5 µl SYBR Green I qPCR side reaction; N = cycle at ¼ maximum fluorescence",
                       MaterialKind.CONDITION)
    remaining = graph.product("remaining_reaction", "remaining 45 µl pre-amplified reaction",
                              MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    qc = graph.product("qpcr_aliquot", "consumed qPCR aliquot", MaterialKind.QC)
    graph.react("choose_cycles", "Use aliquot qPCR to choose the additional cycles",
                StepKind.SAMPLING,
                ((preamp, InputRole.TEMPLATE), (qPCR, InputRole.CONDITION)),
                (remaining, qc), changes_molecule=False)
    graph.mark_terminal(Terminal(qc.id, Disposition.QC))
    amplified = graph.product("amplified", "amplified, complexity-preserving ATAC library",
                              MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    graph.react("amplify_n", "Run the remaining reaction for the additional N cycles", StepKind.PCR,
                ((remaining, InputRole.TEMPLATE), (ad1, InputRole.STRUCTURAL),
                 (ad2, InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)), (amplified,), changes_molecule=False,
                geometry=PCRGeometry.exact(remaining.state,PCRPrimer("Ad1",forward),PCRPrimer("Ad2",reverse),
                    statement="Continued PCR uses the actual preamplified remainder with the same primer-family architectures; physical index members remain unassigned."))
    final_column = graph.input("final_column", "QIAGEN MinElute PCR purification column (cited 2015 workflow)",
                               MaterialKind.SELECTION)
    final = graph.product("library", "purified ATAC sequencing library (20 µl elution)",
                          MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    pcr_waste = graph.product("pcr_waste", "removed primers, polymerase and PCR reagents",
                              MaterialKind.WASTE)
    graph.react("final_cleanup", "Purify the amplified library", StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (final_column, InputRole.SELECTION)),
                (final, pcr_waste))
    graph.discard(pcr_waste, "is removed during the post-PCR column cleanup")
    graph.sequence(final, SequencingPlan(library=FINAL_LIBRARY,
                                         primers=tuple(SEQ_PRIMERS),acquired_roles=tuple(p.role for p in SEQ_PRIMERS),read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope([
        "cells", "pellet_condition", "pellet_cells", "tde1_adapter", "tn5", "digitonin_td",
        "tagment", "cleanup_column", "cleanup", "ad1_primer", "ad2_primer",
        "polymerase", "gap_fill", "pre_amplify", "qpcr_rule", "choose_cycles", "amplify_n",
        "final_column", "final_cleanup", "library"], source)
    return graph.require_publishable()
