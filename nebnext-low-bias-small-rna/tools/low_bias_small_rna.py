"""NEBNext Low-bias Small RNA library preparation, E3420 manual v2.0 5/26."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Segment, revcomp
import illumina as il
from adapter_ligation import AdapterLigationGeometry
from nebnext import STANDARD_PRIMERS, indexed_library
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
TITLE="NEBNext Low-bias Small RNA"
NOTES="01_nebnext-low-bias-small-rna.html"
SOURCE='Commercial protocol: <a href="https://www.neb.com/en/-/media/nebus/files/manuals/manuale3420.pdf">NEB E3420 manual v2.0 (May 2026)</a>.'
RT_SOURCE='NEB <a href="https://www.neb.com/en-us/faqs/when-is-the-rt-primer-annealed-during-the-nebnext-low-bias-small-rna-library-prep-workflow-neb-e3420">E3420 RT-primer FAQ</a>.'
SUMMARY="Small RNAs (< 120 nt) with 5′ phosphate and 3′ hydroxyl are captured by 3′-adaptor ligation, enzymatic removal of excess 3′ adaptor, randomized-splint 5′-adaptor ligation, reverse transcription, a bead selection that removes large cDNA and unligated adaptor, LV UDI PCR and a choice of miRNA-focused size selection or broad cleanup."
CAVEAT="The checked E3420 manual does not print adaptor, splint or UDI primer strands. NEB cites an earlier research protocol with published oligos but says the commercial kit has since been enhanced; equivalence has not been established, and the predecessor 5′ adapter would not give the TruSeq Read 2 trim motif printed by the E3420 manual. The final library is derived from the manual's instruction to trim like standard TruSeq libraries (Read 1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA, Read 2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT): the insert is placed directly between the TruSeq Read 1 and Index 1 arms. Index wells and ligation intermediates remain placeholders."
INPUT=Construct([Segment("5′-phosphorylated small RNA", "R"*28, placeholder=True)], name="small RNA")
THREE=Construct([Segment("small RNA", "R"*28, placeholder=True), Segment("3′ adaptor", "X"*14, "r2", placeholder=True, inferred=True)], name="3′-ligated RNA")
DUAL=Construct([Segment("5′ adaptor", "X"*12, "r1", placeholder=True, inferred=True), Segment("small RNA", "R"*28, placeholder=True), Segment("3′ adaptor", "X"*14, "r2", placeholder=True, inferred=True)], name="two-adaptor small RNA")
CDNA=Construct([Segment("5′ adaptor remnant", "X"*12, "r1", placeholder=True, inferred=True), Segment("small-RNA cDNA", "X"*28, placeholder=True), Segment("3′ adaptor remnant", "X"*12, "r2", placeholder=True, inferred=True)], name="small-RNA cDNA")
THREE_END=AdapterLigationGeometry.unavailable(
    "The E3420 manual identifies the RNA 3′-OH acceptor but does not disclose the 3′-adaptor terminal bases or activation chemistry.")
FIVE_END=AdapterLigationGeometry.unavailable(
    "The E3420 manual identifies a randomized-splint 5′-adaptor ligation but does not disclose the adaptor or splint strands.")
SEQ_PRIMERS=STANDARD_PRIMERS; REQUIRED_ROLES=tuple(p.role for p in SEQ_PRIMERS)
READOUT="The manual states that a single 56-base read is typically sufficient and that longer or paired-end reads can be used; the LV UDI primers carry 8-base i5 and i7 indices (well chosen per sample)."
FINAL_CAPTION="INFERRED — TruSeq-shaped UDI library derived from the E3420 trim motifs (manual, Sequencing and Data Processing); adaptor bases are not printed as oligos, so all arms remain dotted. Sequencing-primer sites are computed from the completed PCR product."
# Verbatim from the E3420 manual v2.0_5/26, "Sequencing and Data Processing".
MANUAL_TRIM={"Read 1":"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA","Read 2":"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"}

def final_library():
    ins=Construct([Segment("small-RNA insert", "X"*22, placeholder=True),
                   Segment("3′ adaptor boundary A", "A", inferred=True,
                           note="first base of the manual's Read 1 trim motif; E3420 adaptor oligo not printed")])
    lib=indexed_library(ins, "NEBNext low-bias small-RNA library",
                        inferred_adaptors=True, dA_junction=False)
    seq=lib.top(); i=seq.index("X"*22); j=i+22
    readthrough={"Read 1":seq[j:j+33],"Read 2":revcomp(seq[:i])[:33]}
    if readthrough!=MANUAL_TRIM:
        raise ValueError(f"final library read-through {readthrough} disagrees with E3420 trim motifs")
    return lib

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    rna=g.input("rna","5-prime-phosphorylated, 3-prime-hydroxyl small RNA",MaterialKind.BIOLOGICAL,state=INPUT,primary=True)
    adapter3=g.input("adapter3","E3420 3-prime RNA adapter",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (bases and activation chemistry not disclosed in checked kit sources)")
    ligase3=g.input("ligase3","3-prime RNA ligase",MaterialKind.CATALYST)
    three=g.product("three_ligated","3-prime-adapter-ligated small RNA",MaterialKind.MOLECULAR,state=THREE,inferred=True,source=SOURCE+" (adapter bases not disclosed)")
    g.react("ligate_three","Ligate the 3-prime adapter",StepKind.LIGATION,((rna,InputRole.TEMPLATE),(adapter3,InputRole.STRUCTURAL),(ligase3,InputRole.CATALYST)),(three,),geometry=THREE_END)
    cleanup=g.input("adapter_cleanup","free-adapter removal enzyme",MaterialKind.CATALYST)
    cleaned=g.product("cleaned","3-prime-ligated RNA after free-adapter removal",MaterialKind.MOLECULAR,state=three.state)
    free=g.product("free_adapter_waste","destroyed excess 3-prime adapter",MaterialKind.WASTE)
    g.react("remove_free_adapter","Remove excess free adapter",StepKind.DIGESTION,((three,InputRole.TEMPLATE),(cleanup,InputRole.CATALYST)),(cleaned,free),changes_molecule=False)
    g.discard(free,"is enzymatically removed before 5-prime adapter ligation")
    adapter5=g.input("adapter5","E3420 5-prime RNA adapter",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (bases not disclosed in checked kit sources)")
    splint=g.input("splint","randomized ligation splint",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (bases not disclosed)")
    ligase5=g.input("ligase5","NEBNext LB 5-prime Ligase",MaterialKind.CATALYST)
    dual=g.product("dual_ligated","two-adapter small RNA with 3-prime-adaptor-derived RT primer",MaterialKind.MOLECULAR,state=DUAL,inferred=True,source=SOURCE+" "+RT_SOURCE+" (adapter boundaries not disclosed)")
    g.react("ligate_five","Ligate the 5-prime adapter and generate the 3-prime-adapter RT primer",StepKind.LIGATION,((cleaned,InputRole.TEMPLATE),(adapter5,InputRole.STRUCTURAL),(splint,InputRole.STRUCTURAL),(ligase5,InputRole.CATALYST)),(dual,),geometry=FIVE_END)
    rt=g.input("rt","NEBNext LB RT Enzyme Mix",MaterialKind.CATALYST)
    cdna=g.product("cdna","small-RNA cDNA",MaterialKind.MOLECULAR,state=CDNA,inferred=True,source=SOURCE+" (adapter remnants not disclosed)")
    g.react("reverse_transcribe","Reverse-transcribe from the 3-prime-adaptor-derived primer",StepKind.REVERSE_TRANSCRIPTION,((dual,InputRole.TEMPLATE),(dual,InputRole.SELF_PRIMER),(rt,InputRole.CATALYST)),(cdna,))
    beads=g.input("beads","NEBNext Sample Purification Beads, 1X then 1X plus isopropanol",MaterialKind.SELECTION)
    selected=g.product("selected","size-selected small-RNA cDNA",MaterialKind.MOLECULAR,state=cdna.state)
    excluded=g.product("excluded_sizes","large cDNA on first-bind beads; unligated adaptor in second-bind supernatant",MaterialKind.WASTE)
    g.react("size_select","Bead size-select the RT reaction (manual step 6)",StepKind.SELECTION,((cdna,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(selected,excluded),changes_molecule=False)
    g.discard(excluded,"large fragments are discarded on the 1X beads and non-ligated adaptor is removed with the second-bind supernatant (manual Figure 1 legend, steps 6.6 and 6.13)")
    fwd=g.input("pcr_forward","NEBNext UDI PCR forward primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (sequence not disclosed)")
    rev=g.input("pcr_reverse","NEBNext UDI PCR reverse primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (sequence not disclosed)")
    pol=g.input("pcr_polymerase","NEBNext Ultra II Q5 Master Mix (9-22 cycles, 62C annealing)",MaterialKind.CATALYST)
    amplified=g.product("amplified","UDI-amplified small-RNA library before bead purification",MaterialKind.MOLECULAR,state=final_library(),inferred=True,source=SOURCE+" (structure derived from the printed TruSeq trim motifs)")
    g.react("pcr","UDI-amplify the selected cDNA",StepKind.PCR,((selected,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),(rev,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(final_library(),reason="The manual does not disclose complete UDI primer sequences."))
    post=g.input("post_pcr_beads","NEBNext Sample Purification Beads: 8A 1.2X/0.8X double-sided size selection or 8B 2X cleanup",MaterialKind.SELECTION)
    library=g.product("library","NEBNext low-bias small-RNA library",MaterialKind.LIBRARY,state=final_library(),inferred=True,source=SOURCE+" (structure derived from the printed TruSeq trim motifs)")
    post_waste=g.product("post_pcr_waste","8A: large fragments on 1.2X beads and short products in 0.8X supernatant; 8B: unbound supernatant",MaterialKind.WASTE)
    g.react("post_pcr_purify","Purify the PCR product: 8A miRNA-focused size selection (< 50 nt inserts) or 8B cleanup (< 120 nt)",StepKind.SELECTION,((amplified,InputRole.TEMPLATE),(post,InputRole.SELECTION)),(library,post_waste),changes_molecule=False)
    g.discard(post_waste,"is removed by the chosen post-PCR bead purification (manual steps 8A.6/8A.11 or 8B.6)")
    g.sequence(library,SequencingPlan(library=final_library(),primers=tuple(SEQ_PRIMERS),acquired_roles=REQUIRED_ROLES,read_lengths={"Read 1":56}))
    g.declare_reviewed_scope(["rna","adapter3","ligase3","ligate_three","adapter_cleanup","remove_free_adapter","adapter5","splint","ligase5","ligate_five","rt","reverse_transcribe","beads","size_select","pcr_forward","pcr_reverse","pcr_polymerase","pcr","amplified","post_pcr_beads","post_pcr_purify","library"],SOURCE)
    return g.require_publishable()
