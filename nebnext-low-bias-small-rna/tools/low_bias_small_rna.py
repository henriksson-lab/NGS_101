"""NEBNext Low-bias Small RNA library preparation, E3420 manual v2.0 9/25."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Row, Segment
from adapter_ligation import AdapterLigationGeometry
from nebnext import STANDARD_PRIMERS, indexed_library
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
TITLE="NEBNext Low-bias Small RNA"
NOTES="01_nebnext-low-bias-small-rna.html"
SOURCE='Commercial protocol: <a href="https://www.neb.com/en-us/products/e3420-nebnext-low-bias-small-rna-library-prep-kit">NEB E3420 manual v2.0 (September 2025)</a>.'
SUMMARY="Small RNAs with 5′ phosphate and 3′ hydroxyl are captured through randomized splint ligation, followed by selective adaptor removal, reverse transcription and UDI PCR."
CAVEAT="The current manual does not disclose adaptor bases. The final panel uses the documented TruSeq-compatible primer-site boundary while rendering every undisclosed base as a placeholder."
INPUT=Construct([Segment("5′-phosphorylated small RNA", "R"*28, placeholder=True)], name="small RNA")
THREE=Construct([Segment("small RNA", "R"*28, placeholder=True), Segment("3′ adaptor", "X"*14, "r2", placeholder=True, inferred=True)], name="3′-ligated RNA")
SPLINT=[Row(chunks=[("5′ adaptor ** small RNA ** 3′ adaptor",None,False)]), Row(chunks=[("       randomized splint pairs across the ligation substrate", "me",True)])]
CDNA=Construct([Segment("5′ adaptor remnant", "X"*12, "r1", placeholder=True, inferred=True), Segment("small-RNA cDNA", "X"*28, placeholder=True), Segment("3′ adaptor remnant", "X"*12, "r2", placeholder=True, inferred=True)], name="small-RNA cDNA")
THREE_END=AdapterLigationGeometry.unavailable(
    "The E3420 manual identifies the RNA 3′-OH acceptor but does not disclose the 3′-adaptor terminal bases or activation chemistry.")
FIVE_END=AdapterLigationGeometry.unavailable(
    "The E3420 manual identifies a randomized-splint 5′-adaptor ligation but does not disclose the adaptor or splint strands.")
STEPS=(("Ligate the 3′ adaptor", THREE, "The acceptor is an RNA 3′ hydroxyl. Excess free adaptor is enzymatically removed before the next ligation.", THREE_END),
       ("Randomized-splint 5′ adaptor ligation", SPLINT, "The splint reduces sequence-dependent ligation bias. Exact adaptor and splint bases are not disclosed in the manual.", FIVE_END),
       ("Modify the 3′ adaptor and reverse-transcribe", CDNA, "During 5′ ligation the existing 3′ adaptor is enzymatically converted into an RT-priming substrate; RT copies the captured RNA."),
       ("Size-select and UDI-amplify", Construct([Segment("selected small-RNA cDNA", "X"*40, placeholder=True)], name="selected cDNA"), "Bead selection can emphasize miRNAs or retain a broader small-RNA range before dual-index PCR."))
SEQ_PRIMERS=STANDARD_PRIMERS; REQUIRED_ROLES=tuple(p.role for p in SEQ_PRIMERS)
READOUT="A 56-base Read 1 is normally sufficient for the small-RNA insert; i5 and i7 are read separately from the NEBNext LV unique-dual-index library."
FINAL_CAPTION="INFERRED — TruSeq-compatible UDI boundary. Proprietary ligated-adaptor remnants remain dotted; the standard sequencing-primer sites are computed from the completed PCR product."
def final_library():
    ins=Construct([Segment("5′ kit remnant", "X"*8, placeholder=True, inferred=True), Segment("small-RNA insert", "X"*28, placeholder=True), Segment("3′ kit remnant", "X"*8, placeholder=True, inferred=True), Segment("3′ adaptor boundary A", "A", inferred=True, note="required by the documented TruSeq-compatible Read 2 site; exact proprietary adaptor sequence is unavailable")])
    return indexed_library(ins, "NEBNext low-bias small-RNA library",
                           inferred_adaptors=True, dA_junction=False)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    rna=g.input("rna","5-prime-phosphorylated, 3-prime-hydroxyl small RNA",MaterialKind.BIOLOGICAL,state=INPUT,primary=True)
    adapter3=g.input("adapter3","proprietary 3-prime RNA adapter",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (bases and activation chemistry not disclosed)")
    ligase3=g.input("ligase3","3-prime RNA ligase",MaterialKind.CATALYST)
    three=g.product("three_ligated","3-prime-adapter-ligated small RNA",MaterialKind.MOLECULAR,state=THREE,inferred=True,source=SOURCE+" (adapter bases not disclosed)")
    g.react("ligate_three","Ligate the 3-prime adapter",StepKind.LIGATION,((rna,InputRole.TEMPLATE),(adapter3,InputRole.STRUCTURAL),(ligase3,InputRole.CATALYST)),(three,),geometry=THREE_END)
    cleanup=g.input("adapter_cleanup","free-adapter removal enzyme",MaterialKind.CATALYST)
    cleaned=g.product("cleaned","3-prime-ligated RNA after free-adapter removal",MaterialKind.MOLECULAR,state=three.state)
    free=g.product("free_adapter_waste","destroyed excess 3-prime adapter",MaterialKind.WASTE)
    g.react("remove_free_adapter","Remove excess free adapter",StepKind.DIGESTION,((three,InputRole.TEMPLATE),(cleanup,InputRole.CATALYST)),(cleaned,free),changes_molecule=False)
    g.discard(free,"is enzymatically removed before 5-prime adapter ligation")
    adapter5=g.input("adapter5","proprietary 5-prime RNA adapter",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (bases not disclosed)")
    splint=g.input("splint","randomized ligation splint",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (bases not disclosed)")
    ligase5=g.input("ligase5","5-prime RNA ligase",MaterialKind.CATALYST)
    dual=g.product("dual_ligated","two-adapter small RNA",MaterialKind.MOLECULAR,state=CDNA,inferred=True,source=SOURCE+" (adapter boundaries not disclosed)")
    g.react("ligate_five","Randomized-splint 5-prime adapter ligation",StepKind.LIGATION,((cleaned,InputRole.TEMPLATE),(adapter5,InputRole.STRUCTURAL),(splint,InputRole.STRUCTURAL),(ligase5,InputRole.CATALYST)),(dual,),geometry=FIVE_END)
    rtprimer=g.input("rt_primer","adapter-complementary RT primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (sequence not disclosed)")
    rt=g.input("rt","reverse transcriptase and adapter-modification enzymes",MaterialKind.CATALYST)
    cdna=g.product("cdna","small-RNA cDNA",MaterialKind.MOLECULAR,state=CDNA,inferred=True,source=SOURCE+" (adapter remnants not disclosed)")
    g.react("reverse_transcribe","Modify the 3-prime adapter and reverse-transcribe",StepKind.REVERSE_TRANSCRIPTION,((dual,InputRole.TEMPLATE),(rtprimer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
    beads=g.input("beads","small-RNA size-selection beads",MaterialKind.SELECTION)
    selected=g.product("selected","size-selected small-RNA cDNA",MaterialKind.MOLECULAR,state=cdna.state)
    excluded=g.product("excluded_sizes","excluded cDNA sizes",MaterialKind.WASTE)
    g.react("size_select","Select the requested small-RNA size range",StepKind.SELECTION,((cdna,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(selected,excluded),changes_molecule=False)
    g.discard(excluded,"falls outside the selected small-RNA library size range")
    fwd=g.input("pcr_forward","NEBNext UDI PCR forward primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (sequence not disclosed)")
    rev=g.input("pcr_reverse","NEBNext UDI PCR reverse primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (sequence not disclosed)")
    pol=g.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=g.product("library","NEBNext low-bias small-RNA library",MaterialKind.LIBRARY,state=final_library(),inferred=True,source=SOURCE+" (proprietary adapter remnants)")
    g.react("pcr","UDI-amplify the selected cDNA",StepKind.PCR,((selected,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),(rev,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(final_library(),reason="The manual does not disclose complete UDI primer sequences."))
    g.sequence(library,SequencingPlan(library=final_library(),primers=tuple(SEQ_PRIMERS),acquired_roles=REQUIRED_ROLES,read_lengths={"Read 1":56}))
    g.declare_reviewed_scope(["rna","adapter3","ligase3","ligate_three","adapter_cleanup","remove_free_adapter","adapter5","splint","ligase5","ligate_five","rt_primer","rt","reverse_transcribe","beads","size_select","pcr_forward","pcr_reverse","pcr_polymerase","pcr","library"],SOURCE)
    return g.require_publishable()
