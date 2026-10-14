"""NEBNext Single Cell/Low Input RNA Library Prep, E6420."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Row, Segment
from adapter_ligation import reported_adapter_end
from endprep import dA_tailed_scene, repair_and_dA_tail
from nebnext import STANDARD_PRIMERS, indexed_library
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
TITLE="NEBNext Single Cell / Low Input RNA"
NOTES="01_nebnext-single-cell-low-input-rna.html"
SOURCE='Commercial protocol: <a href="https://www.neb.com/products/e6420-nebnext-single-cell-low-input-rna-library-prep-kit-for-illumina">NEB E6420 product manual</a>.'
SUMMARY="A published oligo-dT RT primer and template-switching oligo generate full-length amplifiable cDNA; Ultra II FS then fragments and converts that cDNA into a dual-index Illumina library."
RT_PRIMER="AAGCAGTGGTATCAACGCAGAGTAC" + "T"*30 + "V"
TSO="GCTAATCATTGCAAGCAGTGGTATCAACGCAGAGTACATrGrGrG"
CDNA_PCR_PRIMER="AAGCAGTGGTATCAACGCAGAGT"
RT=[Row(chunks=[("mRNA  5′-[transcript]-poly(A)-3′",None,False)]), Row(chunks=[("RT primer  5′-"+RT_PRIMER+"-3′", "r1",True)])]
TS=[Row(chunks=[("TSO  5′-"+TSO+"-3′ → template switch at the mRNA 5′ end", "r2",True)]), Row(chunks=[("TSO-derived handle — full-length first-strand cDNA — RT-primer handle",None,False)])]
WTA=Construct([Segment("TSO-derived handle", "GCTAATCATTGC"+CDNA_PCR_PRIMER, "r2"), Segment("full-length cDNA", "X"*48, placeholder=True), Segment("RT-primer-derived handle", RT_PRIMER[:25], "r1")], name="amplified cDNA")
FS_END_PREP=repair_and_dA_tail(Construct([Segment("FS-fragmented cDNA", "X"*42, placeholder=True)]))
FS_ADAPTOR_END=reported_adapter_end(
    "NEBNext adaptor", kind="cohesive", polarity="3-prime", terminal_sequence="T",
    substrate_scenes=(dA_tailed_scene(FS_END_PREP),),
    statement="Ultra II FS produces 3′-dA cDNA ends for the complementary NEBNext adaptor 3′-dT; the kit-specific adaptor body is not disclosed here.")
STEPS=(("Prime reverse transcription", RT, "The E6420 manual prints the RT primer's 30-T tract and terminal V, as well as its PCR handle."),
       ("Template-switch at the transcript 5′ end", TS, "Template switching installs a second PCR handle, enabling full-length cDNA amplification."),
       ("Amplify full-length cDNA", WTA, "The manual prints the cDNA PCR primer; the transcript-specific insert remains variable."),
       ("Ultra II FS fragmentation and library conversion", Construct([Segment("FS-fragmented cDNA", "X"*42, placeholder=True)], name="FS cDNA fragment"), "Enzymatic fragmentation, end repair and dA-tailing precede NEBNext adaptor ligation and indexed PCR.", FS_ADAPTOR_END))
SEQ_PRIMERS=STANDARD_PRIMERS; REQUIRED_ROLES=tuple(p.role for p in SEQ_PRIMERS)
READOUT="Paired reads cover fragments distributed across the amplified full-length cDNA; the two index reads identify the library, not the cell."
FINAL_CAPTION="Dual-index library after Ultra II FS conversion. The published RT/TSO handles are generally internal to the pre-amplified cDNA and are not asserted at every final fragment end."
def final_library(): return indexed_library(Construct([Segment("cDNA insert", "X"*42, placeholder=True)]), "NEBNext low-input RNA library")

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    rna=g.input("rna","single-cell or low-input polyadenylated RNA",MaterialKind.BIOLOGICAL,state=RT,primary=True)
    rtprimer=g.input("rt_primer","published poly(A)-directed RT primer: "+RT_PRIMER,MaterialKind.OLIGO,source=SOURCE+" (Oligo Sequences table)")
    rt=g.input("rt","template-switching reverse transcriptase",MaterialKind.CATALYST)
    first=g.product("first_strand","nascent first-strand cDNA",MaterialKind.MOLECULAR,state=RT)
    g.react("reverse_transcribe","Prime and synthesize first-strand cDNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(rtprimer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,))
    tso=g.input("tso","published template-switching oligo: "+TSO,MaterialKind.OLIGO,source=SOURCE+" (Oligo Sequences table)")
    switched=g.product("switched","two-handle full-length cDNA",MaterialKind.MOLECULAR,state=WTA)
    g.react("template_switch","Template-switch at the transcript 5-prime end",StepKind.TEMPLATE_SWITCHING,((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,))
    wta_fwd=g.input("wta_forward","cDNA amplification primer at TSO end: "+CDNA_PCR_PRIMER,MaterialKind.OLIGO,source=SOURCE+" (one primer species, Oligo Sequences table)")
    wta_rev=g.input("wta_reverse","same cDNA amplification primer at RT end: "+CDNA_PCR_PRIMER,MaterialKind.OLIGO,source=SOURCE+" (one primer species, Oligo Sequences table)")
    wta_pol=g.input("wta_polymerase","cDNA amplification polymerase",MaterialKind.CATALYST)
    amplified=g.product("amplified","amplified full-length cDNA",MaterialKind.MOLECULAR,state=WTA)
    g.react("wta","Amplify full-length cDNA",StepKind.PCR,((switched,InputRole.TEMPLATE),(wta_fwd,InputRole.STRUCTURAL),(wta_rev,InputRole.STRUCTURAL),(wta_pol,InputRole.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(WTA,reason="The transcript-dependent full-length cDNA insert has no single exact primer-binding coordinate."))
    fs=g.input("fs_mix","Ultra II FS fragmentation/end-prep mix",MaterialKind.CATALYST)
    fragmented=g.product("fragmented","FS-fragmented amplified cDNA",MaterialKind.MOLECULAR,state=FS_END_PREP.insert)
    g.react("fragment","Fragment amplified cDNA",StepKind.FRAGMENTATION,((amplified,InputRole.TEMPLATE),(fs,InputRole.CATALYST)),(fragmented,))
    repaired=g.product("repaired","end-repaired cDNA fragments",MaterialKind.MOLECULAR,state=FS_END_PREP.insert)
    g.react("end_repair","Repair fragmented cDNA ends",StepKind.END_REPAIR,((fragmented,InputRole.TEMPLATE),(fs,InputRole.CATALYST)),(repaired,))
    tailed=g.product("tailed","FS-fragmented, repaired and dA-tailed cDNA",MaterialKind.MOLECULAR,state=dA_tailed_scene(FS_END_PREP))
    g.react("da_tail","Add one 3-prime dA to each repaired cDNA end",StepKind.TAILING,((repaired,InputRole.TEMPLATE),(fs,InputRole.CATALYST)),(tailed,))
    adapter=g.input("adapter","NEBNext library adapter",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (kit-specific body not disclosed)")
    ligase=g.input("ligase","NEBNext adapter ligase",MaterialKind.CATALYST)
    adapted=g.product("adapted","adapter-ligated cDNA fragments",MaterialKind.MOLECULAR,state=final_library(),inferred=True,source=SOURCE+" (kit-specific adapter body not disclosed)")
    g.react("ligate","Ligate the NEBNext library adapter",StepKind.LIGATION,((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(adapted,),geometry=FS_ADAPTOR_END)
    fwd=g.input("index_forward","dual-index PCR forward primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (sequence not disclosed)")
    rev=g.input("index_reverse","dual-index PCR reverse primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (sequence not disclosed)")
    pol=g.input("index_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=g.product("library","dual-index low-input RNA library",MaterialKind.LIBRARY,state=final_library())
    g.react("index_pcr","Complete the dual-index library",StepKind.PCR,((adapted,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),(rev,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(final_library(),reason="The kit-specific index-primer sequences are not disclosed."))
    g.sequence(library,SequencingPlan(library=final_library(),primers=tuple(SEQ_PRIMERS),acquired_roles=REQUIRED_ROLES))
    g.declare_reviewed_scope(["rna","rt_primer","rt","reverse_transcribe","tso","template_switch","wta_forward","wta_reverse","wta_polymerase","wta","fs_mix","fragment","end_repair","da_tail","adapter","ligase","ligate","index_forward","index_reverse","index_polymerase","index_pcr","library"],SOURCE)
    return g.require_publishable()
