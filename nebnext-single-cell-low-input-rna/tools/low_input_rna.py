"""NEBNext Single Cell/Low Input RNA Library Prep, E6420."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Row, Segment, revcomp
from adapter_ligation import reported_adapter_end
from endprep import dA_tailed_scene, repair_and_dA_tail
from nebnext import STANDARD_PRIMERS, indexed_library
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
TITLE="NEBNext Single Cell / Low Input RNA"
NOTES="01_nebnext-single-cell-low-input-rna.html"
SOURCE='Commercial protocol: <a href="https://www.neb.com/products/e6420-nebnext-single-cell-low-input-rna-library-prep-kit-for-illumina">NEB E6420 manual, version 6.0_01/24</a>.'
SUMMARY="A published oligo-dT RT primer and template-switching oligo generate full-length amplifiable cDNA; Ultra II FS then fragments and converts that cDNA into a dual-index Illumina library."
RT_PRIMER="AAGCAGTGGTATCAACGCAGAGTAC" + "T"*30 + "V"
TSO="GCTAATCATTGCAAGCAGTGGTATCAACGCAGAGTACATrGrGrG"
CDNA_PCR_PRIMER="AAGCAGTGGTATCAACGCAGAGT"
RT=[Row(chunks=[("mRNA  5′-[transcript]-poly(A)-3′",None,False)]), Row(chunks=[("RT primer  5′-"+RT_PRIMER+"-3′", "r1",True)])]
TS=[Row(chunks=[("TSO  5′-"+TSO+"-3′ → template switch at the mRNA 5′ end", "r2",True)]), Row(chunks=[("TSO-derived handle — full-length first-strand cDNA — RT-primer handle",None,False)])]
TSO_DNA=TSO[:-len("rGrGrG")]
RT_HANDLE=RT_PRIMER[:25]
# Amplified cDNA in transcript (sense) orientation.  The RT-primer handle appears as its
# reverse complement on this strand, so the one printed cDNA PCR primer lands inward at
# both ends (checked by PCRGeometry.symmetric below).  GGG is copied from the TSO rGrGrG
# and the B base pairs with the primer's terminal V: both are derived, not printed.
def _cdna(name, *, tso_5prime_extension):
    left=[Segment("cDNA PCR primer site (TSO)", CDNA_PCR_PRIMER, "r2"),
          Segment("TSO ACAT", TSO_DNA[len(TSO_DNA)-4:], "r2"),
          Segment("GGG copied from TSO rGrGrG", "GGG", "r2", inferred=True)]
    if tso_5prime_extension:
        left.insert(0, Segment("TSO 5′ extension", TSO_DNA[:12], "r2"))
    return Construct([*left,
                      Segment("full-length cDNA", "X"*48, placeholder=True),
                      Segment("anchor base pairing V", "B", placeholder=True, inferred=True),
                      Segment("poly(A) copied by the T30 tract", "A"*30),
                      Segment("RT-primer handle (reverse complement)", revcomp(RT_HANDLE), "r1")], name=name)
SWITCHED=_cdna("template-switched cDNA (drawn as duplex PCR template)", tso_5prime_extension=True)
WTA=_cdna("amplified cDNA", tso_5prime_extension=False)
WTA_GEOMETRY=PCRGeometry.symmetric(SWITCHED, PCRPrimer("NEBNext Single Cell cDNA PCR Primer", CDNA_PCR_PRIMER),
    statement="The one printed cDNA PCR primer is the 5′ 23 nt of both the TSO DNA body and the RT primer handle, so it primes inward from both cDNA ends; products lack the 12-nt TSO 5′ extension.")
FS_END_PREP=repair_and_dA_tail(Construct([Segment("FS-fragmented cDNA", "X"*42, placeholder=True)]))
FS_ADAPTOR_END=reported_adapter_end(
    "NEBNext adaptor", kind="cohesive", polarity="3-prime", terminal_sequence="T",
    substrate_scenes=(dA_tailed_scene(FS_END_PREP),),
    statement="Ultra II FS produces 3′-dA cDNA ends for the complementary NEBNext adaptor 3′-dT; the kit-specific adaptor body is not disclosed here.")
SEQ_PRIMERS=STANDARD_PRIMERS; REQUIRED_ROLES=tuple(p.role for p in SEQ_PRIMERS)
READOUT="Paired reads (typically 2×75, manual Read Length and Depth) cover fragments distributed across the amplified full-length cDNA. Index reads identify the library, not a cell; the index layout depends on the separately purchased NEBNext oligo kit (the manual also allows single-index kits with a universal i5-side primer), so the dual 8-base index drawn here is one option."
FINAL_CAPTION="INFERRED — dual-index option after Ultra II FS conversion; the adaptor and index primers come from a separately chosen oligo kit not printed in E6420. The published RT/TSO handles are generally internal to the pre-amplified cDNA and are not asserted at every final fragment end."
READ_LENGTHS={"Read 1":75,"Read 2":75}
def final_library(): return indexed_library(Construct([Segment("cDNA insert", "X"*42, placeholder=True)]), "NEBNext low-input RNA library", inferred_adaptors=True)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    rna=g.input("rna","single-cell or low-input polyadenylated RNA",MaterialKind.BIOLOGICAL,state=RT,primary=True)
    rtprimer=g.input("rt_primer","published poly(A)-directed RT primer: "+RT_PRIMER,MaterialKind.OLIGO,source=SOURCE+" (Oligo Sequences table)")
    rt=g.input("rt","template-switching reverse transcriptase",MaterialKind.CATALYST)
    first=g.product("first_strand","nascent first-strand cDNA",MaterialKind.MOLECULAR,state=RT)
    g.react("reverse_transcribe","Prime and synthesize first-strand cDNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(rtprimer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,))
    tso=g.input("tso","published template-switching oligo: "+TSO,MaterialKind.OLIGO,state=Construct([Segment("TSO DNA body", TSO_DNA, "r2"), Segment("3′ riboguanosines rGrGrG", "GGG", "r2")],name="TSO"),source=SOURCE+" (Oligo Sequences table)")
    switched=g.product("switched","template-switched full-length first-strand cDNA",MaterialKind.MOLECULAR,state=SWITCHED)
    g.react("template_switch","Template-switch at the transcript 5-prime end",StepKind.TEMPLATE_SWITCHING,((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,))
    wta_primer=g.input("wta_primer","NEBNext Single Cell cDNA PCR Primer used at both cDNA ends: "+CDNA_PCR_PRIMER,MaterialKind.OLIGO,state=Construct([Segment("cDNA PCR primer", CDNA_PCR_PRIMER, "r2")],name="cDNA PCR primer"),source=SOURCE+" (Oligo Sequences table)")
    wta_pol=g.input("wta_polymerase","NEBNext Single Cell cDNA PCR Master Mix (11-21 cycles, 62C anneal, 3 min extension)",MaterialKind.CATALYST)
    amplified=g.product("amplified","amplified full-length cDNA",MaterialKind.MOLECULAR,state=WTA)
    g.react("wta","Amplify full-length cDNA",StepKind.PCR,((switched,InputRole.TEMPLATE),(wta_primer,InputRole.STRUCTURAL),(wta_pol,InputRole.CATALYST)),(amplified,),geometry=WTA_GEOMETRY)
    cdna_beads=g.input("cdna_beads","SPRIselect/AMPure beads 0.6X, then re-bound with Bead Reconstitution Buffer",MaterialKind.SELECTION)
    cdna_clean=g.product("cdna_clean","purified amplified cDNA",MaterialKind.MOLECULAR,state=WTA)
    cdna_waste=g.product("cdna_cleanup_waste","supernatants of the two cDNA bead cleanups (primers, short products)",MaterialKind.WASTE)
    g.react("cdna_cleanup","Double bead cleanup of amplified cDNA (manual 1.6)",StepKind.CLEANUP,((amplified,InputRole.TEMPLATE),(cdna_beads,InputRole.SELECTION)),(cdna_clean,cdna_waste),changes_molecule=False)
    g.discard(cdna_waste,"supernatants are discarded at steps 1.6.5 and 1.6.14")
    fs=g.input("fs_mix","Ultra II FS fragmentation/end-prep mix",MaterialKind.CATALYST)
    fragmented=g.product("fragmented","FS-fragmented amplified cDNA",MaterialKind.MOLECULAR,state=FS_END_PREP.insert)
    g.react("fragment","Fragment amplified cDNA (Ultra II FS, 25 min 37C then 30 min 65C)",StepKind.FRAGMENTATION,((cdna_clean,InputRole.TEMPLATE),(fs,InputRole.CATALYST)),(fragmented,))
    repaired=g.product("repaired","end-repaired cDNA fragments",MaterialKind.MOLECULAR,state=FS_END_PREP.insert)
    g.react("end_repair","Repair fragmented cDNA ends",StepKind.END_REPAIR,((fragmented,InputRole.TEMPLATE),(fs,InputRole.CATALYST)),(repaired,))
    tailed=g.product("tailed","FS-fragmented, repaired and dA-tailed cDNA",MaterialKind.MOLECULAR,state=dA_tailed_scene(FS_END_PREP))
    g.react("da_tail","Add one 3-prime dA to each repaired cDNA end",StepKind.TAILING,((repaired,InputRole.TEMPLATE),(fs,InputRole.CATALYST)),(tailed,))
    adapter=g.input("adapter","NEBNext library adapter",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (kit-specific body not disclosed)")
    ligase=g.input("ligase","NEBNext Ultra II Ligation Master Mix plus Ligation Enhancer (20C 15 min)",MaterialKind.CATALYST)
    adapted=g.product("adapted","adapter-ligated cDNA fragments",MaterialKind.MOLECULAR,state=final_library(),inferred=True,source=SOURCE+" (kit-specific adapter body not disclosed)")
    g.react("ligate","Ligate the NEBNext library adapter",StepKind.LIGATION,((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(adapted,),geometry=FS_ADAPTOR_END)
    user=g.input("user","USER Enzyme from the oligo kit (3 ul into the ligation, 37C 15 min)",MaterialKind.CATALYST)
    user_treated=g.product("user_treated","USER-treated adapter-ligated cDNA fragments",MaterialKind.MOLECULAR,state=final_library(),inferred=True,source=SOURCE+" (step 1.9.6 requires USER for NEBNext Adaptors but does not state the adaptor structure or what USER cleaves)")
    g.react("user_treat","Add USER Enzyme to the ligation (manual 1.9.6-1.9.7)",StepKind.DIGESTION,((adapted,InputRole.TEMPLATE),(user,InputRole.CATALYST)),(user_treated,),changes_molecule=False)
    lig_beads=g.input("ligation_beads","SPRIselect/AMPure beads 0.8X",MaterialKind.SELECTION)
    lig_clean=g.product("ligation_clean","purified adapter-ligated cDNA fragments",MaterialKind.MOLECULAR,state=final_library(),inferred=True,source=SOURCE+" (adapter body not disclosed)")
    lig_waste=g.product("ligation_cleanup_waste","ligation supernatant (free adaptor, short products)",MaterialKind.WASTE)
    g.react("ligation_cleanup","Bead cleanup of the ligation (manual 1.10)",StepKind.CLEANUP,((user_treated,InputRole.TEMPLATE),(lig_beads,InputRole.SELECTION)),(lig_clean,lig_waste),changes_molecule=False)
    g.discard(lig_waste,"supernatant is discarded at step 1.10.5")
    fwd=g.input("index_forward","dual-index PCR forward primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (sequence not disclosed)")
    rev=g.input("index_reverse","dual-index PCR reverse primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (sequence not disclosed)")
    pol=g.input("index_polymerase","NEBNext Ultra II Q5 Master Mix (typically 8 cycles, 65C anneal/extension)",MaterialKind.CATALYST)
    pcr_product=g.product("pcr_product","indexed low-input RNA library before cleanup",MaterialKind.MOLECULAR,state=final_library(),inferred=True,source=SOURCE+" (oligo-kit sequences not printed)")
    g.react("index_pcr","Complete the indexed library",StepKind.PCR,((lig_clean,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),(rev,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(pcr_product,),geometry=PCRGeometry.unavailable(final_library(),reason="The kit-specific index-primer sequences are not disclosed."))
    pcr_beads=g.input("pcr_beads","SPRIselect/AMPure beads 0.9X",MaterialKind.SELECTION)
    library=g.product("library","indexed low-input RNA library",MaterialKind.LIBRARY,state=final_library(),inferred=True,source=SOURCE+" (oligo-kit sequences not printed)")
    pcr_waste=g.product("pcr_cleanup_waste","PCR supernatant (primers; ~80 bp primer and 128-140 bp adaptor-dimer peaks prompt a repeat cleanup)",MaterialKind.WASTE)
    g.react("pcr_cleanup","Bead cleanup of the PCR (manual 1.12)",StepKind.CLEANUP,((pcr_product,InputRole.TEMPLATE),(pcr_beads,InputRole.SELECTION)),(library,pcr_waste),changes_molecule=False)
    g.discard(pcr_waste,"supernatant is discarded at step 1.12.5")
    g.sequence(library,SequencingPlan(library=final_library(),primers=tuple(SEQ_PRIMERS),acquired_roles=REQUIRED_ROLES,read_lengths=READ_LENGTHS))
    g.declare_reviewed_scope(["rna","rt_primer","rt","reverse_transcribe","tso","template_switch","wta_primer","wta_polymerase","wta","cdna_beads","cdna_cleanup","fs_mix","fragment","end_repair","da_tail","adapter","ligase","ligate","user","user_treat","ligation_beads","ligation_cleanup","index_forward","index_reverse","index_polymerase","index_pcr","pcr_beads","pcr_cleanup","library"],SOURCE)
    return g.require_publishable()
