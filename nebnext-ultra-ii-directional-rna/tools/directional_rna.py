"""NEBNext Ultra II Directional RNA library preparation, E7760."""
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
TITLE="NEBNext Ultra II Directional RNA"
NOTES="01_nebnext-ultra-ii-directional-rna.html"
SOURCE='Commercial protocol: <a href="https://www.neb.com/en-us/products/next-generation-sequencing-library-preparation/library-preparation-for-illumina/rna-library-prep-for-illumina">NEB E7760 product manual and official RNA workflow</a>.'
SUMMARY="Fragmented RNA is copied into double-stranded cDNA with dUTP in strand two. USER selectively removes that marked strand before PCR, preserving transcript orientation."
RNA=Construct([Segment("fragmented RNA", "R"*42, placeholder=True)], name="RNA fragment")
FIRST=Construct([Segment("first-strand cDNA", "X"*42, placeholder=True)], name="RNA:cDNA product")
SECOND=[Row(chunks=[("first-strand cDNA  3′-xxxxxxxxxxxxxxxx-5′",None,False)]), Row(chunks=[("second strand      5′-xxxUxxxxUxxxxxxx-3′", "w1",False)])]
USER=[Row(chunks=[("USER cuts the dUTP-labelled second strand → it cannot seed library PCR",None,False)]), Row(chunks=[("the unmarked first strand defines the sequenced orientation",None,False)])]
END_PREP=repair_and_dA_tail(Construct([Segment("dUTP-marked cDNA duplex", "X"*42, placeholder=True)]))
ADAPTOR_END=reported_adapter_end(
    "NEBNext adaptor", kind="cohesive", polarity="3-prime", terminal_sequence="T",
    substrate_scenes=(dA_tailed_scene(END_PREP),),
    statement="The Ultra II directional workflow dA-tails the cDNA before ligating the complementary 3′-dT NEBNext adaptor; the selected oligo-set body varies by kit configuration.")
STEPS=(("Fragment RNA and synthesize the first cDNA strand", FIRST, "Random priming copies enriched or rRNA-depleted RNA fragments."),
       ("Make a dUTP-marked second strand", SECOND, "dUTP substitutes for dT during second-strand synthesis; the mark belongs only to strand two."),
       ("End-prep, ligate adaptor and remove strand two", USER, "USER excises the dU-containing strand. Strand-specificity reagent reinforces the same directional selection.", ADAPTOR_END),
       ("PCR-complete the surviving strand", Construct([Segment("orientation-preserving cDNA insert", "X"*42, placeholder=True)], name="directional insert"), "NEBNext dual-index primers amplify molecules derived from the surviving first strand."))
SEQ_PRIMERS=STANDARD_PRIMERS; REQUIRED_ROLES=tuple(p.role for p in SEQ_PRIMERS)
READOUT="Read 1 and Read 2 traverse opposite ends of the cDNA insert; their orientation is interpretable because the dUTP-labelled second strand was removed before amplification."
FINAL_CAPTION="Dual-index directional RNA library. The molecule has ordinary TruSeq-compatible sequencing sites; directionality was imposed earlier by selective second-strand destruction."
def final_library(): return indexed_library(Construct([Segment("directional cDNA insert", "X"*42, placeholder=True)]), "Ultra II Directional RNA library")

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    rna=g.input("rna","fragmented enriched or rRNA-depleted RNA",MaterialKind.BIOLOGICAL,state=RNA,primary=True)
    rtprimer=g.input("rt_primer","random RT primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (sequence not disclosed)")
    rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
    first=g.product("first_strand","RNA:first-strand cDNA",MaterialKind.MOLECULAR,state=FIRST)
    g.react("reverse_transcribe","Synthesize first-strand cDNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(rtprimer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,))
    second_primer=g.input("second_primer","second-strand synthesis primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (sequence not disclosed)")
    second_pol=g.input("second_polymerase","dUTP second-strand polymerase",MaterialKind.CATALYST)
    duplex=g.product("duplex","dUTP-marked cDNA duplex",MaterialKind.MOLECULAR,state=END_PREP.insert)
    g.react("second_strand","Synthesize a dUTP-marked second strand",StepKind.EXTENSION,((first,InputRole.TEMPLATE),(second_primer,InputRole.STRUCTURAL),(second_pol,InputRole.CATALYST)),(duplex,))
    endmix=g.input("endprep_mix","end-repair and dA-tailing enzyme mix",MaterialKind.CATALYST)
    repaired=g.product("repaired","end-repaired cDNA duplex",MaterialKind.MOLECULAR,state=END_PREP.insert)
    g.react("end_repair","Repair cDNA duplex ends",StepKind.END_REPAIR,((duplex,InputRole.TEMPLATE),(endmix,InputRole.CATALYST)),(repaired,))
    tailed=g.product("tailed","end-repaired, dA-tailed cDNA duplex",MaterialKind.MOLECULAR,state=dA_tailed_scene(END_PREP))
    g.react("da_tail","Add one 3-prime dA to each cDNA end",StepKind.TAILING,((repaired,InputRole.TEMPLATE),(endmix,InputRole.CATALYST)),(tailed,))
    adapter=g.input("adapter","NEBNext adapter",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (adapter body varies by kit)")
    ligase=g.input("ligase","NEBNext adapter ligase",MaterialKind.CATALYST)
    adapted=g.product("adapted","adapter-ligated dUTP-marked cDNA",MaterialKind.MOLECULAR,state=final_library(),inferred=True,source=SOURCE+" (adapter body varies by kit)")
    g.react("ligate","Ligate NEBNext adapters",StepKind.LIGATION,((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(adapted,),geometry=ADAPTOR_END)
    user=g.input("user","USER enzyme",MaterialKind.CATALYST)
    selected=g.product("selected","orientation-preserving adapter-ligated first strand",MaterialKind.MOLECULAR,state=adapted.state)
    destroyed=g.product("destroyed_second","excised dUTP-labelled second strand",MaterialKind.WASTE)
    g.react("destroy_second","Excise the dUTP-labelled second strand",StepKind.DIGESTION,((adapted,InputRole.TEMPLATE),(user,InputRole.CATALYST)),(selected,destroyed))
    g.discard(destroyed,"cannot seed directional library PCR after USER excision")
    fwd=g.input("pcr_forward","NEBNext dual-index PCR forward primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (kit oligo sequence not disclosed)")
    rev=g.input("pcr_reverse","NEBNext dual-index PCR reverse primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (kit oligo sequence not disclosed)")
    pol=g.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=g.product("library","directional dual-index RNA library",MaterialKind.LIBRARY,state=final_library())
    g.react("pcr","PCR-complete the surviving strand",StepKind.PCR,((selected,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),(rev,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(final_library(),reason="The product geometry is known but the manual does not disclose the complete kit primer sequences."))
    g.sequence(library,SequencingPlan(library=final_library(),primers=tuple(SEQ_PRIMERS),acquired_roles=REQUIRED_ROLES))
    g.declare_reviewed_scope(["rna","rt_primer","rt","reverse_transcribe","second_primer","second_polymerase","second_strand","endprep_mix","end_repair","da_tail","adapter","ligase","ligate","user","destroy_second","pcr_forward","pcr_reverse","pcr_polymerase","pcr","library"],SOURCE)
    return g.require_publishable()
