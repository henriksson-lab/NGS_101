"""NEBNext Enzymatic Methyl-seq v2, E8015 manual v2.0 6/25."""
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
TITLE="NEBNext Enzymatic Methyl-seq v2"
NOTES="01_nebnext-em-seq-v2.html"
SOURCE='Commercial protocol: <a href="https://www.neb.com/en-us/products/e8015-nebnext-enzymatic-methyl-seq-v2-kit">NEB E8015 manual v2.0 (June 2025)</a>.'
SUMMARY="Adapters are ligated before enzymatic conversion. TET2 and T4-BGT protect 5mC/5hmC; APOBEC converts unmodified C to U; Q5U then amplifies the converted library."
OUTCOME={"C":"U","5mC":"C","5hmC":"C"}
DNA=Construct([Segment("fragmented DNA", "X"*42, placeholder=True)], name="fragmented DNA")
END_PREP=repair_and_dA_tail(DNA)
ADAPTOR_END=reported_adapter_end(
    "EM-seq adaptor", kind="cohesive", polarity="3-prime", terminal_sequence="T",
    substrate_scenes=(dA_tailed_scene(END_PREP),),
    statement="The E8015 workflow specifies dA-tailed DNA and its complementary 3′-dT conversion-resistant adaptor end; the remaining adaptor body is proprietary.")
PROTECTION=[Row(chunks=[("C → C     5mC → 5caC     5hmC → 5gmC",None,False)])]
DEAMINATION=[Row(chunks=[("unmodified C → U     protected 5caC / 5gmC → C in sequence",None,False)])]
def final_library(): return indexed_library(Construct([Segment("enzymatically converted insert", "X"*42, placeholder=True)]), "EM-seq v2 library")
STEPS=(("End-repair and add one 3-prime dA to each end", dA_tailed_scene(END_PREP), "The dA-tailed product is shown before the conversion-resistant adaptor is introduced."),
       ("Ligate the EM-seq adaptors", final_library(), "The conversion-resistant EM-seq adaptor is attached before base conversion.", ADAPTOR_END),
       ("Protect modified cytosines", PROTECTION, "TET2 oxidizes 5mC through 5caC; T4-BGT glucosylates 5hmC. Both products resist APOBEC."),
       ("Deaminate unmodified cytosine", DEAMINATION, "APOBEC changes unprotected C to U while protected derivatives remain readable as C."),
       ("Amplify converted strands with Q5U", Construct([Segment("converted insert: T where C was unmethylated", "X"*46, placeholder=True)], name="converted library"), "Q5U copies uracil-containing templates; dual-index primers complete the flow-cell arms."))
SEQ_PRIMERS=STANDARD_PRIMERS; REQUIRED_ROLES=tuple(p.role for p in SEQ_PRIMERS)
READOUT="Paired reads sequence the converted insert; C calls represent protected 5mC/5hmC and T calls at genomic C positions represent unmodified cytosine. i5 and i7 are separate index reads."
FINAL_CAPTION="PCR-completed, dual-index EM-seq library. The manual prints TruSeq-compatible trimming sequences, so primer binding is computed from the canonical arms."

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    dna=g.input("dna","fragmented genomic DNA",MaterialKind.BIOLOGICAL,state=DNA,primary=True)
    endmix=g.input("endprep_mix","end-repair and dA-tailing enzyme mix",MaterialKind.CATALYST)
    repaired=g.product("repaired","end-repaired fragmented DNA",MaterialKind.MOLECULAR,state=END_PREP.insert)
    g.react("end_repair","Repair fragmented DNA ends",StepKind.END_REPAIR,((dna,InputRole.TEMPLATE),(endmix,InputRole.CATALYST)),(repaired,))
    tailed=g.product("tailed","end-repaired, dA-tailed DNA",MaterialKind.MOLECULAR,state=dA_tailed_scene(END_PREP))
    g.react("da_tail","Add one 3-prime dA to each repaired end",StepKind.TAILING,((repaired,InputRole.TEMPLATE),(endmix,InputRole.CATALYST)),(tailed,))
    adapter=g.input("adapter","conversion-resistant EM-seq adapter",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (adapter body proprietary)")
    ligase=g.input("ligase","EM-seq adapter ligase",MaterialKind.CATALYST)
    adapted=g.product("adapted","adapter-ligated DNA",MaterialKind.MOLECULAR,state=final_library(),inferred=True,source=SOURCE+" (adapter body proprietary)")
    g.react("ligate","Ligate conversion-resistant adapters",StepKind.LIGATION,((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(adapted,),geometry=ADAPTOR_END)
    protect=g.input("protection_enzymes","TET2 and T4-BGT",MaterialKind.CATALYST)
    protected=g.product("protected","5mC/5hmC-protected adapter-ligated DNA",MaterialKind.MOLECULAR,state=adapted.state)
    g.react("protect","Protect modified cytosines",StepKind.CHEMICAL_CONVERSION,((adapted,InputRole.TEMPLATE),(protect,InputRole.CATALYST)),(protected,))
    apobec=g.input("apobec","APOBEC deaminase",MaterialKind.CATALYST)
    converted=g.product("converted","DNA with unmodified C converted to U",MaterialKind.MOLECULAR,state=protected.state)
    g.react("deaminate","Deaminate unmodified cytosine",StepKind.CHEMICAL_CONVERSION,((protected,InputRole.TEMPLATE),(apobec,InputRole.CATALYST)),(converted,))
    fwd=g.input("pcr_forward","NEBNext dual-index PCR forward primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (kit oligo body proprietary)")
    rev=g.input("pcr_reverse","NEBNext dual-index PCR reverse primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (kit oligo body proprietary)")
    q5u=g.input("q5u","Q5U polymerase",MaterialKind.CATALYST)
    library=g.product("library","PCR-completed EM-seq library",MaterialKind.LIBRARY,state=final_library())
    g.react("pcr","Amplify converted strands with Q5U",StepKind.PCR,((converted,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),(rev,InputRole.STRUCTURAL),(q5u,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(final_library(),reason="The manual does not disclose the complete kit primer sequences."))
    g.sequence(library,SequencingPlan(library=final_library(),primers=tuple(SEQ_PRIMERS),acquired_roles=REQUIRED_ROLES))
    g.declare_reviewed_scope(["dna","endprep_mix","end_repair","da_tail","adapter","ligase","ligate","protection_enzymes","protect","apobec","deaminate","pcr_forward","pcr_reverse","q5u","pcr","library"],SOURCE)
    return g.require_publishable()
