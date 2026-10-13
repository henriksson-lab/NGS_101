from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import joined_scene,seg,truseq_library
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct,MolecularState,Row,Segment,Workflow,circle_rows
from circular import circularize_ssdna,rolling_circle
from pcr import PCRGeometry
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind
from protocol_model import molecule,role

TITLE="eccDNA Circle-Seq — enrichment of native DNA circles"
NOTES="01_eccdna-circle-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1073/pnas.1508825112">Møller et al., <i>PNAS</i> (2015)</a>.'
SUMMARY="Native extrachromosomal circles survive exhaustive restriction/exonuclease depletion of linear DNA, are amplified by phi29 rolling-circle amplification, then fragmented for indexed Illumina sequencing."
CAVEAT="The protocol enriches pre-existing eccDNA; it does not circularize genomic DNA during library preparation. The Apollo 324 adapter bases were not printed, so standard Illumina arms are shown as inferred chemistry."
UNIT=Construct([Segment("native eccDNA unit","ACGTTGCAACGATCGT",placeholder=False)],name="representative eccDNA")
CIRCLE=circularize_ssdna(UNIT,five_prime_phosphate=True)
RCA=rolling_circle(CIRCLE,"TGCAACGT",copies=3)
FINAL_LIBRARY,SEQ_PRIMERS=truseq_library([seg("phi29-amplified eccDNA fragment","X"*48,placeholder=True)],"eccDNA Circle-Seq library",dual_index=False,inferred_adapters=True)
FINAL_CAPTION="Sonicated phi29 product receives the Apollo 324 platform's indexed Illumina adapters; their unpublished bases are marked inferred."
SEQUENCING_INTRO="The study used single-end HiSeq 2000 sequencing. The complete inferred library geometry also shows the opposite standard primer site without implying that Read 2 was collected."

def sections(): return [
 ("Purify circular DNA from cells",circle_rows(UNIT,"native eccDNA junction"),"The molecular input is already circular in vivo; column purification enriches the circular fraction."),
 ("Cut linear chromosomes and digest every free-ended molecule",[Row(chunks=[("native circle survives  |  NotI-cut linear DNA -- Plasmid-Safe DNase (130 h) -> degraded",None,False)])],"NotI increases accessible ends and repeated ATP/DNase additions deplete linear DNA by many orders of magnitude."),
 ("Amplify surviving circles with phi29",[Row(chunks=[(f"rolling-circle concatemer: {RCA.sequence}","r1",False)])],"Random-primed strand-displacing synthesis makes tandem copies of each surviving circular template."),
 ("Sonicate the concatemer and attach indexed adapters",AdapterLigationGeometry.unavailable("Møller et al. identify Apollo 324 library preparation but do not print its adapter oligos or ligation-facing overhang geometry.").rows()+joined_scene([seg("indexed Illumina left arm","X"*20,"r1",placeholder=True,inferred=True),seg("sonicated phi29 concatemer fragment","X"*40,placeholder=True),seg("indexed Illumina right arm","X"*20,"r2",placeholder=True,inferred=True)],(("indexed Illumina left arm","sonicated phi29 concatemer fragment","adapter ligation"),("sonicated phi29 concatemer fragment","indexed Illumina right arm","adapter ligation")),label="INFERRED adapter-ligated eccDNA fragment").rows(),"The amplified product is made into a conventional single-end sequencing library."),
]

def workflow():
 s=sections();w=Workflow(MolecularState(s[0][0],tuple(s[0][1])));w.react(s[1][0],s[1][1],note=s[1][2]);w.react(s[2][0],s[2][1],note=s[2][2])
 g=AdapterLigationGeometry.unavailable("Møller et al. identify Apollo 324 library preparation but do not print its adapter oligos or ligation-facing overhang geometry.")
 w.ligate_adapter(s[3][0],g,s[3][1],note=s[3][2]);return w

def protocol_graph():
 graph=ProtocolGraph(TITLE,authority=SOURCE)
 dna=graph.input("cellular_dna","cellular DNA containing native eccDNA circles and linear chromosomes",MaterialKind.BIOLOGICAL,state=molecule("cellular DNA",(role("native eccDNA circles"),role("linear chromosomes"))),primary=True)
 purification=graph.input("circle_purification","circular-DNA purification condition",MaterialKind.CONDITION)
 fraction=graph.product("circle_enriched_fraction","circle-enriched cellular DNA fraction",MaterialKind.MOLECULAR,state=dna.state)
 graph.react("purify_circles","Purify the circular-DNA-enriched fraction",StepKind.EXTRACTION,((dna,InputRole.TEMPLATE),(purification,InputRole.CONDITION)),(fraction,))
 noti=graph.input("noti","NotI restriction endonuclease",MaterialKind.CATALYST)
 cut=graph.product("noti_treated","NotI-treated circle-enriched DNA",MaterialKind.MOLECULAR,state=fraction.state)
 graph.react("cut_linear","Cut susceptible linear chromosomes with NotI",StepKind.DIGESTION,((fraction,InputRole.TEMPLATE),(noti,InputRole.CATALYST)),(cut,))
 exonuclease=graph.input("plasmid_safe_dnase","Plasmid-Safe ATP-dependent DNase",MaterialKind.CATALYST)
 circles=graph.product("native_circles","exonuclease-resistant native eccDNA circles",MaterialKind.MOLECULAR,state=CIRCLE)
 linear_waste=graph.product("linear_dna_waste","degraded free-ended chromosomal DNA",MaterialKind.WASTE)
 graph.react("digest_linear","Exhaustively digest free-ended DNA",StepKind.DIGESTION,((cut,InputRole.TEMPLATE),(exonuclease,InputRole.CATALYST)),(circles,linear_waste))
 graph.discard(linear_waste,"has free ends and is depleted during native-circle enrichment")
 phi29=graph.input("phi29","phi29 strand-displacing polymerase",MaterialKind.CATALYST)
 amplified=graph.product("rca_product","rolling-circle-amplified eccDNA concatemers",MaterialKind.MOLECULAR,state=RCA)
 graph.react("rolling_circle","Amplify surviving native circles with phi29",StepKind.EXTENSION,((circles,InputRole.TEMPLATE),(phi29,InputRole.CATALYST)),(amplified,))
 sonication=graph.input("sonication","concatemer-sonication condition",MaterialKind.CONDITION)
 fragments=graph.product("fragments","sonicated phi29 concatemer fragments",MaterialKind.MOLECULAR,state=amplified.state)
 graph.react("fragment","Sonicate the rolling-circle products",StepKind.FRAGMENTATION,((amplified,InputRole.TEMPLATE),(sonication,InputRole.CONDITION)),(fragments,))
 adapter=graph.input("apollo_adapter","Apollo 324 indexed Illumina adapter",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (adapter oligos not printed)")
 ligase=graph.input("adapter_ligase","library-adapter DNA ligase",MaterialKind.CATALYST)
 ligated=graph.product("adapter_ligated","adapter-ligated eccDNA fragments",MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True,source=SOURCE+" (adapter bases unavailable)")
 graph.react("adapter_ligation","Attach Apollo 324 indexed adapters",StepKind.LIGATION,((fragments,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(ligated,),geometry=AdapterLigationGeometry.unavailable("Møller et al. identify Apollo 324 library preparation but do not print its adapter oligos or ligation-facing geometry."))
 forward=graph.input("pcr_forward","Apollo library PCR forward primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (primer sequence not printed)")
 reverse=graph.input("pcr_reverse","Apollo library PCR reverse primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (primer sequence not printed)")
 polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
 library=graph.product("sequencing_library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY,inferred=True,source=SOURCE+" (compatible Illumina shell inferred)")
 graph.react("library_pcr","Amplify the indexed eccDNA library",StepKind.PCR,((ligated,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),(reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(FINAL_LIBRARY,reason="The Apollo 324 PCR primers are not printed in the defining source."))
 graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),acquired_roles=("Read 1",)))
 graph.declare_reviewed_scope(("cellular_dna","circle_purification","purify_circles","noti","cut_linear","plasmid_safe_dnase","digest_linear","phi29","rolling_circle","sonication","fragment","apollo_adapter","adapter_ligase","adapter_ligation","pcr_forward","pcr_reverse","pcr_polymerase","library_pcr","sequencing_library"),SOURCE)
 return graph.require_publishable()
