"""iMARGI in-situ RNA–DNA bridge ligation."""
from multibranch import illumina_branch
from protocol_extensions import proximity_product
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState,Row,Scene,Segment,Workflow
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)
TITLE="iMARGI — in situ RNA–genome interaction mapping"
NOTES="01_imargi.html"
SOURCE='Defining source: <a href="https://doi.org/10.1073/pnas.1819788116">Yan et al., <i>PNAS</i> (2019)</a>; detailed protocol <a href="https://doi.org/10.1038/s41596-019-0229-4">Nat. Protocols (2019)</a>.'
SUMMARY="A bivalent linker is ligated in situ first to fragmented chromatin-associated RNA and then to nearby genomic DNA; reverse transcription and Illumina preparation preserve the RNA–linker–DNA order."
FINAL_LIBRARY,SEQ_PRIMERS=illumina_branch("iMARGI","RNA–linker–DNA contact chimera",umi=8)
FINAL_CAPTION="Paired reads recover the RNA-derived and DNA-derived ends separated by the diagnostic bridge."
SEQUENCING_INTRO="Declared Illumina primers bind the completed iMARGI library."
def rna_ligation_geometry():
    rna=Scene(); rna.strand("RNA arm",[Segment("preadenylated RNA-ligation end","[RNA-arm end]",placeholder=True),Segment("annealed bivalent-linker body","[bivalent linker]",placeholder=True)],label="iMARGI linker")
    rna.ligation_end("RNA arm","preadenylated RNA-ligation end",kind="single-strand",polarity="5-prime",owner="adapter",text="preadenylated 5′ RNA-ligation end")
    return AdapterLigationGeometry.known("single-strand",rna,statement="The detailed protocol identifies a preadenylated RNA-compatible end on the annealed bivalent linker; this is a single-strand end ligation, not a duplex overhang.")
def dna_ligation_geometry():
    return AdapterLigationGeometry.unavailable("The sources establish the second linker-to-restriction-fragment ligation but do not provide a source-supported adapter-overhang sequence or polarity.")
def sections():
    wf=Workflow(MolecularState("Crosslinked and fragmented chromatin",(Row(chunks=[("RNA … crosslinked proximity … restriction-cut genomic DNA",None,False)]),)))
    wf.ligate_adapter("Ligate the RNA end of the bivalent linker in situ",rna_ligation_geometry(),[Row(chunks=[("chromatin RNA ** RNA arm—annealed bridge","me",False)])],note="The linker end joins fragmented RNA.")
    wf.ligate_adapter("Ligate the bridge DNA end to nearby genomic DNA",dna_ligation_geometry(),proximity_product(left="chromatin RNA",right="nearby genomic DNA",bridge="iMARGI bivalent linker").rows(),note="The second ligation makes the diagnostic RNA–linker–DNA chimera in the nucleus.")
    wf.react("Reverse-transcribe, purify and make the paired-end library",[Row(chunks=[("RNA-derived cDNA—linker—genomic DNA → paired-end adapters","r1",False)])],note="The linker orientation distinguishes RNA and DNA read ends.")
    return [(r.action,r.after.rows,r.note) for r in wf.reactions]

def protocol_graph():
 source="Defining source: https://doi.org/10.1073/pnas.1819788116; detailed protocol: https://doi.org/10.1038/s41596-019-0229-4"
 graph=ProtocolGraph("iMARGI",authority=source)
 chromatin=graph.input("chromatin","cross-linked RNA–chromatin contacts",
                       MaterialKind.BIOLOGICAL,primary=True)
 nuclease=graph.input("restriction_enzyme","chromatin restriction enzyme",MaterialKind.CATALYST)
 fragmented=graph.product("fragmented_contacts","fragmented RNA and restriction-cut chromatin",
                          MaterialKind.MOLECULAR)
 graph.react("fragmentation","Fragment cross-linked chromatin",StepKind.DIGESTION,
             ((chromatin,InputRole.TEMPLATE),(nuclease,InputRole.CATALYST)),(fragmented,))
 linker=graph.input("bivalent_linker","annealed iMARGI bivalent linker",MaterialKind.ADAPTER)
 rna_ligase=graph.input("rna_ligase","RNA-end ligase",MaterialKind.CATALYST)
 rna_joined=graph.product("rna_joined","RNA-end-ligated bivalent linker contacts",
                         MaterialKind.MOLECULAR)
 graph.react("rna_ligation","Ligate the RNA end of the bivalent linker in situ",
             StepKind.LIGATION,
             ((fragmented,InputRole.TEMPLATE),(linker,InputRole.STRUCTURAL),
              (rna_ligase,InputRole.CATALYST)),(rna_joined,),geometry=rna_ligation_geometry())
 dna_ligase=graph.input("dna_ligase","DNA-end ligase",MaterialKind.CATALYST)
 chimera=graph.product("contact_chimera","RNA–linker–DNA contact chimera",
                       MaterialKind.MOLECULAR,
                       state=proximity_product(left="chromatin RNA",right="nearby genomic DNA",bridge="iMARGI bivalent linker"))
 graph.react("dna_ligation","Ligate the bridge DNA end to nearby genomic DNA",
             StepKind.LIGATION,
             ((rna_joined,InputRole.TEMPLATE),(linker,InputRole.STRUCTURAL),
              (dna_ligase,InputRole.CATALYST)),(chimera,),geometry=dna_ligation_geometry())
 rt_primer=graph.input("rt_primer","iMARGI reverse-transcription primer",MaterialKind.OLIGO)
 rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
 cdna=graph.product("cdna_chimera","cDNA–linker–genomic DNA chimera",MaterialKind.MOLECULAR)
 graph.react("reverse_transcription","Reverse-transcribe the RNA side of the chimera",
             StepKind.REVERSE_TRANSCRIPTION,
             ((chimera,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
              (rt,InputRole.CATALYST)),(cdna,))
 p5=graph.input("p5_primer","P5 library PCR primer",MaterialKind.OLIGO)
 p7=graph.input("p7_primer","indexed P7 library PCR primer",MaterialKind.OLIGO)
 polymerase=graph.input("polymerase","library PCR polymerase",MaterialKind.CATALYST)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 graph.react("library_pcr","Purify and PCR-complete the paired-end library",StepKind.PCR,
             ((cdna,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
              (p7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The compact Illumina branch records the completed library rather than separate PCR landing objects."))
 graph.sequence(library,SequencingPlan(
     library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
     acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
     read_uses={p.role:("demultiplex" if "Index" in p.role else "include") for p in SEQ_PRIMERS}))
 graph.declare_reviewed_scope((
     "chromatin","restriction_enzyme","fragmentation","bivalent_linker","rna_ligase",
     "rna_ligation","dna_ligase","dna_ligation","rt_primer","reverse_transcriptase",
     "reverse_transcription","p5_primer","p7_primer","polymerase","library_pcr","library"),source)
 return graph.require_publishable()
