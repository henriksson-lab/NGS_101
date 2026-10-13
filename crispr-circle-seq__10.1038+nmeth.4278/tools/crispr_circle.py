from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import joined_scene,seg,truseq_library
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct,MolecularState,Row,Scene,Segment,Workflow,circle_rows
from circular import CircularizationGeometry,circularize_ssdna
from pcr import PCRGeometry
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind

TITLE="CRISPR CIRCLE-seq — in-vitro off-target cleavage"
NOTES="01_crispr-circle-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/nmeth.4278">Tsai et al., <i>Nature Methods</i> (2017)</a>.'
SUMMARY="Random genomic fragments are converted to covalently closed circles and residual linear DNA is destroyed; Cas9–guide cleavage alone re-creates ends that can accept sequencing adapters."
CAVEAT="This CIRCLE-seq profiles CRISPR cleavage in vitro. It is unrelated to eccDNA Circle-Seq, which enriches naturally circular DNA and amplifies it with phi29."
STEM_LOOP="CGGTGGACCGATGATCUATCGGTCCACCGT"
LINEAR=Construct([Segment("genomic fragment","X"*48,placeholder=True),Segment("USER-opened stem-loop scar","N"*8,placeholder=True)],name="circularization substrate")
CIRCLE=circularize_ssdna(LINEAR,five_prime_phosphate=True)
FINAL_LIBRARY,SEQ_PRIMERS=truseq_library([seg("Cas9-cleaved circle flank A","X"*35,placeholder=True),seg("cut site","N"*6,placeholder=True),seg("Cas9-cleaved circle flank B","X"*35,placeholder=True)],"CRISPR CIRCLE-seq library")
FINAL_CAPTION="Cas9 linearization places both sides of one cut site in the same paired-end library molecule."
SEQUENCING_INTRO="Paired-end 150-bp reads approach the single in-vitro Cas9 cleavage site from its two newly adapter-ligated ends."

def stem_loop_geometry():
 sc=Scene();sc.strand("prepared genomic DNA",[seg("genomic duplex","X"*28,placeholder=True),seg("DNA 3′ dA","A")],label="dA-tailed fragment")
 sc.overhang("prepared genomic DNA","DNA 3′ dA",polarity="3-prime",owner="substrate",text="substrate 3′-A")
 sc.strand("oSQT1288",[seg("USER stem-loop body",STEM_LOOP[:-1],"r1"),seg("stem-loop 3′ dT",STEM_LOOP[-1],"r1")],label="5′-phosphorylated oSQT1288")
 sc.overhang("oSQT1288","stem-loop 3′ dT",polarity="3-prime",owner="adapter",text="adapter 3′-T")
 return AdapterLigationGeometry.known("cohesive",sc,statement="The supplementary oSQT1288 sequence ends in 3′-T, complementary to the explicitly dA-tailed genomic substrate.")

def sections(): return [
 ("Shear genomic DNA, repair and add 3-prime dA ends",Scene.duplex([seg("3-prime dA","A"),seg("genomic fragment","X"*40,placeholder=True),seg("3-prime dA","A")],label="dA-tailed genomic fragment").rows(),"The cohesive ends are shown before stem-loop ligation."),
 ("Ligate the USER stem-loop",joined_scene([seg("genomic fragment","X"*40,placeholder=True),seg("oSQT1288 USER stem-loop",STEM_LOOP,"r1",placeholder=True)],(("genomic fragment","oSQT1288 USER stem-loop","adapter ligation"),),label="stem-loop-ligated genomic fragment",duplex=False).rows(),"The exact oSQT1288 oligo carries a 5-prime phosphate, one deoxyuridine and a terminal phosphorothioate bond."),
 ("Open the stem-loop and circularize each fragment",circle_rows(LINEAR,"intramolecular T4 ligase closure"),"USER and T4 PNK expose compatible palindromic ends; dilute ligation favors one-fragment circles."),
 ("Destroy every remaining linear molecule",[Row(chunks=[("closed genomic circle survives  |  linear DNA -- Plasmid-Safe DNase -> degraded",None,False)])],"Exonuclease selection makes a free end a specific reporter of the next cleavage step."),
 ("Cut circles with Cas9–guide complexes",[Row(chunks=[("closed circle -- Cas9 cut at target/off-target -> one linear molecule with two fresh ends",None,False)])],"Each cleaved circle carries both flanks of a single cleavage site."),
 ("Add 3-prime dA to the Cas9-cleaved ends",Scene.duplex([seg("3-prime dA","A"),seg("Cas9-cleaved circle","X"*52,placeholder=True),seg("3-prime dA","A")],label="dA-tailed Cas9-cleaved circle").rows(),"The cleavage-created ends are shown before adapter ligation."),
 ("Ligate USER hairpin adapters, open and PCR",joined_scene([seg("left Illumina arm","X"*20,"r1",placeholder=True),seg("flank A and Cas9 cut","X"*32,placeholder=True),seg("flank B","X"*32,placeholder=True),seg("right Illumina arm","X"*20,"r2",placeholder=True)],(("left Illumina arm","flank A and Cas9 cut","adapter ligation"),("flank B","right Illumina arm","adapter ligation")),label="Cas9-cleavage library").rows(),"The cleavage-created ends are selectively converted into a paired-end library."),
]

def workflow():
 s=sections();w=Workflow(MolecularState(s[0][0],tuple(s[0][1])))
 w.ligate_adapter(s[1][0],stem_loop_geometry(),s[1][1],note=s[1][2])
 for action,rows,note in s[2:6]:w.react(action,rows,note=note)
 final_geometry=AdapterLigationGeometry.unavailable("CIRCLE-seq describes USER hairpin adapter conversion after Cas9 cleavage, but the cut-facing adapter strands have not yet been verified from the supplementary source; their overhang bases are not asserted.")
 w.ligate_adapter(s[6][0],final_geometry,s[6][1],note=s[6][2]);return w

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 dna=g.input("genomic_dna","purified genomic DNA",MaterialKind.BIOLOGICAL,primary=True)
 shear=g.input("shearing","mechanical shearing condition",MaterialKind.CONDITION)
 fragments=g.product("fragments","random genomic DNA fragments",MaterialKind.MOLECULAR)
 g.react("shear","Shear genomic DNA",StepKind.FRAGMENTATION,((dna,InputRole.TEMPLATE),(shear,InputRole.CONDITION)),(fragments,))
 repair=g.input("repair_mix","end-repair enzymes",MaterialKind.CATALYST)
 repaired=g.product("repaired","end-repaired genomic fragments",MaterialKind.MOLECULAR)
 g.react("end_repair","Repair genomic fragment ends",StepKind.END_REPAIR,((fragments,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(repaired,))
 tail=g.input("tailing_polymerase","dA-tailing polymerase",MaterialKind.CATALYST)
 datp=g.input("datp","dATP tailing nucleotide",MaterialKind.MOLECULAR)
 tailed=g.product("da_tailed","3′-dA-tailed genomic fragments",MaterialKind.MOLECULAR,state=sections()[0][1])
 g.react("da_tail","Add 3′ dA ends",StepKind.TAILING,((repaired,InputRole.TEMPLATE),(datp,InputRole.STRUCTURAL),(tail,InputRole.CATALYST)),(tailed,))
 stem=g.input("stem_loop","published oSQT1288 USER stem-loop",MaterialKind.ADAPTER)
 ligase=g.input("stem_ligase","stem-loop DNA ligase",MaterialKind.CATALYST)
 stemmed=g.product("stem_ligated","USER-stem-loop-ligated genomic fragments",MaterialKind.MOLECULAR,state=LINEAR)
 g.react("ligate_stem","Ligate the USER stem-loop",StepKind.LIGATION,((tailed,InputRole.TEMPLATE),(stem,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(stemmed,),geometry=stem_loop_geometry())
 user=g.input("user_enzyme","USER enzyme",MaterialKind.CATALYST)
 opened=g.product("opened_stem","USER-opened circularization substrate",MaterialKind.MOLECULAR,state=LINEAR)
 g.react("open_stem","Open the uracil-containing stem-loop",StepKind.DIGESTION,((stemmed,InputRole.TEMPLATE),(user,InputRole.CATALYST)),(opened,))
 pnkligase=g.input("circle_ligase","T4 PNK and circle-closing DNA ligase",MaterialKind.CATALYST)
 circle=g.product("circles","covalently closed genomic circles",MaterialKind.MOLECULAR,state=CIRCLE)
 g.react("circularize","Circularize each prepared genomic fragment",StepKind.CIRCULARIZATION,((opened,InputRole.TEMPLATE),(pnkligase,InputRole.CATALYST)),(circle,),geometry=CircularizationGeometry.known("double-strand intramolecular ligation",statement="USER-opened, PNK-prepared palindromic ends are dilute-ligated into one-fragment circles."))
 exonuclease=g.input("exonuclease","Plasmid-Safe linear-DNA exonuclease",MaterialKind.CATALYST)
 selected=g.product("selected_circles","exonuclease-resistant genomic circles",MaterialKind.MOLECULAR,state=CIRCLE)
 linear_waste=g.product("linear_waste","digested residual linear DNA",MaterialKind.WASTE)
 g.react("digest_linear","Destroy every remaining linear molecule",StepKind.DIGESTION,((circle,InputRole.TEMPLATE),(exonuclease,InputRole.CATALYST)),(selected,linear_waste));g.discard(linear_waste,"is not covalently closed")
 guide=g.input("guide_rna","target or off-target interrogating guide RNA",MaterialKind.OLIGO)
 cas9=g.input("cas9","Cas9 nuclease",MaterialKind.CATALYST)
 cut=g.product("cas9_linearized","Cas9-linearized circle with two fresh ends",MaterialKind.MOLECULAR)
 g.react("cas9_cut","Cut circles with Cas9–guide complexes",StepKind.DIGESTION,((selected,InputRole.TEMPLATE),(guide,InputRole.STRUCTURAL),(cas9,InputRole.CATALYST)),(cut,))
 tail2=g.input("post_cut_tailing","post-cleavage dA-tailing polymerase",MaterialKind.CATALYST)
 datp2=g.input("post_cut_datp","post-cleavage dATP",MaterialKind.MOLECULAR)
 cut_tailed=g.product("cut_da_tailed","dA-tailed Cas9-cleaved circles",MaterialKind.MOLECULAR,state=sections()[5][1])
 g.react("tail_cut","Add 3′ dA to Cas9-created ends",StepKind.TAILING,((cut,InputRole.TEMPLATE),(datp2,InputRole.STRUCTURAL),(tail2,InputRole.CATALYST)),(cut_tailed,))
 adapters=g.input("library_adapters","USER hairpin sequencing adapters",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (cut-facing strands not verified)")
 lib_ligase=g.input("library_ligase","library-adapter DNA ligase",MaterialKind.CATALYST)
 ligated=g.product("adapter_ligated","adapter-ligated cleavage products",MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True,source=SOURCE+" (adapter geometry unavailable)")
 g.react("ligate_library","Ligate USER hairpin sequencing adapters",StepKind.LIGATION,((cut_tailed,InputRole.TEMPLATE),(adapters,InputRole.STRUCTURAL),(lib_ligase,InputRole.CATALYST)),(ligated,),geometry=AdapterLigationGeometry.unavailable("The cut-facing USER hairpin adapter strands have not been verified from the supplementary source."))
 f=g.input("pcr_forward","library PCR forward primer",MaterialKind.OLIGO);r=g.input("pcr_reverse","library PCR reverse primer",MaterialKind.OLIGO);pol=g.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
 library=g.product("sequencing_library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 g.react("library_pcr","Open and PCR-amplify the paired-end library",StepKind.PCR,((ligated,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(FINAL_LIBRARY,reason="The verified model does not yet contain the cut-facing USER hairpin primer geometry."))
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),acquired_roles=tuple(p.role for p in SEQ_PRIMERS)))
 reviewed=("genomic_dna","shearing","shear","repair_mix","end_repair","tailing_polymerase","datp","da_tail","stem_loop","stem_ligase","ligate_stem","user_enzyme","open_stem","circle_ligase","circularize","exonuclease","digest_linear","guide_rna","cas9","cas9_cut","post_cut_tailing","post_cut_datp","tail_cut","library_adapters","library_ligase","ligate_library","pcr_forward","pcr_reverse","pcr_polymerase","library_pcr","sequencing_library")
 g.declare_reviewed_scope(reviewed,SOURCE);return g.require_publishable()
