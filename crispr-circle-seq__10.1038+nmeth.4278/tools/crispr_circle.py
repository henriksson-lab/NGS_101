from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import joined_scene,seg,truseq_library
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct,MolecularState,Row,Scene,Segment,Workflow,circle_rows
from circular import circularize_ssdna

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
