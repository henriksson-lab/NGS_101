from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import joined_scene,seg,truseq_library
from chemdraw import Construct,Row,Segment,circle_rows
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

def sections(): return [
 ("Shear genomic DNA, repair, dA-tail and ligate the USER stem-loop",joined_scene([seg("genomic fragment","X"*40,placeholder=True),seg("oSQT1288 USER stem-loop",STEM_LOOP,"r1",placeholder=True)],(("genomic fragment","oSQT1288 USER stem-loop","adapter ligation"),),label="stem-loop-ligated genomic fragment",duplex=False).rows(),"The exact oSQT1288 oligo carries a 5-prime phosphate, one deoxyuridine and a terminal phosphorothioate bond."),
 ("Open the stem-loop and circularize each fragment",circle_rows(LINEAR,"intramolecular T4 ligase closure"),"USER and T4 PNK expose compatible palindromic ends; dilute ligation favors one-fragment circles."),
 ("Destroy every remaining linear molecule",[Row(chunks=[("closed genomic circle survives  |  linear DNA -- Plasmid-Safe DNase -> degraded",None,False)])],"Exonuclease selection makes a free end a specific reporter of the next cleavage step."),
 ("Cut circles with Cas9–guide complexes",[Row(chunks=[("closed circle -- Cas9 cut at target/off-target -> one linear molecule with two fresh ends",None,False)])],"Each cleaved circle carries both flanks of a single cleavage site."),
 ("dA-tail, ligate USER hairpin adapter, open and PCR",joined_scene([seg("left Illumina arm","X"*20,"r1",placeholder=True),seg("flank A and Cas9 cut","X"*32,placeholder=True),seg("flank B","X"*32,placeholder=True),seg("right Illumina arm","X"*20,"r2",placeholder=True)],(("left Illumina arm","flank A and Cas9 cut","adapter ligation"),("flank B","right Illumina arm","adapter ligation")),label="Cas9-cleavage library").rows(),"The cleavage-created ends are selectively converted into a paired-end library."),
]
