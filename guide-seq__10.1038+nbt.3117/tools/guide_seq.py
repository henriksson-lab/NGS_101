from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import joined_scene,seg,truseq_library
from chemdraw import Row,Scene,Segment,feature

TITLE="GUIDE-seq — cellular capture of nuclease-cut sites"
NOTES="01_guide-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/nbt.3117">Tsai et al., <i>Nature Biotechnology</i> (2015)</a>.'
SUMMARY="A protected, phosphorylated blunt dsODN integrates into nuclease-created breaks in living cells; single-tail adapter/tag PCR then recovers genomic DNA immediately beside either dsODN end."
CAVEAT="The two STAT-PCR orientations are parallel libraries from the two dsODN strands. The schematic shows one orientation and names the mirrored reaction explicitly."
UMI=feature("guide_seq_umi","umi","random")
FINAL_LIBRARY,SEQ_PRIMERS=truseq_library([seg("8-nt molecular barcode","U"*8,"umi",placeholder=True,feature=UMI),seg("integrated dsODN end","D"*18,placeholder=True),seg("adjacent genomic DNA","X"*42,placeholder=True)],"GUIDE-seq STAT-PCR library")
FINAL_CAPTION="STAT-PCR joins one dsODN-specific end to a single-tail adapter, while an 8-nt molecular barcode distinguishes captured molecules."
SEQUENCING_INTRO="Read 1 enters the dsODN/genome junction; the opposite read and index sites are supplied by nested library PCR."

def integration_scene():
    left=[Segment("genome left","X"*24,placeholder=True),Segment("dsODN","D"*34,"cbc",placeholder=True),Segment("genome right","X"*24,placeholder=True)]
    sc=Scene.duplex(left,label="NHEJ-tagged break")
    sc.junction("top","genome left","dsODN","NHEJ capture")
    sc.junction("top","dsODN","genome right","NHEJ capture")
    sc.labels("top"); return sc

def sections(): return [
 ("Integrate the protected blunt dsODN at cellular DSBs",integration_scene().rows(),"The 34-bp dsODN is 5-prime phosphorylated and terminally phosphorothioate-protected; NHEJ captures it at on- and off-target breaks."),
 ("Shear genomic DNA and ligate a single-tail adapter",joined_scene([seg("single-tail adapter","X"*20,"r1",placeholder=True),seg("random genomic shear fragment","X"*34,placeholder=True),seg("integrated dsODN","D"*18,"cbc",placeholder=True)],(("single-tail adapter","random genomic shear fragment","adapter ligation"),),label="single-tail capture molecule").rows(),"The adapter has only one amplifiable tail, suppressing background from fragments lacking the dsODN tag."),
 ("Run dsODN-to-adapter STAT-PCR in both orientations",[Row(chunks=[("dsODN strand A primer ---> flank A    |    dsODN strand B primer ---> flank B", "r1",False)])],"Separate reactions recover genomic sequence from both sides of the integrated tag."),
 ("Nested PCR adds UMI, index and complete Illumina arms",[Row(chunks=[("P5 -- UMI -- dsODN -- genomic flank -- P7",None,False)])],"The random 8-bp molecular barcode supports PCR-bias correction before sequencing."),
]
