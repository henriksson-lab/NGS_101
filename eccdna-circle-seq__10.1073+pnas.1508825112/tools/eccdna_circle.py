from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import joined_scene,seg,truseq_library
from chemdraw import Construct,Row,Segment,circle_rows
from circular import circularize_ssdna,rolling_circle

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
 ("Sonicate the concatemer and attach indexed adapters",joined_scene([seg("indexed Illumina left arm","X"*20,"r1",placeholder=True,inferred=True),seg("sonicated phi29 concatemer fragment","X"*40,placeholder=True),seg("indexed Illumina right arm","X"*20,"r2",placeholder=True,inferred=True)],(("indexed Illumina left arm","sonicated phi29 concatemer fragment","adapter ligation"),("sonicated phi29 concatemer fragment","indexed Illumina right arm","adapter ligation")),label="INFERRED adapter-ligated eccDNA fragment").rows(),"The amplified product is made into a conventional single-end sequencing library."),
]
