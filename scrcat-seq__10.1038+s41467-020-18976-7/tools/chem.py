"""scRCAT-seq cap-and-tail endpoint libraries."""
from chemdraw import Row
from multibranch import illumina_branch
TITLE="scRCAT-seq — single-cell RNA cap and tail sequencing"
NOTES="01_scrcat-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-020-18976-7">Liu et al., <i>Nature Communications</i> (2020)</a>.'
SUMMARY="Barcoded, UMI-tagged full-length cDNA is tagmented; targeted PCR separately enriches fragments adjacent to the transcript 5′ and 3′ boundaries."
CAVEAT="The RT tag, barcode lengths and selective end PCR are source-exact. The compact completed Nextera shells are marked inferred pending a full base-by-base reconstruction of both selected end products."
TSS,TP=illumina_branch("scRCAT-seq TSS library","5-prime cDNA boundary",architecture="nextera",cell=16,umi=8,inferred_adapters=True)
TES,EP=illumina_branch("scRCAT-seq TES library","3-prime cDNA boundary",architecture="nextera",cell=16,umi=8,inferred_adapters=True)
FINAL_LIBRARIES=(("TSS-end library",TSS,TP,"Tag-adjacent fragments identify the template-switched transcript 5′ end.","Nextera-compatible primers read the TSS branch."),("TES-end library",TES,EP,"Poly(A)-adjacent fragments identify the transcript 3′ end.","Nextera-compatible primers read the TES branch."))
def oligos(): return ["<b>RT primer architecture</b> — 5′-AAGCAGTGGTATCAACGCAGAGTNNNNNNNN[16-bp cell barcode]T₃₀VN-3′"]
def sections(): return [
 ("Tag both boundaries during reverse transcription",[Row(chunks=[("5′ tag—UMI—cell barcode ** full-length cDNA ** poly(A)-side tag","cbc",False)])],"The RT primer and template switch place related tags at the cDNA ends."),
 ("Amplify full-length cDNA",[Row(chunks=[("tagged cDNA → modified Smart-seq2 PCR","tso",False)])],"Plate position and the 16-base cell identifier preserve cell identity."),
 ("Tagment amplified cDNA",[Row(chunks=[("end-tagged cDNA ** internal Tn5 cuts","me",False)])],"Tagmentation creates short fragments while retaining end-adjacent molecules."),
 ("INFERRED — complete the selected end libraries",[Row(chunks=[("tag primer + Tn5 primer → TSS and TES libraries","r2",False)])],"Targeted selection is source-established; the compact outer adapter completion is inferred."),]
