"""HiRES joint chromatin-contact and RNA profiling."""
from chemdraw import Row
from multibranch import illumina_branch
from protocol_extensions import analyte_split, proximity_product
TITLE="HiRES — single-cell Hi-C and RNA-seq"
NOTES="01_hires.html"
SOURCE='Defining source: <a href="https://doi.org/10.1126/science.adg3797">Liu et al., <i>Science</i> (2023)</a>.'
SUMMARY="Fixed nuclei undergo in situ proximity ligation and RNA capture with coordinated combinatorial indexes, producing separately amplified contact and transcript libraries from the same cell."
CELL_PARTS=(("round 1 index",8),("round 2 index",8),("round 3 index",8))
CONTACT,CP=illumina_branch("HiRES contacts","DNA–DNA ligation junction",cell_parts=CELL_PARTS)
RNA,RP=illumina_branch("HiRES RNA","RNA-derived cDNA",cell_parts=CELL_PARTS,umi=8)
FINAL_LIBRARIES=(("Chromatin-contact library",CONTACT,CP,"Paired genomic inserts flank an in situ ligation junction.","Illumina primers on the contact branch."),("RNA library",RNA,RP,"Transcript cDNA carries the matching combinatorial cell identity.","Illumina primers on the RNA branch."))
def sections():
    return [("Digest fixed chromatin and ligate nearby DNA ends",proximity_product(left="contact locus A",right="contact locus B",bridge="Hi-C junction").rows(),"In situ ligation records spatially proximal genomic loci."),("Capture and reverse-transcribe RNA in the same nucleus",[Row(chunks=[("cellular RNA → indexed cDNA", "r1", False)])],"The RNA route remains distinct from the DNA–DNA contact junction."),("Apply coordinated combinatorial indexes",[Row(chunks=[("Hi-C product ** cell-index rounds", "cbc", False)]),Row(chunks=[("RNA cDNA ** matching cell-index rounds", "cbc", False)])],"The index combination links the two modalities to one cell."),("Amplify branch-specific libraries",analyte_split("chromatin-contact library","RNA library"),"Separate primer selections retain the shared cellular provenance."),]
