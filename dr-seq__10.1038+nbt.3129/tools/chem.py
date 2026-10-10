"""DR-seq co-amplification and downstream branch split."""
from chemdraw import Row
from multibranch import illumina_branch
from protocol_extensions import analyte_split, coamplified_analytes
TITLE="DR-seq — genome and transcriptome without physical separation"
NOTES="01_dr-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/nbt.3129">Dey et al., <i>Nature Biotechnology</i> (2015)</a>.'
SUMMARY="A common quasilinear amplification copies genomic DNA and poly(A)-primed cDNA together; aliquots are then converted into distinct genome and transcriptome libraries."
DNA,DP=illumina_branch("DR-seq genome","genomic insert",architecture="nextera")
RNA,RP=illumina_branch("DR-seq transcriptome","RNA-derived insert",umi=8)
FINAL_LIBRARIES=(("Genome library",DNA,DP,"Nextera library from the shared amplification product.","Nextera sequencing primers on the genome branch."),("Transcriptome library",RNA,RP,"Directional RNA-derived library from a separate aliquot.","Illumina sequencing primers on the RNA branch."))
def sections():
    return [("Reverse-transcribe poly(A) RNA in the intact lysate",[Row(chunks=[("mRNA poly(A) || tagged oligo(dT) → cDNA", "r1", False)])],"Genomic DNA remains in the same tube."),("Quasilinearly amplify DNA and cDNA together",coamplified_analytes("genomic DNA","RNA-derived cDNA"),"MALBAC-derived primers give both analytes a common amplification handle."),("Split the amplified pool",analyte_split("Nextera genome library","directional transcript library"),"Only after shared amplification are the modalities separated into library reactions."),]
