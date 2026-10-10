"""scM&T-seq physical RNA/DNA separation."""
from chemdraw import Row
from multibranch import illumina_branch
from protocol_extensions import analyte_split
TITLE="scM&T-seq — transcriptome and DNA methylome"
NOTES="01_scm-t-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/nmeth.3728">Angermueller et al., <i>Nature Methods</i> (2016)</a>.'
SUMMARY="Biotinylated oligo(dT) beads remove a cell's mRNA for Smart-seq2 while the genomic-DNA supernatant enters scBS-seq."
RNA,RP=illumina_branch("scM&T-seq RNA","Smart-seq2 cDNA",architecture="nextera")
DNA,DP=illumina_branch("scM&T-seq DNA","bisulfite-converted genomic DNA")
FINAL_LIBRARIES=(("RNA library",RNA,RP,"Smart-seq2 cDNA after Nextera XT.","Nextera primers declared for the RNA branch."),("DNA methylome library",DNA,DP,"Bisulfite-converted genomic fragments after scBS-seq.","Illumina primers declared for the DNA branch."))
def sections():
    return [("Capture mRNA on biotinylated oligo(dT) beads",[Row(chunks=[("mRNA poly(A) || biotin—oligo(dT) ● bead", "r1", False)])],"The bead-bound fraction is physically removed from genomic DNA."),("Separate the analytes",analyte_split("bead-bound RNA → Smart-seq2","supernatant DNA → scBS-seq"),"Wash fractions are retained with the DNA to reduce loss."),("Convert and amplify each branch",[Row(chunks=[("RNA: RT → template switch → cDNA PCR → Nextera", "tso", False)]),Row(chunks=[("DNA: bisulfite → random priming → adapter PCR", "r2", False)])],"The branch chemistries remain distinct after separation."),]
