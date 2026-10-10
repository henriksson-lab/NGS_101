"""scTrio-seq cytoplasm/nucleus split and its RNA and DNA branches."""
from chemdraw import Row
from multibranch import illumina_branch
from protocol_extensions import analyte_split
TITLE="scTrio-seq — genome copy number, methylome and transcriptome"
NOTES="01_sctrio-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/cr.2016.23">Hou et al., <i>Cell Research</i> (2016)</a>.'
SUMMARY="A single cell is mechanically separated into cytoplasm for RNA-seq and nucleus for scRRBS; the DNA library jointly reports copy number and CpG methylation."
RNA,RP=illumina_branch("scTrio-seq RNA","cytoplasmic cDNA")
DNA,DP=illumina_branch("scTrio-seq DNA","MspI reduced-representation insert")
FINAL_LIBRARIES=(("Transcriptome library",RNA,RP,"Cytoplasmic RNA-derived library.","Illumina primers on the RNA branch."),("Genome and methylome library",DNA,DP,"One scRRBS library supports both copy-number and methylation inference.","Illumina primers on the DNA branch."))
def sections():
    return [("Separate cytoplasm from nucleus",analyte_split("cytoplasm: poly(A) RNA","nucleus: genomic DNA"),"The three measurements arise from two physical branches."),("Prepare the cytoplasmic transcriptome",[Row(chunks=[("poly(A) RNA → RT → cDNA amplification → Illumina library", "r1", False)])],"The RNA workflow follows the single-cell Tang lineage."),("Digest nuclear DNA with MspI and ligate methylated adapters",[Row(chunks=[("C^CGG genomic ends ** methylated sequencing adapter", "me", False)])],"MspI enriches CpG-containing fragments before bisulfite conversion."),("Bisulfite convert and amplify",[Row(chunks=[("scRRBS library → methylation calls + copy-number bins", "r2", False)])],"Both DNA measurements derive from this one sequencing library."),]
