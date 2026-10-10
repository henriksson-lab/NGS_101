"""scDam&T-seq DamID and transcriptome branches."""
from chemdraw import Row
from multibranch import illumina_branch
from protocol_extensions import analyte_split
TITLE="scDam&T-seq — protein–DNA contacts and transcriptome"
NOTES="01_scdam-t-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41587-019-0150-y">Rooijers et al., <i>Nature Biotechnology</i> (2019)</a>; detailed protocol <a href="https://doi.org/10.1038/s41596-020-0314-8">Nat. Protocols (2020)</a>.'
SUMMARY="Dam fusion proteins methylate nearby GATC sites in vivo; mRNA is removed on oligo(dT) beads, leaving genomic DNA for DpnI-selective DamID amplification."
RNA,RP=illumina_branch("scDam&T-seq RNA","CEL-seq2 transcript insert",umi=6)
DNA,DP=illumina_branch("scDam&T-seq DamID","methylated-GATC genomic insert")
FINAL_LIBRARIES=(("Transcriptome library",RNA,RP,"Bead-bound RNA is processed by CEL-seq2.","Illumina primers on the RNA branch."),("DamID library",DNA,DP,"DpnI-selective fragments report protein-proximal adenine methylation.","Illumina primers on the DamID branch."))
def sections():
    return [("Mark protein-proximal GATC sites in vivo",[Row(chunks=[("DNA — G^mATC — DNA   ← Dam fusion at the chromatin target", "w1", False)])],"Dam deposition is the molecular record of the protein–DNA contact."),("Capture poly(A) RNA on magnetic beads",analyte_split("bead-bound RNA → CEL-seq2","supernatant DNA → DamID"),"The plate position preserves cell identity across branches."),("Select methylated genomic GATC",[Row(chunks=[("DpnI cuts G^mATC → adapter ligation ** → selective PCR", "me", False)])],"DpnI creates amplifiable ends only at methylated GATC sites."),]
