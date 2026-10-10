"""SNARE-seq RNA and accessible-chromatin branches."""
from batch_ngs import seg
from chemdraw import Row
from multibranch import illumina_branch
from protocol_extensions import analyte_split, indexed_molecule

TITLE="SNARE-seq — transcriptome and chromatin accessibility"
NOTES="01_snare-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41587-019-0290-0">Chen et al., <i>Nature Biotechnology</i> (2019)</a>.'
SUMMARY="A poly(dT)-bearing splint couples pre-tagmented nuclear DNA to Drop-seq beads, so RNA and accessible DNA inherit the same droplet cell barcode before separate libraries are made."

RNA, RP = illumina_branch("SNARE-seq RNA", "3-prime cDNA", cell=12, umi=8)
ATAC, AP = illumina_branch("SNARE-seq accessibility", "tagmented genomic DNA",
                           architecture="nextera", cell=12)
FINAL_LIBRARIES=(
    ("RNA library",RNA,RP,"Cell barcode and UMI precede transcript-derived sequence.","Declared Illumina primers on the RNA branch."),
    ("Accessibility library",ATAC,AP,"The matching cell barcode accompanies accessible genomic DNA.","Declared Nextera primers on the DNA branch."),)

def oligos():
    return ["<b>Nextera-R1-rc-polyA</b> — 5′-/5Phos/CTGTCTCTTATACACATCTGACGCTGCCGACGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA-3′",
            "<b>TSO</b> — 5′-AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG-3′"]

def sections():
    return [
      ("Tagment accessible chromatin in nuclei",[Row(chunks=[("accessible genomic DNA ** Nextera end", "me", False)])],"Tn5 fragments open chromatin before droplet partitioning."),
      ("Bridge tagmented DNA to poly(dT)",[Row(chunks=[("Nextera end ** phosphorylated splint — poly(A)", "me", False)])],"The splint makes the DNA product capturable by a Drop-seq bead oligo."),
      ("Co-encapsulate one nucleus and one barcode bead",indexed_molecule(payload="RNA- or DNA-derived molecule",barcode_parts=(("12-nt cell barcode",12),),umi=8).rows(),"Reverse transcription transfers bead identity; the UMI is informative for RNA."),
      ("Separate the molecular branches",analyte_split("amplify transcript-derived cDNA","amplify accessible-DNA products"),"Branch-specific PCR produces the two sequencing libraries."),]
