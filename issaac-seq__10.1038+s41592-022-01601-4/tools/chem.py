"""ISSAAC-seq plate and droplet multimodal chemistry."""
from chemdraw import Row
from multibranch import illumina_branch
from protocol_extensions import analyte_split, indexed_molecule

TITLE="ISSAAC-seq — accessibility and RNA from the same cell"
NOTES="01_issaac-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41592-022-01601-4">Xu et al., <i>Nature Methods</i> (2022)</a>.'
SUMMARY="Biotinylated oligo(dT) captures RNA, while Tn5 tags accessible DNA; FACS or droplet indexing supplies a shared cell identity before the RNA and ATAC branches are amplified separately."
RNA,RP=illumina_branch("ISSAAC-seq RNA","3-prime cDNA",cell=16,umi=10)
ATAC,AP=illumina_branch("ISSAAC-seq ATAC","accessible genomic DNA",architecture="nextera",cell=16)
FINAL_LIBRARIES=(("RNA library",RNA,RP,"TruSeq handle, UMI and cell identity lead into transcript sequence.","Illumina primers declared for the RNA library."),("ATAC library",ATAC,AP,"Nextera-ended accessible DNA carries the matching cell identity.","Nextera primers declared for the ATAC library."))
def oligos():
    return ["<b>TruSeqR1-oligo(dT)</b> — 5′-biotin-CTACACGACGCTCTTCCGATCTNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN-3′",
            "<b>ME bottom</b> — 5′-Phos-C*T*G*T*C*T*C*T*T*A*T*A*C*A*C*A*T*C*/iInvdT/-3′"]
def sections():
    return [("Prime poly(A) RNA with the biotinylated UMI oligo",[Row(chunks=[("mRNA poly(A) || UMI—oligo(dT)—biotin", "umi", False)])],"The RT oligo identifies molecules and supplies an Illumina arm."),("Tagment accessible chromatin",[Row(chunks=[("accessible DNA ** mosaic end", "me", False)])],"Tn5 creates the ATAC branch in the same fixed cell or nucleus."),("Add cell identity",indexed_molecule(payload="RNA- or ATAC-derived product",barcode_parts=(("cell barcode",16),),umi=10).rows(),"The implementation may use FACS wells or 10x ATAC droplets; both preserve the branch identity."),("Amplify separately",analyte_split("RNA library","ATAC library"),"Branch-selective primer pairs retain the common cellular index."),]
