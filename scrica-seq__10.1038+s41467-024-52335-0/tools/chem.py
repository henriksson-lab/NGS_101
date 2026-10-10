"""scRICA-seq circular cDNA isoform reconstruction plus ATAC."""
from multibranch import illumina_branch
from protocol_extensions import analyte_split, indexed_molecule
from chemdraw import Row
TITLE="scRICA-seq — RNA isoforms and chromatin accessibility"
NOTES="01_scrica-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-024-52335-0">Jiang et al., <i>Nature Communications</i> (2024)</a>.'
SUMMARY="Full-length cDNA receives end UMIs, is circularized and amplified so short Tn5 fragments can be grouped back into isoforms; a matched ATAC branch reports chromatin accessibility from the same cell."
RNA,RP=illumina_branch("scRICA-seq isoform","RCA-derived full-length-cDNA fragment",architecture="nextera",cell=16,umi=12)
ATAC,AP=illumina_branch("scRICA-seq ATAC","accessible genomic DNA",architecture="nextera",cell=16)
FINAL_LIBRARIES=(("Isoform library",RNA,RP,"Tagmented RCA copies share the UMI of one original full-length cDNA.","Nextera primers on the isoform branch."),("Accessibility library",ATAC,AP,"ATAC fragments carry the matching cell barcode.","Nextera primers on the accessibility branch."))
def sections():
    return [("Capture RNA and accessible DNA from the same cell",analyte_split("full-length RNA/cDNA branch","ATAC branch"),"Low-throughput separation or the high-throughput microfluidic implementation preserves cellular pairing."),("Label both cDNA ends with molecule identity",indexed_molecule(payload="full-length cDNA",barcode_parts=(("cell barcode",16),),umi=12).rows(),"End labels let fragments from one original transcript be regrouped."),("Circularize and amplify full-length cDNA",[Row(chunks=[("end-labelled cDNA ** end-to-end ligation → circle → rolling-circle copies", "w1", False)])],"Circular amplification makes repeated copies of the same full-length molecule."),("Tagment the amplified circle and build both libraries",[Row(chunks=[("RCA concatemer → random Tn5 fragments → isoform library", "me", False)]),Row(chunks=[("accessible genomic DNA → ATAC library", "me", False)])],"Short isoform fragments sharing a UMI are integrated computationally across TSS, TES and exons."),]
