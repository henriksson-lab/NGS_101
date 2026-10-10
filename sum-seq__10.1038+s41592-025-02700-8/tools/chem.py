"""SUM-seq sample pre-indexing and overloaded droplet RNA/ATAC profiling."""
from multibranch import illumina_branch
from protocol_extensions import analyte_split, indexed_molecule
from chemdraw import Row
TITLE="SUM-seq — sample-indexed, overloaded single-cell RNA and ATAC"
NOTES="01_sum-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41592-025-02700-8">Llorens-Rico et al., <i>Nature Methods</i> (2025)</a>.'
SUMMARY="ATAC fragments and RNA are separately sample-indexed in fixed nuclei, then a second 10x droplet barcode resolves overloaded nuclei; the sample-index–droplet-index pair identifies each cell across both modalities."
CELL_PARTS=(("sample index",8),("droplet barcode",16))
RNA,RP=illumina_branch("SUM-seq RNA","RNA-derived cDNA",cell_parts=CELL_PARTS,umi=10)
ATAC,AP=illumina_branch("SUM-seq ATAC","accessible genomic DNA",architecture="nextera",cell_parts=CELL_PARTS)
FINAL_LIBRARIES=(("RNA library",RNA,RP,"RT sample index plus droplet barcode identifies RNA molecules.","Illumina primers on the RNA branch."),("ATAC library",ATAC,AP,"Barcoded-Tn5 sample index plus the same droplet barcode identifies ATAC fragments.","Nextera primers on the ATAC branch."))
def sections():
    return [("Install modality-specific sample indexes",[Row(chunks=[("RNA: oligo(dT) RT ** sample index", "cbc", False)]),Row(chunks=[("ATAC: accessible DNA ** barcoded Tn5 sample index", "cbc", False)])],"Each nucleus receives matching sample provenance on its two analytes."),("Tagment RNA/cDNA hybrids",[Row(chunks=[("RNA–cDNA hybrid → Tn5 primer site", "me", False)])],"This makes the RNA branch compatible with the later microfluidic barcode transfer."),("Overload nuclei into 10x droplets",indexed_molecule(payload="sample-indexed RNA or ATAC molecule",barcode_parts=(("sample index",8),("droplet barcode",16),),umi=10).rows(),"The second barcode resolves multiple nuclei sharing a droplet when combined with sample index."),("Break droplets and split the library",analyte_split("RNA-specific amplification","ATAC-specific amplification"),"Branch-specific amplification retains the shared sample × droplet identity."),]
