"""UDA-seq universal droplet-first combinatorial indexing."""
from multibranch import illumina_branch
from protocol_extensions import analyte_split, overload_then_unpack
TITLE="UDA-seq — droplet-first combinatorial indexing"
NOTES="01_uda-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41592-024-02586-y">Wang et al., <i>Nature Methods</i> (2025)</a>.'
SUMMARY="Fixed cells or nuclei are intentionally overloaded into droplets for the first barcode, recovered intact, redistributed to PCR wells and given a second index; the barcode pair resolves cells for several droplet-compatible modalities."
CELL_PARTS=(("droplet barcode",16),("well barcode",8))
RNA,RP=illumina_branch("UDA-seq RNA","RNA-derived cDNA",cell_parts=CELL_PARTS,umi=10)
ATAC,AP=illumina_branch("UDA-seq ATAC","accessible genomic DNA",architecture="nextera",cell_parts=CELL_PARTS)
FINAL_LIBRARIES=(("Representative RNA library",RNA,RP,"Droplet barcode × well barcode identifies each recovered cell.","Illumina primers on an RNA UDA-seq branch."),("Representative ATAC library",ATAC,AP,"The same indexing logic can wrap a Multiome ATAC branch.","Nextera primers on an ATAC UDA-seq branch."))
def sections():
    return [("Fix cells or nuclei and prepare the chosen assay",overload_then_unpack(droplet_barcode=16,well_barcode=8),"The molecular target may be RNA, ATAC, VDJ or guide capture; UDA changes the indexing wrapper."),("Overload the microfluidic droplets",overload_then_unpack(droplet_barcode=16,well_barcode=8)[:1],"Multiple intact cells share one droplet barcode by design."),("Break the emulsion and recover intact cells",overload_then_unpack(droplet_barcode=16,well_barcode=8)[1:],"Recovered cells are randomly redistributed among 96 or 384 wells for index PCR."),("Build modality-specific sublibraries",analyte_split("RNA / VDJ / guide branch as configured","DNA / ATAC branch as configured"),"The final construct depends on the wrapped assay; the two representative branches show the published Multiome use."),]
