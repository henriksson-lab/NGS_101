"""OAK overloading-and-unpacking combinatorial indexing."""
from multibranch import illumina_branch
from protocol_extensions import analyte_split, overload_then_unpack
from chemdraw import Row
TITLE="OAK — overloading and unpacking droplet barcodes"
NOTES="01_oak.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-024-53227-z">Darmanis et al., <i>Nature Communications</i> (2024)</a>.'
SUMMARY="Fixed cells or nuclei are overloaded into Chromium droplets for in situ first-round barcoding, recovered intact and redistributed into aliquots for a second index; RNA, ATAC and antibody-tag branches can share the resulting cell identity."
CELL_PARTS=(("droplet barcode",16),("aliquot index",8))
RNA,RP=illumina_branch("OAK RNA","3-prime cDNA",cell_parts=CELL_PARTS,umi=10)
ATAC,AP=illumina_branch("OAK ATAC","accessible genomic DNA",architecture="nextera",cell_parts=CELL_PARTS)
FINAL_LIBRARIES=(("RNA library",RNA,RP,"Droplet and aliquot indexes jointly identify a cell's transcripts.","Illumina primers on the OAK RNA branch."),("ATAC library",ATAC,AP,"The matching two-round identity accompanies Multiome ATAC fragments.","Nextera primers on the OAK ATAC branch."))
def sections():
    return [("Fix cells or nuclei so they survive droplet recovery",[Row(chunks=[("fixed intact cell/nucleus containing retained indexed molecules", None, False)])],"The cell itself becomes the carrier between indexing rounds."),("Overload Chromium droplets and perform first-round barcoding",overload_then_unpack(droplet_barcode=16,well_barcode=8)[:1],"Many cells use each droplet's bead barcode."),("Break droplets, mix cells and redistribute into aliquots",overload_then_unpack(droplet_barcode=16,well_barcode=8)[1:],"Aliquot-specific PCR supplies the orthogonal secondary index."),("Separate configured modalities",analyte_split("RNA / antibody-tag library","ATAC library when Multiome is used"),"The demonstration supports RNA alone and joint RNA–ATAC; the cell identity is shared."),]
