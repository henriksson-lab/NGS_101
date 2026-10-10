"""Slide-tags molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import nextera_library,seg,truseq_library
from chemdraw import Row,feature
from multimodal_spatial import modality_split

TITLE="Slide-tags — spatial barcodes transferred to intact nuclei"
NOTES="01_slide-tags.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41586-023-06837-4">Russell et al., <i>Nature</i> (2024)</a>.'
SUMMARY="Photocleaved oligos from a sequenced bead array diffuse into overlying tissue nuclei; 10x single-nucleus RNA or Multiome processing then links spatial tags to each nucleus barcode."
CAVEAT="The public page shows the published TAGS/SLAC transfer logic and the multiome branch. Commercial 10x bead regions remain molecular-role placeholders."
CELL=feature("slidetags_cell","cell_barcode","whitelist",whitelist="10x Chromium bead barcode")
SPATIAL=feature("slidetags_spatial","spatial_barcode","whitelist",whitelist="sequenced Slide-tags puck")
UMI=feature("slidetags_umi","umi","random")
RNA,RNA_P=truseq_library([seg("cell barcode","C"*16,"cbc",placeholder=True,feature=CELL),seg("RNA UMI","U"*12,"umi",placeholder=True,feature=UMI),seg("nuclear cDNA","X"*40,placeholder=True)],"Slide-tags RNA library")
ATAC,ATAC_P=nextera_library([seg("cell barcode association","C"*16,"cbc",placeholder=True,feature=CELL),seg("accessible chromatin","X"*40,placeholder=True)],"Slide-tags ATAC library")
TAG,TAG_P=truseq_library([seg("cell barcode","C"*16,"cbc",placeholder=True,feature=CELL),seg("tag UMI","U"*10,"umi",placeholder=True,feature=UMI),seg("spatial barcode","S"*20,"cbc",placeholder=True,feature=SPATIAL)],"Slide-tags spatial-tag library")
FINAL_LIBRARIES=(
 ("RNA library",RNA,RNA_P,"10x 3-prime or Multiome gene-expression library from tagged nuclei.","The droplet cell barcode and UMI identify each transcript."),
 ("ATAC library (Multiome branch)",ATAC,ATAC_P,"Tagmented accessibility library from the same spatially tagged nuclei.","Paired genomic reads use Nextera sites and retain shared cell identity."),
 ("Spatial-tag library",TAG,TAG_P,"Feature-barcode-style amplification links released puck tags to 10x nucleus barcodes.","The cell barcode groups spatial-tag UMIs; the tag sequence maps to bead coordinates."),
)
READ_LENGTHS={"RNA library":{"Read 1":28,"Read 2":90},"ATAC library (Multiome branch)":{"Read 1":50,"Read 2":50},"Spatial-tag library":{"Read 1":28,"Read 2":50}}

def sections(): return [
 ("Place tissue on a spatially sequenced bead array",[Row(chunks=[("bead — photocleavable linker — spatial barcode — 10x-capture sequence","cbc",False)])],"Each puck barcode already has a measured x/y coordinate."),
 ("Photocleave spatial oligos into the overlying tissue",[Row(chunks=[("UV cleavage → local spatial tags enter nearby nuclei","w1",False)])],"Short diffusion transfers several local bead identities into each intact nucleus."),
 ("Dissociate tagged nuclei and partition with 10x beads",[Row(chunks=[("spatial tag + one GEM cell barcode", "cbc",False)])],"Feature-barcode capture associates puck tags with a nucleus barcode."),
 ("Prepare assay and spatial-tag libraries",modality_split("RNA library","ATAC library in Multiome experiments","spatial-tag library"),"The separate tag library supplies coordinates for the same cell barcodes found in molecular assay libraries."),
]
