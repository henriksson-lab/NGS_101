"""TEA-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from antibody_tags import feature_library
from batch_ngs import nextera_library,seg
from chemdraw import Row,feature
from tenx_rna import SEQ_PRIMERS,transcript_library
CELL=feature("tea_cell","cell_barcode","whitelist",whitelist="10x Multiome bead barcode")
RNA=transcript_library("TEA-seq transcript library")
ATAC,ATAC_P=nextera_library([seg("shared cell barcode","C"*16,"cbc",placeholder=True,feature=CELL),seg("accessible DNA","X"*42,placeholder=True)],"TEA-seq ATAC library")
ADT,ADT_P=feature_library("TEA-seq antibody-tag library",barcode_id="tea_antibody",cell_id="tea_cell")
TITLE="TEA-seq — transcripts, epitopes and accessibility"
NOTES="01_tea-seq.html"
SOURCE='<a href="https://doi.org/10.7554/eLife.63632">Swanson et al., <i>eLife</i> (2021)</a>.'
SUMMARY="Oligo-tagged antibodies are combined with a 10x Single Cell Multiome workflow so RNA, surface epitopes and accessible chromatin share one cell identity but become three sequencing libraries."
CAVEAT="Commercial bead and TotalSeq-A sequences are represented structurally. The page does not claim proprietary bases."
FINAL_LIBRARIES=(("RNA library",RNA,SEQ_PRIMERS,"Poly(A) RNA captured by the Multiome bead.","Read 1 carries cell barcode and UMI; Read 2 reports cDNA."),("ATAC library",ATAC,ATAC_P,"Tagmented accessible chromatin carrying the same cell identity.","Paired genomic reads use Nextera sites."),("Antibody-tag library",ADT,ADT_P,"Polyadenylated antibody tags selectively amplified from the shared product.","Barcode and UMI identify the epitope observation."))
READ_LENGTHS={"RNA library":{"Read 1":28,"Read 2":90},"ATAC library":{"Read 1":51,"Index 1 (i7)":8,"Index 2 (i5)":16,"Read 2":51},"Antibody-tag library":{"Read 1":28,"Read 2":50}}
def sections():return [("Stain permeabilized cells",[Row(chunks=[("epitope — antibody — feature barcode — poly(A)","cbc",False)])],"TotalSeq-A oligos encode surface epitopes."),("Tagment and partition",[Row(chunks=[("accessible DNA + poly(A) RNA + antibody tags → one barcoded GEM","me",False)])],"The Multiome workflow associates all three analytes with the same cell."),("Split by analyte",[Row(chunks=[("shared product ├─ RNA\n               ├─ ATAC\n               └─ antibody tags",None,False)])],"Modality-specific PCRs create three libraries.")]
