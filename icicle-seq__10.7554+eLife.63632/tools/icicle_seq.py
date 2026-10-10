"""ICICLE-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from antibody_tags import feature_library
from chromatin_epigenetics import droplet_targeted_tagment_library
from chemdraw import Row
ATAC,ATAC_P=droplet_targeted_tagment_library("accessible DNA",assay_name="poly(dT)/R2N capture branch",assay_barcode_len=8,assay_is_feature=False)
ADT,ADT_P=feature_library("ICICLE-seq antibody-tag library",barcode_id="icicle_antibody",cell_id="cell_barcode")
TITLE="ICICLE-seq — chromatin accessibility and epitopes"
NOTES="01_icicle-seq.html"
SOURCE='<a href="https://doi.org/10.7554/eLife.63632">Swanson et al., <i>eLife</i> (2021)</a>.'
SUMMARY="A custom Tn5 transposome with poly(dT)- or R2N-bearing ends makes ATAC fragments capturable by 10x 3′ gel beads alongside polyadenylated antibody tags."
CAVEAT="The displayed capture branch is structural: the paper's supplement defines several custom oligos, not a biological feature barcode."
FINAL_LIBRARIES=(("ATAC library",ATAC,ATAC_P,"Single-ended ATAC insert carrying a 10x cell barcode and UMI.","Read 1 is the 28-cycle cell-barcode/UMI read; Read 2 is the ATAC insert."),("Antibody-tag library",ADT,ADT_P,"ADT barcode from the same GEM is selectively amplified.","Read 1 carries cell barcode and UMI; Read 2 identifies the antibody."))
READ_LENGTHS={"ATAC library":{"Read 1":28,"Index 1 (i7)":8,"Read 2":100},"Antibody-tag library":{"Read 1":28,"Index 1 (i7)":8,"Read 2":100}}
def sections():return [("Stain cells with TotalSeq-A antibodies",[Row(chunks=[("antibody — ADT barcode — poly(A)","cbc",False)])],"Antibody-derived tags remain associated with each cell."),("Tagment with custom Tn5",[Row(chunks=[("poly(dT) or R2N capture end ** accessible DNA","me",False)])],"The custom transposome makes genomic fragments compatible with the 10x 3′ bead capture reaction."),("Barcode in GEMs",[Row(chunks=[("10x cell barcode + UMI → ATAC fragments and ADTs","cbc",False)])],"Both analytes receive one cell identity."),("Split and amplify",[Row(chunks=[("barcoded product ├─ ATAC PCR\n                 └─ ADT PCR",None,False)])],"The two libraries are sequenced separately.")]
