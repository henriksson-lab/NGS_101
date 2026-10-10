"""DOGMA-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from antibody_tags import feature_library
from batch_ngs import nextera_library, seg, truseq_library
from chemdraw import Row, feature
from multimodal_spatial import modality_split

TITLE = "DOGMA-seq — RNA, ATAC and protein from one cell"
NOTES = "01_dogma-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41587-021-00927-2">Mimitou et al., <i>Nature Biotechnology</i> (2021)</a>.'
SUMMARY = "A 10x Multiome-compatible workflow jointly captures poly(A) RNA, accessible chromatin and antibody-derived tags from each permeabilized cell."
CAVEAT = "The three outputs share cell identity but are distinct sequencing libraries. Proprietary bead and antibody-tag regions remain structural placeholders."
CELL = feature("dogma_cell", "cell_barcode", "whitelist", whitelist="10x Chromium Next GEM bead barcode")
UMI = feature("dogma_rna_umi", "umi", "random")
RNA, RNA_P = truseq_library([seg("cell barcode", "C"*16,"cbc",placeholder=True,feature=CELL),seg("RNA UMI","U"*12,"umi",placeholder=True,feature=UMI),seg("3-prime cDNA","X"*40,placeholder=True)], "DOGMA-seq RNA library")
ATAC, ATAC_P = nextera_library([seg("cell barcode association", "C"*16,"cbc",placeholder=True,feature=CELL),seg("accessible genomic DNA","X"*40,placeholder=True)], "DOGMA-seq ATAC library")
ADT, ADT_P = feature_library("DOGMA-seq antibody-tag library", barcode_id="dogma_antibody",cell_id="dogma_cell")
FINAL_LIBRARIES = (
 ("RNA library",RNA,RNA_P,"Poly(dT) capture, GEM barcoding and cDNA library construction recover gene expression.","Read 1 carries cell barcode/UMI; the insert read reports cDNA."),
 ("ATAC library",ATAC,ATAC_P,"Tagmented chromatin is barcoded and PCR-completed as the accessibility library.","Paired genomic reads use Nextera sites."),
 ("Antibody-tag library",ADT,ADT_P,"Antibody-derived oligos are selectively amplified from the same emulsion product.","Feature barcode and UMI are read separately from genomic and cDNA inserts."),
)
READ_LENGTHS={"RNA library":{"Read 1":28,"Read 2":90},"ATAC library":{"Read 1":50,"Read 2":50},"Antibody-tag library":{"Read 1":28,"Read 2":50}}

def sections(): return [
 ("Stain cells with DNA-barcoded antibodies",[Row(chunks=[("antibody — feature barcode — UMI — capture handle","cbc",False)])],"Antibody tags encode the detected epitope."),
 ("Tagment accessible chromatin",[Row(chunks=[("Tn5 inserts adapters at accessible DNA","me",False)])],"Tagmentation occurs while RNA and antibody tags remain associated with the cell."),
 ("Partition cells with barcoded gel beads",[Row(chunks=[("one GEM barcode → RNA + ATAC + antibody tag","cbc",False)])],"The shared bead identity links the three measurements computationally."),
 ("Recover three modality-specific libraries",modality_split("RNA library","ATAC library","antibody-tag library"),"Separate enrichments preserve a common cell barcode while using modality-specific primers."),
]
