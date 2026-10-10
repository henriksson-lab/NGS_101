"""SEC-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row
from tenx_rna import SEQ_PRIMERS,poly_a_feature_library,transcript_library

TITLE="SEC-seq — secretion-linked single-cell RNA sequencing"
NOTES="01_sec-seq.html"
SOURCE='<a href="https://doi.org/10.1038/s41467-023-39367-8">Cheng et al., <i>Nature Communications</i> (2023)</a>.'
SUMMARY="A cell is held in an antibody-coated nanovial that captures its secreted protein. An oligo-barcoded detection antibody and the cell transcriptome are then read together with 10x Chromium."
CAVEAT="The publication uses commercial TotalSeq oligos; their feature identity and capture architecture are public, but the product sequence is proprietary."
RNA=transcript_library("SEC-seq gene-expression library")
SEC=poly_a_feature_library("secretion-antibody barcode",feature_len=15,feature_id="sec_antibody",name="SEC-seq secretion-tag library")
FINAL_LIBRARIES=(("Gene-expression library",RNA,SEQ_PRIMERS,"10x 3′ v3.1 transcript library from the nanovial-bound cell.","Read 1 records cell barcode and UMI; Read 2 reports cDNA."),("Secretion-tag library",SEC,SEQ_PRIMERS,"Feature-barcode library from the oligo-labeled antibody bound to captured secretion.","Read 1 records the same cell barcode and UMI; Read 2 identifies the secretion antibody."))
READ_LENGTHS={"Gene-expression library":{"Read 1":28,"Index 1 (i7)":10,"Read 2":90},"Secretion-tag library":{"Read 1":28,"Index 1 (i7)":10,"Read 2":90}}
def sections(): return [
 ("Capture one cell in a nanovial",[Row(chunks=[("antibody-coated cavity ← cell",None,False)])],"The cavity retains the cell and concentrates secreted protein on its surface."),
 ("Label captured secretion",[Row(chunks=[("nanovial — secreted IgG — detection antibody — DNA feature barcode","cbc",False)])],"An oligo-barcoded secondary antibody converts secretion into a sequenceable feature."),
 ("Partition cell-loaded nanovials",[Row(chunks=[("one GEM barcode → cellular mRNA + nanovial antibody tags","cbc",False)])],"The nanovial fits through the Chromium microfluidic channel and does not replace the standard bead chemistry."),
 ("Amplify two libraries",[Row(chunks=[("barcoded cDNA ├─ gene expression\n              └─ secretion feature tags",None,False)])],"The paper pools 80% RNA and 20% antibody-tag library for sequencing."),]
