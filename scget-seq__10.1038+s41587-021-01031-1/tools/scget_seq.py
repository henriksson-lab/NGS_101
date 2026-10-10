"""scGET-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row
from chromatin_epigenetics import droplet_targeted_tagment_library

TITLE="scGET-seq — open and heterochromatic DNA in the same cell"
NOTES="01_scget-seq.html"
SOURCE='<a href="https://doi.org/10.1038/s41587-021-01031-1">Tedesco et al., <i>Nature Biotechnology</i> (2022)</a>.'
SUMMARY="Ordinary barcoded Tn5 marks accessible chromatin while HP1α–Tn5 carrying a distinct adaptor barcode marks H3K9me3-associated heterochromatin; Chromium ATAC supplies a shared cell identity."
CAVEAT="The molecule shows one HP1α–Tn5 modality barcode. Parallel ordinary-Tn5 products use a different published barcode and share the same droplet cell barcode."
FINAL_LIBRARY,SEQ_PRIMERS=droplet_targeted_tagment_library("accessible or H3K9me3-associated DNA",assay_name="Tn5-versus-TnH barcode")
FINAL_CAPTION="Representative scGET-seq product. The inline transposase barcode identifies open-chromatin Tn5 versus HP1α-targeted TnH; the 10x barcode identifies the cell."
SEQUENCING_INTRO="Paired genomic reads identify the insertion locus. A protocol barcode distinguishes the two transposases and the Chromium barcode joins them to one cell."
SEQUENCING_UNAVAILABLE="The complete dedicated 10x cell-barcode sequencing-primer oligo is vendor-supplied; genomic Read 1, Index 1 and Read 2 landing sites are shown."
def sections(): return [
 ("Tag open chromatin",[Row(chunks=[("ordinary barcoded Tn5 → accessible DNA",None,False)])],"A modality-specific adaptor records accessible insertions."),
 ("Tag heterochromatin",[Row(chunks=[("HP1α–Tn5 (TnH) → H3K9me3-associated DNA",None,False)])],"The HP1α fusion supplies targeting and carries a distinct adaptor barcode."),
 ("Encapsulate nuclei in Chromium ATAC droplets",[Row(chunks=[("both modality tags + one 10x cell barcode", "cbc",False)])],"The same gel-bead barcode links both chromatin modalities to one nucleus."),
 ("PCR-complete and sequence",[Row(chunks=[("P5 — cell barcode — modality barcode — ME — genomic insert — ME — P7",None,False)])],"PCR produces the paired-end library."),]
