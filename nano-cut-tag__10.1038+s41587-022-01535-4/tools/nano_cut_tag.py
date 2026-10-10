"""nano-CUT&Tag molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row,Segment,oligo
from chromatin_epigenetics import antibody_enzyme_rows,droplet_targeted_tagment_library
BARCODE_A="ACGCTATAGCCT"
FINAL_LIBRARY,SEQ_PRIMERS=droplet_targeted_tagment_library("nanobody-targeted chromatin",assay_barcode_len=12,assay_name="nano-Tn5 modality barcode",connector_seq="GCGATCGAGGACGGC")
TITLE="nano-CUT&Tag — multimodal nanobody–Tn5 chromatin profiling"
NOTES="01_nano-cut-tag.html"
SOURCE='<a href="https://doi.org/10.1038/s41587-022-01535-4">Bartosovic et al., <i>Nature Biotechnology</i> (2023)</a>.'
SUMMARY="Mouse- and rabbit-specific secondary nanobody–Tn5 fusions transfer distinct MeA barcodes at antibody-bound chromatin. Chromium linear amplification adds cell identity, then MeB-loaded Tn5 installs the opposite library end."
FINAL_CAPTION="Representative nano-CT molecule: one published nano-Tn5 modality barcode plus a 10x cell barcode and the second-round MeB end."
SEQUENCING_INTRO="Paired genomic reads locate the targeted fragment; inline MeA barcodes identify the nanobody–Tn5 modality and the Chromium barcode identifies the cell."
SEQUENCING_UNAVAILABLE="The vendor gel-bead primer is not printed in full; the contiguous S5-side genomic, i7 and S7-side genomic primer sites are shown."
def oligos(): yield oligo("Tn5_P5_MeA_BcdA_0N",[Segment("S5","TCGTCGGCAGCGTC","s5"),Segment("connector","TCC"),Segment("barcode A",BARCODE_A,"cbc"),Segment("connector","GCGATCGAGGACGGC"),Segment("ME","AGATGTGTATAAGAGACAG","me")])
def sections(): return [
 ("Bind primary antibodies and nanobody–Tn5",antibody_enzyme_rows("nanobody–Tn5","MeA/barcode transfer"),"Nanobody specificity links each barcoded Tn5 to a primary-antibody species."),
 ("Deterministic targeted tagmentation",[Row(chunks=[("MeA — modality barcode — ME ** target DNA",None,False)])],"Distinct loaded fusions mark up to two histone features; an optional ordinary Tn5 marks accessibility."),
 ("Chromium linear amplification",[Row(chunks=[("10x cell barcode → repeated copies of each tagged fragment", "cbc",False)])],"Droplet barcoding attaches a shared cell identity before endpoint PCR."),
 ("Second MeB tagmentation and PCR",[Row(chunks=[("P5 — cell barcode — modality barcode — insert — MeB — P7",None,False)])],"A second Tn5 reaction supplies the opposite adapter and PCR completes the library."),]
