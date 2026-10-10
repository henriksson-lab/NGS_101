"""NTT-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import seqprimers as sp
from chemdraw import Row,Segment,oligo
from chromatin_epigenetics import antibody_enzyme_rows,droplet_targeted_tagment_library
MEDSA1="TCGTCGGCAGCGTCGGATTGCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG"
CUSTOM_R1="GCGATCGAGGACGGCAGATGTGTATAAGAGACAG"
CUSTOM_I5="CTGTCTCTTATACACATCTGCCGTCCTCGATCGC"
FINAL_LIBRARY,_BASE=droplet_targeted_tagment_library("nanobody-targeted chromatin",assay_name="nanobody–Tn5 barcode",connector_seq="GCGATCGAGGACGGC")
SEQ_PRIMERS=(sp.custom("Read 1","NTT custom Read 1",CUSTOM_R1,"Extended Data Table 2"),_BASE[1],_BASE[2])
TITLE="NTT-seq — nanobody-tethered transposition with single-cell indexing"
NOTES="01_ntt-seq.html"
SOURCE='<a href="https://doi.org/10.1038/s41587-022-01588-5">Stuart et al., <i>Nature Biotechnology</i> (2023)</a>.'
SUMMARY="Secondary nanobody–Tn5 fusions recognize primary-antibody species and transfer distinct 8-nt adaptor barcodes at several chromatin targets before 10x Chromium ATAC cell indexing."
FINAL_CAPTION="Representative NTT-seq chromatin molecule. The transposome barcode records the antibody/nanobody modality and the droplet barcode records the cell."
SEQUENCING_INTRO="Custom Read 1 enters through the connector and ME; i5 is 38 cycles to recover cell and assay structure, i7 is 8 cycles, and paired genomic reads are 60 cycles (50 cycles in the stated PBMC run without surface proteins)."
SEQUENCING_UNAVAILABLE="The public page places the printed custom Read 1 and genomic Read 2 sites. The complete 10x cell-barcode primer remains vendor-supplied."
READ_LENGTHS={"Read 1":60,"Index 1 (i7)":8,"Read 2":60}
def oligos():
 yield oligo("MEDSA_1",[Segment("S5","TCGTCGGCAGCGTC","s5"),Segment("8-nt modality barcode","GGATTGCT","cbc"),Segment("connector","GCGATCGAGGACGGC"),Segment("ME","AGATGTGTATAAGAGACAG","me")])
 yield oligo("Custom R1",[Segment("binding region",CUSTOM_R1,"r1")])
 yield oligo("Custom i5",[Segment("binding region",CUSTOM_I5,"r2")])
def sections(): return [
 ("Bind primary-antibody mixture",[Row(chunks=[("chromatin targets ← species/subclass-specific primary antibodies",None,False)])],"Each target is assigned to a compatible secondary nanobody."),
 ("Recruit barcoded nanobody–Tn5",antibody_enzyme_rows("nanobody–Tn5","barcode-specific tagmentation"),"Each fusion carries one of the printed MEDSA barcode adaptors."),
 ("Tagment and enter Chromium ATAC",[Row(chunks=[("target DNA ** modality barcode → 10x droplet cell barcode", "cbc",False)])],"Bulk targeted insertions retain their modality before cell-specific linear amplification."),
 ("PCR and custom sequencing",[Row(chunks=[("i5 38 bp | i7 8 bp | Read 1 60 bp | Read 2 60 bp",None,False)])],"Custom primers recover the nonstandard inserted barcode architecture."),]
