"""sci-Plex nuclear-hashing model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg,truseq_library
from chemdraw import Row,feature

FINAL_LIBRARY,SEQ_PRIMERS=truseq_library([
 seg("sci-RNA cell barcode combination","C"*20,"cbc",placeholder=True,feature=feature("sciplex_cell","cell_barcode","combinatorial",group="cell",part="sci-RNA rounds")),
 seg("UMI","U"*8,"umi",placeholder=True,feature=feature("sciplex_umi","umi","random")),
 seg("hash barcode or cDNA insert","X"*36,placeholder=True),
],"sci-Plex library",dual_index=True)
TITLE="sci-Plex — nuclear hashing before combinatorial RNA indexing"
NOTES="01_sci-plex.html"
SOURCE='<a href="https://doi.org/10.1126/science.aax6234">Srivatsan et al., <i>Science</i> (2020)</a>.'
SUMMARY="Each treated sample is labeled with a well-specific polyadenylated ssDNA hash oligo. Pooled nuclei then undergo sci-RNA-seq, so hash oligos and endogenous transcripts receive the same combinatorial cell indexes."
FINAL_CAPTION="A sci-RNA-seq product representing either endogenous cDNA or a captured nuclear-hash oligo; both retain the same cell barcode combination and UMI."
SEQUENCING_INTRO="The combinatorial index combination identifies a nucleus; hash-derived reads recover treatment/sample identity and transcript-derived reads measure expression."
def sections(): return [
 ("Treat samples in separate wells",[Row(chunks=[("sample 1 | sample 2 | … | sample 4992",None,False)])],"Perturbation identity initially exists as physical well identity."),
 ("Hash permeabilized nuclei",[Row(chunks=[("well-specific ssDNA hash barcode — poly(A) → nucleus", "cbc",False)])],"Unmodified polyadenylated DNA associates with exposed nuclei and is fixed in place."),
 ("Pool and perform sci-RNA-seq",[Row(chunks=[("hash oligo + endogenous mRNA → same RT/ligation/PCR index combination",None,False)])],"The poly(A) tail makes the hash compatible with the transcript workflow."),]
