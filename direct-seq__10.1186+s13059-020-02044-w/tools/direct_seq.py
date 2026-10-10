"""Direct-seq programmed-guide capture model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row,Segment,oligo
from tenx_rna import SEQ_PRIMERS,poly_a_feature_library,transcript_library

CAPTURE="AAAAAAAAGAAAAAAAGAAAAAAAGAAAAA"
RNA_LIB=transcript_library("Direct-seq gene-expression library")
GUIDE_LIB=poly_a_feature_library("gRNA spacer",feature_len=20,feature_role="guide_barcode",feature_id="direct_guide",name="Direct-seq enriched guide library")
FINAL_LIBRARIES=(
 ("Gene-expression library",RNA_LIB,SEQ_PRIMERS,"Standard 10x 3′ v3 transcript library.","Read 1 records cell barcode and UMI; Read 2 enters endogenous cDNA."),
 ("Enriched guide library",GUIDE_LIB,SEQ_PRIMERS,"Nested PCR enriches cDNA copied directly from the programmed guide transcript.","Read 1 records the same cell barcode and UMI; Read 2 identifies the guide spacer."),)
TITLE="Direct-seq — poly(dT) capture of programmed guide RNA"
NOTES="01_direct-seq.html"
SOURCE='<a href="https://doi.org/10.1186/s13059-020-02044-w">Song et al., <i>Genome Biology</i> (2020)</a>.'
SUMMARY="A 30-nt A/G capture tract is inserted into the sgRNA scaffold. Standard poly(dT) reverse transcription therefore copies guide transcripts alongside mRNA, after which nested PCR produces a guide-enriched library sharing the cell barcode."
READ_LENGTHS={"Gene-expression library":{"Read 1":28},"Enriched guide library":{"Read 1":28}}
def oligos(): yield oligo("8A8G capture tract",[Segment("A/G mixed capture",CAPTURE,"tso")])
def sections(): return [
 ("Program the sgRNA scaffold",[Row(chunks=[("sgRNA scaffold — 8A8G capture tract — guide spacer",None,False)])],"The mixed tract avoids long consecutive A runs while remaining accessible to poly(dT)."),
 ("Capture mRNA and guide RNA in the same GEM",[Row(chunks=[("gel-bead cell barcode — UMI — poly(dT) || mRNA or programmed sgRNA", "cbc",False)])],"No guide-specific reverse-transcription primer is added."),
 ("Pre-amplify and split",[Row(chunks=[("shared barcoded cDNA ├─ standard mRNA library\n                     └─ nested guide PCR",None,False)])],"Both outputs retain the same cell barcode."),]
