"""ProBac-seq probe-capture model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row
from tenx_rna import SEQ_PRIMERS,poly_a_feature_library
FINAL_LIBRARY=poly_a_feature_library("gene-specific probe identity",feature_len=50,feature_role="feature_barcode",feature_id="probac_probe",name="ProBac-seq probe library")
TITLE="ProBac-seq — probe-based bacterial single-cell transcript counting"
NOTES="01_probac-seq.html"
SOURCE='<a href="https://doi.org/10.1038/s41596-024-01002-1">Samanta et al., <i>Nature Protocols</i> (2024)</a>.'
SUMMARY="Amplified gene-specific DNA probes receive a 12-nt UMI and poly(A)30, hybridize to fixed bacterial transcripts, and are PCR-copied with a 10x gel-bead primer so each bound probe gains a cell barcode."
FINAL_CAPTION="Representative captured probe. Probe identity reports the target gene; the probe-specific UMI was installed before hybridization and the 10x cell barcode in the GEM."
SEQUENCING_INTRO="The protocol uses an 8-cycle i7 read and a 119-cycle Read 1 on NextSeq 2000 with 30% PhiX; the drawn canonical endpoint exposes the cell barcode, UMI and gene-specific probe region."
READ_LENGTHS={"Read 1":119,"Index 1 (i7)":8}
def sections(): return [
 ("Design and RCA-amplify gene-specific probes",[Row(chunks=[("PCR handle — 50-nt transcript-binding region — HindIII extender",None,False)])],"Three RCA rounds produce enough probe in the orientation complementary to bacterial mRNA."),
 ("Add probe UMI and poly(A)",[Row(chunks=[("probe identity — random N12 UMI — poly(A)30", "umi",False)])],"Isothermal extension gives every probe molecule a counting UMI and a 10x-compatible tail."),
 ("Hybridize in fixed, permeabilized bacteria",[Row(chunks=[("probe binding region || target mRNA",None,False)])],"Unbound and mismatched probes are stringently washed away."),
 ("Encapsulate and barcode",[Row(chunks=[("10x cell barcode — probe UMI — target-probe identity", "cbc",False)])],"A six-cycle in-GEM PCR replaces reverse transcription and copies bound DNA probes."),
 ("Scale-up PCR and indexed library PCR",[Row(chunks=[("PCR-2 enrichment → PCR-3 adds Illumina i5/i7 arms",None,False)])],"The final product is an Illumina-compatible probe-counting library."),]
