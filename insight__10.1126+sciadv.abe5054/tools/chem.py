"""INSIGHT barcoded NASBA-to-NGS chemistry."""
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, feature, revcomp

TITLE="INSIGHT — barcoded NASBA followed by pooled sequencing"
NOTES="01_insight.html"
SOURCE='Defining source: <a href="https://doi.org/10.1126/sciadv.abe5054">Wu et al., <i>Science Advances</i> (2021)</a>.'
SUMMARY="Paired 5-nt sample barcodes enter a SARS-CoV-2 amplicon during isothermal NASBA; completed reactions are pooled and receive flow-cell adapters and run indexes in one PCR."
CUSTOM_R1="AGCCAGCTCTGGAGAATTCTAATACGACTCACTATAGGGAGAAGG"
TARGET_F="CCAGCAACTGTTTGTGGACCTA"; TARGET_R="ATGGTTTAACAGGCACAGGTGT"
FINAL_LIBRARY=Construct([
 seg("P5",il.P5,"p5"),seg("i5","J"*8,"cbc",placeholder=True,feature=feature("run_i5","sample_index","random")),
 seg("custom T7-containing read arm",CUSTOM_R1,"r1"),
 seg("right sample barcode","B"*5,"cbc",placeholder=True,feature=feature("nasba_right","sample_index","combinatorial",group="nasba_sample",part="right")),
 seg("viral target amplicon","X"*46,placeholder=True),
 seg("left sample barcode","C"*5,"cbc",placeholder=True,feature=feature("nasba_left","sample_index","combinatorial",group="nasba_sample",part="left")),
 seg("TruSeq Read 2 arm reverse complement",revcomp(il.TRUSEQ_READ2),"r2"),
 seg("i7 reverse complement","I"*8,"cbc",placeholder=True,feature=feature("run_i7","sample_index","random")),
 seg("P7 reverse complement",il.P7_RC,"p7")],name="INSIGHT sequencing amplicon")
R1=sp.custom("Read 1","INSIGHT T7-containing sequencing primer",CUSTOM_R1,"Wu et al. 2021, Table 2")
SEQ_PRIMERS=(R1,sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
FINAL_CAPTION="The NASBA barcode pair identifies the original sample; the outer PCR supplies run indexes and flow-cell ends."
SEQUENCING_INTRO="INSIGHT uses its T7-containing custom primer from the P5 side and a standard TruSeq primer from the P7 side."
def oligos(): return [
 "<b>Barcoded NASBA forward primer</b> — 5′-TGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNCCAGCAACTGTTTGTGGACCTA-3′",
 "<b>Barcoded NASBA reverse primer</b> — 5′-AATTCTAATACGACTCACTATAGGGAGAAGGNNNNNACACCTGTGCCTGTTAAACCAT-3′",
 "<b>Custom sequencing primer</b> — 5′-AGCCAGCTCTGGAGAATTCTAATACGACTCACTATAGGGAGAAGG-3′"]
def sections(): return [
 ("Reverse transcribe the viral RNA",[Row(chunks=[("viral RNA ← barcoded forward primer","cbc",False)])],"The forward primer contributes one sample barcode and an Illumina handle."),
 ("RNase H and second-strand synthesis create a T7 template",[Row(chunks=[("Illumina arm—left barcode—viral target—right barcode—T7 promoter","cbc",False)])],"The reverse primer supplies the second barcode and T7 promoter."),
 ("T7 transcription and repeated NASBA cycling",[Row(chunks=[("dsDNA template → many ssRNA amplicons → dsDNA templates","r1",False)])],"Reverse transcriptase, RNase H and T7 RNA polymerase amplify isothermally."),
 ("Pool completed reactions and perform universal PCR",[Row(chunks=[("P5—run index—T7 read arm—barcode pair—target—TruSeq arm—run index—P7","cbc",False)])],"Residual NASBA primers are removed before the minimal-cycle pooled PCR.")]
