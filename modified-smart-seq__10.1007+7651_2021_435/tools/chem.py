"""Xi Chen laboratory Modified SMART-Seq chemistry."""
import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, feature

TITLE="Modified SMART-Seq — 3′ counting with well barcodes and UMIs"
NOTES="01_modified-smart-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1007/7651_2021_435">Zheng et al., <i>Methods in Molecular Biology</i> (2022)</a>.'
SUMMARY="A barcoded oligo(dT) installs an 8-nt well barcode and 10-nt UMI during reverse transcription; template switching and cDNA PCR precede tagmentation and indexed library PCR."
RT="CTACACGACGCTCTTCCGATCT"
TSO="AAGCAGTGGTATCAACGCAGAGTACATGGG"
FINAL_LIBRARY=Construct([
    seg("P5",il.P5,"p5"),seg("Read 1 arm",RT,"r1"),
    seg("well barcode","B"*8,"cbc",placeholder=True,feature=feature(
        "well_barcode","cell_barcode","whitelist",
        whitelist="Zheng et al. 2022, Table 1 well-barcode set")),
    seg("UMI","U"*10,"umi",placeholder=True,feature=feature("umi","umi","random")),
    seg("transcript cDNA","X"*36,placeholder=True),seg("TSO-derived end",TSO),
    seg("mosaic end reverse complement",nx.ME_RC,"me"),seg("S7 reverse complement",nx.S7_RC,"s7"),
    seg("plate index reverse complement","I"*8,"cbc",placeholder=True,feature=feature(
        "plate_index","sample_index","whitelist",
        whitelist="Zheng et al. 2022, Table 1 plate-barcode set")),
    seg("P7 reverse complement",il.P7_RC,"p7")],name="Modified SMART-Seq library")
SEQ_PRIMERS=(sp.TRUSEQ["R1"],sp.NEXTERA["I1"],sp.NEXTERA["R2"])
FINAL_CAPTION="Read 1 begins with the well barcode and UMI; Read 2 reads transcript sequence from the tagmented end."
SEQUENCING_INTRO="The protocol specifies at least 18 Read 1 cycles, an 8-cycle plate-index read and at least 50 Read 2 cycles."
READ_LENGTHS={"Read 1":18,"Index 1 (i7)":8,"Read 2":75}
SHOW_READ_CONTENT=True
def oligos(): return [
 "<b>Barcoded RT primer</b> — 5′-CTACACGACGCTCTTCCGATCT[8-bp well barcode]NNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTVN-3′",
 "<b>Template-switching oligo</b> — 5′-AAGCAGTGGTATCAACGCAGAGTACATrGrGrG-3′"]
def sections(): return [
 ("Reverse transcription installs well barcode and UMI",[Row(chunks=[("poly(A) RNA ← oligo(dT)—UMI—well barcode—Read 1 arm","umi",False)])],"Each well receives its own RT primer before cells are sorted."),
 ("Template switch at the cDNA 5′ end",[Row(chunks=[("first-strand cDNA ** TSO-derived PCR handle","r1",False)])],"The terminal riboguanosines support template switching."),
 ("Pool wells after cDNA PCR and Exonuclease I cleanup",[Row(chunks=[("well-tagged full-length cDNA pool",None,False)])],"Primer removal precedes pooling and suppresses cross-well mispriming."),
 ("Tagment and add the plate index by library PCR",[Row(chunks=[("P5—well barcode—UMI—cDNA ** mosaic end—plate index—P7","cbc",False)])],"Tn5 fragments amplified cDNA; the Nextera-style index primer identifies the plate.")]
