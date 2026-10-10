"""easySHARE-seq dual-modality, three-part combinatorial indexing."""
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, feature
TITLE="easySHARE-seq — joint RNA and chromatin accessibility"
NOTES="01_easyshare-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1101/2024.02.26.581705">Ma et al., bioRxiv (2024)</a>; later published in <i>eLife</i>.'
SUMMARY="Fixed nuclei receive RNA RT handles and ATAC Tn5 handles, then two 7-nt ligation barcodes plus a third PCR index identify the same nucleus in both libraries."
def _library(name,payload):
    lib=Construct([seg("P5",il.P5,"p5"),seg("cell barcode 3","C"*8,"cbc",placeholder=True,feature=feature("cell_bc3","cell_barcode","whitelist",group="cell_barcode",part="round 3",whitelist="easySHARE-seq 96-index set")),seg("Read 1 arm",il.TRUSEQ_READ1,"r1"),seg(payload,"X"*42,placeholder=True),seg("dA junction","A"),seg("Index 1 / Read 2 arm",il.INDEX1_PRIMER,"r2"),seg("cell barcode 1","B"*7,"cbc",placeholder=True,feature=feature("cell_bc1","cell_barcode","whitelist",group="cell_barcode",part="round 1",whitelist="easySHARE-seq 192-barcode set")),seg("3-nt ligation overhang","L"*3,placeholder=True),seg("cell barcode 2","D"*7,"cbc",placeholder=True,feature=feature("cell_bc2","cell_barcode","whitelist",group="cell_barcode",part="round 2",whitelist="easySHARE-seq 192-barcode set")),seg("P7 reverse complement",il.P7_RC,"p7")],name=name)
    primers=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["I2"],sp.TRUSEQ["R2"])
    problems=sp.verify(lib,primers)
    if problems: raise ValueError(name+": "+"; ".join(problems))
    return lib,primers
RNA,RP=_library("easySHARE-seq RNA library","cDNA insert")
ATAC,AP=_library("easySHARE-seq ATAC library","accessible chromatin insert")
FINAL_LIBRARIES=(("RNA library",RNA,RP,"Transcript-derived inserts carry the same three-part cellular identity.","Standard Illumina primers read insert and barcode/index reads."),("ATAC library",ATAC,AP,"Accessible-DNA inserts carry the same three-part cellular identity.","Standard Illumina primers read insert and barcode/index reads."))
def oligos(): return ["<b>Tn5-B2S</b> — 5′-/5Phos/GGGCTCGGAGATGTGTATAAGAGACAG-3′"]
def sections(): return [
 ("Mark RNA and accessible DNA in fixed nuclei",[Row(chunks=[("RNA → RT handle     chromatin ** Tn5 handles","me",False)])],"Both analytes remain associated with each nucleus."),
 ("Ligate round-1 barcodes",[Row(chunks=[("RNA/ATAC handle ** 7-nt BC1 ** 3-nt overhang","cbc",False)])],"The first 192-way split-pool ligation starts the cell identity."),
 ("Ligate round-2 barcodes",[Row(chunks=[("BC1 ** 3-nt junction ** 7-nt BC2","cbc",False)])],"The second 192-way ligation completes the 17-cycle Index 1 structure."),
 ("Reverse crosslinks and split modalities",[Row(chunks=[("barcoded pool → cDNA branch | accessible-DNA branch","r1",False)])],"RNA and ATAC material are enriched separately."),
 ("Add round-3 barcode and sequencing ends by PCR",[Row(chunks=[("BC3 × BC1 × BC2 = cellular identity","r2",False)])],"A 96-way PCR index supplies the third cell-barcode component and permits standard sequencing."),]
