"""easySHARE-seq dual-modality, three-part combinatorial indexing."""
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Row, Scene, Workflow, feature
TITLE="easySHARE-seq — joint RNA and chromatin accessibility"
NOTES="01_easyshare-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1101/2024.02.26.581705">Ma et al., bioRxiv (2024)</a>; later published in <i>eLife</i>.'
SUMMARY="Fixed nuclei receive RNA RT handles and ATAC Tn5 handles, then two 7-nt ligation barcodes plus a third PCR index identify the same nucleus in both libraries."
R1_001="GGTTAGCGTCTCGT"
R1_001_SCAFFOLD="CCGAGCCCACGAGACGCTAACCTCT"
R2_001="CAAGCAGAAGACGGCATACGAGATTCTGGCAAGA"
R2_001_SCAFFOLD="TGCCAGAATCTCGTATGCCGTCTTCTGCTTG"

def barcode_ligation_geometry():
    sc=Scene()
    sc.strand("round-1 scaffold",[seg("scaffold duplex arm",R1_001_SCAFFOLD[:-3]),seg("TCT cohesive end",R1_001_SCAFFOLD[-3:])],label="Round1_001 scaffold")
    sc.strand("round-1 barcode",[seg("round-1 barcode oligo",R1_001,"cbc")],label="5′-phosphorylated Round1_001")
    sc.strand("round-2 barcode",[seg("round-2 barcode oligo",R2_001,"cbc")],label="Round2_001")
    sc.strand("round-2 scaffold",[seg("round-2 scaffold",R2_001_SCAFFOLD)],label="Round2_001 scaffold")
    sc.overhang("round-1 scaffold","TCT cohesive end",polarity="3-prime",owner="adapter",text="published 3′-TCT ligation overhang")
    sc.labels("round-1 scaffold")
    return AdapterLigationGeometry.known("cohesive",sc,statement="The easySHARE-seq supplementary oligo workbook prints each barcode with its scaffold and states that the scaffold leaves single-stranded ligation overhangs; Round1_001 exposes the terminal TCT shown here.")
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
 ("Ligate round-1 and round-2 barcodes",barcode_ligation_geometry().rows()+[Row(chunks=[("BC1 ** TCT ** BC2","cbc",False)])],"The supplementary workbook prints the oligo/scaffold pairs. TCT is the three-base linker retained between the two seven-base barcode reads."),
 ("Reverse crosslinks and split modalities",[Row(chunks=[("barcoded pool → cDNA branch | accessible-DNA branch","r1",False)])],"RNA and ATAC material are enriched separately."),
 ("Add round-3 barcode and sequencing ends by PCR",[Row(chunks=[("BC3 × BC1 × BC2 = cellular identity","r2",False)])],"A 96-way PCR index supplies the third cell-barcode component and permits standard sequencing."),]

def workflow():
    g=barcode_ligation_geometry()
    w=Workflow(MolecularState("RNA and ATAC handles",tuple(sections()[0][1])))
    w.ligate_adapter("Seal the two barcode rounds",g,[Row(chunks=[("BC1 ** TCT ** BC2","cbc",False)])],note="T4 DNA ligase seals the published cohesive barcode/scaffold system.")
    return w
