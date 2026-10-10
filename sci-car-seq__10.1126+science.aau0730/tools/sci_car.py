"""Cao et al. sci-CAR RNA/ATAC co-assay."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import joined_scene, nextera_library, seg, truseq_library
from chemdraw import Row, feature
from multimodal_spatial import modality_split

TITLE="sci-CAR-seq — joint combinatorial RNA and ATAC profiling"
NOTES="01_sci-car-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1126/science.aau0730">Cao et al., <i>Science</i> (2018)</a>.'
SUMMARY="One intact nucleus receives a barcoded poly(dT) RT primer and a separately barcoded Tn5 adapter in the first plate; FACS redistribution supplies a second well index while RNA and ATAC aliquots become separate libraries."
CAVEAT="RNA and ATAC first-round indexes are chemically distinct but map to the same first well; their second-round PCR indexes likewise map to the same FACS well. That correspondence—not a shared physical barcode molecule—links the two libraries."
RNA1=feature("rna_rt_well","cell_barcode","combinatorial",group="cell_id",part="RT well",whitelist="published sci-CAR RT primer plate")
ATAC1=feature("atac_tn5_well","cell_barcode","combinatorial",group="cell_id",part="Tn5 well",whitelist="published sci-CAR Tn5 plate")
RNA2=feature("rna_pcr_well","cell_barcode","combinatorial",group="cell_id",part="FACS/PCR well",whitelist="published sci-CAR PCR plate")
ATAC2=feature("atac_pcr_well","cell_barcode","combinatorial",group="cell_id",part="FACS/PCR well",whitelist="published sci-CAR PCR plate")
UMI=feature("rna_umi","umi","random")
RNA,RNA_P=truseq_library([seg("RNA PCR well index","B"*10,"cbc",placeholder=True,feature=RNA2),seg("RT well index","A"*10,"cbc",placeholder=True,feature=RNA1),seg("UMI","U"*8,"umi",placeholder=True,feature=UMI),seg("3-prime cDNA tag","X"*38,placeholder=True)],"sci-CAR RNA library")
ATAC,ATAC_P=nextera_library([seg("Tn5 well index","A"*10,"cbc",placeholder=True,feature=ATAC1),seg("accessible DNA","X"*38,placeholder=True),seg("PCR well index","B"*10,"cbc",placeholder=True,feature=ATAC2)],"sci-CAR ATAC library")
FINAL_LIBRARIES=(("RNA library",RNA,RNA_P,"Barcoded poly(dT) RT, second-strand synthesis, unindexed tagmentation and indexed PCR produce 3-prime cDNA tags.","Read/index cycles recover the RT barcode, UMI, PCR-well barcode and cDNA."),("ATAC library",ATAC,ATAC_P,"Barcoded Tn5 insertions are completed by a second well-specific PCR index.","Paired genomic reads are associated with the first Tn5 well and second PCR well."))
READ_LENGTHS={"RNA library":{"Read 1":25,"Read 2":50},"ATAC library":{"Read 1":50,"Read 2":50}}

def sections(): return [
 ("Index RNA by in-situ reverse transcription",[Row(chunks=[("poly(dT) -- UMI -- RT-well barcode ---> poly(A) RNA", "cbc",False)])],"A well-specific poly(dT) primer copies mRNA and installs the first RNA index plus UMI."),
 ("Index accessible DNA by in-situ tagmentation",joined_scene([seg("left barcoded Tn5 adapter","B"*18,"cbc",placeholder=True),seg("accessible genomic DNA","X"*34,placeholder=True),seg("right barcoded Tn5 adapter","B"*18,"cbc",placeholder=True)],(("left barcoded Tn5 adapter","accessible genomic DNA","Tn5 transfer"),("accessible genomic DNA","right barcoded Tn5 adapter","Tn5 transfer")),label="barcoded accessible-DNA fragment").rows(),"A separately barcoded transposome installs the first ATAC index in the same intact nuclei."),
 ("Pool nuclei and redistribute by FACS",[Row(chunks=[("RT/Tn5 well identity + FACS well identity = cell identity", "cbc",False)])],"Most nuclei traverse a unique pair of wells; matched well coordinates link modalities."),
 ("Make second-strand cDNA, lyse and split each well",modality_split("RNA-dedicated lysate","ATAC-dedicated lysate"),"Only after second-strand synthesis is each lysate divided into the two library branches."),
 ("Complete the two indexed libraries",modality_split("unindexed Tn5 on cDNA + RNA-index PCR","ATAC-index PCR"),"The second well-specific index is added by modality-specific PCR before separate pooling and sequencing."),
]
