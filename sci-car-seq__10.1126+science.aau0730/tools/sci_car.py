"""Cao et al. sci-CAR RNA/ATAC co-assay."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import joined_scene, nextera_library, seg, truseq_library
from chemdraw import Row, feature
from multimodal_spatial import finish_pcr_branch, modality_split
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

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

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 nucleus=g.input("nucleus","intact nucleus containing accessible chromatin and poly(A) RNA",MaterialKind.BIOLOGICAL,state=molecule("intact nucleus",(role("accessible chromatin"),role("poly(A) RNA"))),primary=True)
 rtprimer=g.input("rt_primer","well-indexed poly(dT)/UMI RT primer",MaterialKind.OLIGO)
 rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
 copied=g.product("copied","nucleus retaining first-well-indexed cDNA and chromatin",MaterialKind.MOLECULAR,state=molecule("RT-indexed nucleus",(role("RT-well barcode"),role("UMI"),role("first-strand cDNA"),role("accessible chromatin"))))
 g.react("reverse_transcribe","Index RNA by in-situ reverse transcription",StepKind.REVERSE_TRANSCRIPTION,((nucleus,InputRole.TEMPLATE),(rtprimer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(copied,))
 transposome=g.input("atac_transposome","separately well-indexed ATAC transposome",MaterialKind.ADAPTER)
 tn5=g.input("atac_tn5","Tn5 transposase",MaterialKind.CATALYST)
 tagged=g.product("tagged","nucleus carrying indexed cDNA and tagmented chromatin",MaterialKind.MOLECULAR,state=copied.state)
 twaste=g.product("tn5_waste","non-library tagmentation products",MaterialKind.WASTE)
 g.react("tagment_atac","Index accessible DNA by in-situ tagmentation",StepKind.TAGMENTATION,((copied,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,twaste)); g.discard(twaste,"does not retain an amplifiable accessible-DNA fragment")
 facs=g.input("facs_wells","FACS redistribution and second-well indexing condition",MaterialKind.CONDITION)
 redistributed=g.product("redistributed","uniquely redistributed indexed nuclei",MaterialKind.MOLECULAR,state=tagged.state)
 excluded=g.product("excluded","empty or colliding FACS-well material",MaterialKind.WASTE)
 g.react("redistribute","Pool nuclei and redistribute them by FACS",StepKind.PARTITION,((tagged,InputRole.TEMPLATE),(facs,InputRole.CONDITION)),(redistributed,excluded),changes_molecule=False); g.discard(excluded,"does not preserve a unique pair of well coordinates")
 second_primer=g.input("second_primer","RNA second-strand synthesis primer",MaterialKind.OLIGO)
 second_pol=g.input("second_polymerase","RNA second-strand polymerase",MaterialKind.CATALYST)
 ds=g.product("double_stranded","nuclei containing double-stranded cDNA and ATAC fragments",MaterialKind.MOLECULAR,state=redistributed.state)
 g.react("second_strand","Make second-strand cDNA before lysis",StepKind.EXTENSION,((redistributed,InputRole.TEMPLATE),(second_primer,InputRole.STRUCTURAL),(second_pol,InputRole.CATALYST)),(ds,))
 rna=g.product("rna_branch","RNA-dedicated lysate",MaterialKind.MOLECULAR,state=ds.state); atac=g.product("atac_branch","ATAC-dedicated lysate",MaterialKind.MOLECULAR,state=ds.state)
 g.react("split","Lyse and split each well into RNA and ATAC branches",StepKind.SPLIT,((ds,InputRole.TEMPLATE),),(rna,atac),changes_molecule=False)
 rna_transposome=g.input("rna_transposome","unindexed cDNA transposome",MaterialKind.ADAPTER)
 rna_tn5=g.input("rna_tn5","cDNA Tn5 transposase",MaterialKind.CATALYST)
 tagged_rna=g.product("tagged_rna","tagmented cell-indexed cDNA",MaterialKind.MOLECULAR,state=RNA)
 rwaste=g.product("rna_tn5_waste","non-library cDNA tagmentation products",MaterialKind.WASTE)
 g.react("tagment_rna","Tagment second-strand cDNA",StepKind.TAGMENTATION,((rna,InputRole.TEMPLATE),(rna_transposome,InputRole.STRUCTURAL),(rna_tn5,InputRole.CATALYST)),(tagged_rna,rwaste)); g.discard(rwaste,"does not retain an amplifiable cDNA fragment")
 reviewed=["nucleus","rt_primer","rt","reverse_transcribe","atac_transposome","atac_tn5","tagment_atac","facs_wells","redistribute","second_primer","second_polymerase","second_strand","split","rna_transposome","rna_tn5","tagment_rna"]
 finish_pcr_branch(g,"rna",tagged_rna,RNA,RNA_P,reviewed,label="sci-CAR RNA",read_lengths=READ_LENGTHS["RNA library"])
 finish_pcr_branch(g,"atac",atac,ATAC,ATAC_P,reviewed,label="sci-CAR ATAC",read_lengths=READ_LENGTHS["ATAC library"])
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
