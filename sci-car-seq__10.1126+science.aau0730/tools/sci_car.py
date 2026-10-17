"""Cao et al. sci-CAR RNA/ATAC co-assay."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import joined_scene, seg
from chemdraw import Construct, Row, Segment, feature, revcomp
import illumina as il
import nextera as nx
import seqprimers as sp
from multimodal_spatial import finish_pcr_branch, modality_split
from protocol_graph import Disposition, Terminal
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

TITLE="sci-CAR-seq — joint combinatorial RNA and ATAC profiling"
NOTES="01_sci-car-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1126/science.aau0730">Cao et al., <i>Science</i> (2018)</a>.'
SUMMARY="One intact nucleus receives a barcoded poly(dT) RT primer and a separately barcoded Tn5 adapter in the first plate; FACS redistribution supplies a second well index while RNA and ATAC aliquots become separate libraries."
CAVEAT=("<b>Source scope and inferred adapter sequences.</b> Step order, reagents, conditions, the RT primer (Read 1 partial arm, 8-N UMI, 10-bp RT index, T30VN), the RNA and ATAC P5/P7 PCR primers, the cycle numbers and the RNA read recipe come from the published Supplementary Materials (Materials and Methods, &#8216;Sample processing for cultured cells&#8217;, pp. 16&#8211;18); the article&#8217;s Fig. 1 legend assigns each identifier to a read. The supplementary table workbook (Tables S1&#8211;S13) holds analysis results, not oligos. Dotted segments remain inferred: the mosaic end of the &#8216;i7 only TDE1&#8217; adaptor on the RNA cDNA end (the P7 primer establishes only S7), the 10-nt PCR-index lengths (only the RNA run states 10 index cycles), and the interior of the &#8216;indexed TDE1&#8217; ATAC transposons, which follows the Amini 2014 Table 4 architecture (connector, 8-nt index, extended mosaic end) because the article attributes transposase-complex preparation to that reference; the ATAC P5/P7 primers establish only the S5/S7 ends. ATAC libraries were read 51 + 51 with a custom recipe and &#8216;custom primers as previously described&#8217; (Cusanovich 2018, not read here); the ATAC sequencing primers drawn are therefore inferred and ATAC index-read cycles are not stated. "
        "RNA and ATAC first-round indexes are chemically distinct but map to the same first well; their second-round PCR indexes (carried by both PCR primers of each branch) likewise map to the same FACS well. "
        "<b>Variants not drawn:</b> PFA-fixed mouse-kidney nuclei are reverse-crosslinked (SDS, proteinase K, 65 &#176;C 16 h) after sorting and split before second-strand synthesis, which is then done on the AMPure-purified RNA portion only; &#8216;sci-RNA-seq only&#8217; and &#8216;sci-ATAC-seq only&#8217; control branches. "
        "<b>Published predecessor correction:</b> <a href=\"https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fng.3119/MediaObjects/41588_2014_BFng3119_MOESM2_ESM.xlsx\">Amini 2014 Supplementary Table 4</a>, Correction Notice and corrected rows 1–5, says the first five indexed Tn5 oligo rows were omitted from the original online table and restored on 20 May 2015. Their corrected 8-base i5 values are TATAGCCT, ATAGAGGC, CCTATCCT, GGCTCTGA and AGGCGAAG. This model uses generic index placeholders and does not assign those Amini members to sci-CAR; the inspected Cao 2018 supplementary materials print no barcode lists.")
RNA1=feature("rna_rt_well","cell_barcode","combinatorial",group="cell_id",part="RT well",whitelist="published sci-CAR RT primer plate")
ATAC1=feature("atac_tn5_well","cell_barcode","combinatorial",group="cell_id",part="Tn5 well (P5 side)",whitelist="published sci-CAR Tn5 plate")
UMI=feature("rna_umi","umi","random")
# Fig. 1 legend (PMC6571013 JATS line 475): RNA -- index 2 and read 1 cover the i5
# index, UMI and RT barcode; index 1 and read 2 cover the i7 index and cDNA.
# ATAC -- reads 1 and 2 cover genomic DNA; indices 1 and 2 cover the Tn5 and PCR
# barcodes.  Line 78: both PCR primers of each branch carry the second (PCR-well)
# index.  The published Supplementary Materials (NIHMS1033869-supplement-1.pdf,
# 'Sample processing for cultured cells', pp. 16-18; text twin lines 62-243) print
# the same primers, enzymes and cycles as the earlier-read GEO text, plus the steps
# GEO truncated (RNA/ATAC cleanups, ATAC sequencing recipe).  GEO GSE117089
# sample extraction protocol (author-deposited; identical in GSM3271040-45,
# geo-GSE117089-all-full.txt line 344488) prints the RT primer
# ACGACGCTCTTCCGATCT-N8-[10bp index]-T30VN, RNA P5 primer
# P5-[i5]-ACACTCTTTCCCTACACGACGCTCTTCCGATCT, ATAC P5 primer P5-[i5]-TCGTCGGCAGCGTC,
# P7 primer CAAGCAGAAGACGGCATACGAGAT-[i7]-GTCTCGTGGGCTCGG (both branches),
# 'i7 only TDE1' for cDNA tagmentation and the cycle numbers below.
GEO="GEO GSE117089 sample extraction protocol (author-deposited; GSM3271040-GSM3271045)"
SM="Cao 2018 Supplementary Materials, Materials and Methods, 'Sample processing for cultured cells' (NIHMS1033869-supplement-1.pdf pp. 16-18)"
AMINI = "Amini 2014 Supplementary Table 4 (doi:10.1038/ng.3119), cited by Cao 2018 line 140 for transposase-complex preparation"
# Amini Table 4 transposon architecture: S5 + connector A + 8-nt index + extension A + ME;
# S7 + connector B + 8-nt index + extension B + ME.
CONNECTOR_A="TCCACGC"; EXT_A="GCGATCGAGGACGGC"
CONNECTOR_B="CTGTCCCTGTCC"; EXT_B="CACCGTCTCCGCCTC"
ATAC_R1=EXT_A+nx.ME
ATAC_R2=EXT_B+nx.ME
ATAC_PRIMERS=(
 sp.custom("Read 1","extended mosaic end A (inferred)",ATAC_R1,AMINI),
 sp.custom("Index 1 (i7)","reverse complement of extended mosaic end B (inferred)",revcomp(ATAC_R2),AMINI),
 sp.custom("Index 2 (i5)","reverse complement of extended mosaic end A (derived, inferred)",revcomp(ATAC_R1),AMINI+"; derived as the reverse complement of Read 1"),
 sp.custom("Read 2","extended mosaic end B (inferred)",ATAC_R2,AMINI),
)

def _pcr_index(name,feat):
 return seg(name,"J"*10,"cbc",placeholder=True,inferred=True,feature=feat)

def rna_library():
 i5=feature("rna_pcr_i5","cell_barcode","combinatorial",group="cell_id",part="FACS/PCR well (i5)",whitelist="published sci-CAR PCR plate")
 i7=feature("rna_pcr_i7","cell_barcode","combinatorial",group="cell_id",part="FACS/PCR well (i7)",whitelist="published sci-CAR PCR plate")
 lib=Construct([seg("P5",il.P5,"p5"),_pcr_index("RNA PCR i5",i5),
  seg("TruSeq Read 1 arm (P5 primer tail + RT primer 18 nt)",il.TRUSEQ_READ1,"r1"),
  seg("UMI","U"*8,"umi",placeholder=True,feature=UMI),
  seg("RT well index","A"*10,"cbc",placeholder=True,feature=RNA1),
  seg("oligo(dT) and 3-prime cDNA tag","X"*38,placeholder=True),
  seg("unindexed Tn5 mosaic end",nx.ME_RC,"me",inferred=True),
  seg("Tn5 adaptor S7 (P7 primer 3' end)",nx.S7_RC,"s7"),
  Segment("RNA PCR i7'","I"*10,"cbc",True,inferred=True,feature=i7),
  seg("P7'",il.P7_RC,"p7")],name="sci-CAR RNA library")
 primers=(sp.TRUSEQ["R1"],sp.NEXTERA["I1"],sp.TRUSEQ["I2"],sp.NEXTERA["R2"])
 errors=sp.verify(lib,primers)
 if errors: raise ValueError("invalid sci-CAR RNA library: "+"; ".join(errors))
 return lib,primers

def atac_library():
 i5=feature("atac_pcr_i5","cell_barcode","combinatorial",group="cell_id",part="FACS/PCR well (i5)",whitelist="published sci-CAR PCR plate")
 i7=feature("atac_pcr_i7","cell_barcode","combinatorial",group="cell_id",part="FACS/PCR well (i7)",whitelist="published sci-CAR PCR plate")
 tn5_b=feature("atac_tn5_well_b","cell_barcode","combinatorial",group="cell_id",part="Tn5 well (P7 side)",whitelist="published sci-CAR Tn5 plate")
 lib=Construct([seg("P5",il.P5,"p5"),_pcr_index("ATAC PCR i5",i5),
  seg("S5 (ATAC P5 primer 3' end)",nx.S5,"s5"),seg("connector A",CONNECTOR_A,"r2",inferred=True),
  seg("Tn5 well index (P5 side)","A"*8,"cbc",placeholder=True,inferred=True,feature=ATAC1),
  seg("extension A",EXT_A,"r1",inferred=True),seg("ME",nx.ME,"me",inferred=True),
  seg("accessible DNA","X"*38,placeholder=True),
  seg("ME'",nx.ME_RC,"me",inferred=True),seg("extension B'",revcomp(EXT_B),"r1",inferred=True),
  Segment("Tn5 well index (P7 side)'","B"*8,"cbc",True,inferred=True,feature=tn5_b),
  seg("connector B'",revcomp(CONNECTOR_B),"r3",inferred=True),seg("S7' (P7 primer 3' end)",nx.S7_RC,"s7"),
  Segment("ATAC PCR i7'","I"*10,"cbc",True,inferred=True,feature=i7),
  seg("P7'",il.P7_RC,"p7")],name="sci-CAR ATAC library")
 errors=sp.verify(lib,ATAC_PRIMERS)
 if errors: raise ValueError("invalid sci-CAR ATAC library: "+"; ".join(errors))
 return lib,ATAC_PRIMERS

RNA,RNA_P=rna_library()
ATAC,ATAC_P=atac_library()
FINAL_LIBRARIES=(("RNA library",RNA,RNA_P,"INFERRED — the cDNA-end mosaic end (vendor 'i7 only TDE1' adaptor) and the 10-nt PCR-index lengths; P5/P7 primers, Read 1 arm, 8-N UMI and 10-bp RT index follow the published Supplementary Materials (pp. 16-17).","Fig. 1 legend: Read 1 covers the UMI and RT barcode, Index 2 the i5 PCR index, Index 1 the i7 PCR index and Read 2 the cDNA."),("ATAC library",ATAC,ATAC_P,"INFERRED — transposon interior (connector, Tn5 index, extended mosaic end) follows Amini 2014 Table 4 (the cited transposase-complex method); the supplementary materials call it 'indexed TDE1 enzyme (Illumina)' without bases. S5/S7 ends follow the printed PCR primers (Supplementary Materials p. 18).","Fig. 1 legend: Reads 1 and 2 cover genomic DNA; Index 1 and Index 2 each cover a Tn5 barcode and a PCR barcode."))
# Supplementary Materials p. 17 (and GEO GSE117089): RNA "Read 1: 18 cycles, Read 2: 52 cycles, Index 1: 10
# cycles, Index 2: 10 cycles" (NextSeq V2 75-cycle kit); ATAC p. 18 "paired end 51 bp
# reads with index reads that covered both the Tn5 barcode and the library amplification
# barcode", custom recipe and primers from ref 29 (index cycles not stated).
READ_LENGTHS={"RNA library":{"Read 1":18,"Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":52},"ATAC library":{"Read 1":51,"Read 2":51}}

def sections(): return [
 ("Index RNA by in-situ reverse transcription",[Row(chunks=[("poly(dT) -- UMI -- RT-well barcode ---> poly(A) RNA", "cbc",False)])],"A well-specific anchored oligo(dT) primer (Read 1 partial arm, 8-N UMI, 10-bp index, T30VN; SuperScript IV, 4&#8211;55 &#176;C gradient) copies mRNA and installs the first RNA index plus UMI."),
 ("Index accessible DNA by in-situ tagmentation",joined_scene([seg("left barcoded Tn5 adapter","B"*18,"cbc",placeholder=True),seg("accessible genomic DNA","X"*34,placeholder=True),seg("right barcoded Tn5 adapter","B"*18,"cbc",placeholder=True)],(("left barcoded Tn5 adapter","accessible genomic DNA","Tn5 transfer"),("accessible genomic DNA","right barcoded Tn5 adapter","Tn5 transfer")),label="barcoded accessible-DNA fragment").rows(),"A separately barcoded transposome installs the first ATAC index in the same intact nuclei."),
 ("Pool nuclei and redistribute by FACS",[Row(chunks=[("RT/Tn5 well identity + FACS well identity = cell identity", "cbc",False)])],"Most nuclei traverse a unique pair of wells; matched well coordinates link modalities."),
 ("Make second-strand cDNA, purify and split each well",modality_split("RNA-dedicated lysate","ATAC-dedicated lysate"),"After second-strand synthesis (16 &#176;C, 3 h) each well is stopped with DNA binding buffer, purified on AMPure XP beads and the eluate split 6 &#181;L + 6 &#181;L into an RNA plate and an ATAC plate."),
 ("Complete the two indexed libraries",modality_split("unindexed Tn5 on cDNA + RNA-index PCR","ATAC-index PCR"),"The second well-specific index is added by modality-specific PCR before separate pooling and sequencing."),
]

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 nucleus=g.input("nucleus","intact nucleus containing accessible chromatin and poly(A) RNA",MaterialKind.BIOLOGICAL,state=molecule("intact nucleus",(role("accessible chromatin"),role("poly(A) RNA"))),primary=True)
 rtprimer=g.input("rt_primer","well-indexed anchored oligo(dT) RT primer (Read 1 partial arm, 8-N UMI, 10-bp index, T30VN)",MaterialKind.OLIGO,source=SM,
  state=Construct([seg("TruSeq Read 1 partial arm","ACGACGCTCTTCCGATCT","r1"),seg("UMI","N"*8,"umi",placeholder=True,feature=UMI),
   seg("RT well index","A"*10,"cbc",placeholder=True,feature=RNA1),seg("oligo(dT)30","T"*30),seg("VN anchor","VN",placeholder=True)],name="sci-CAR RT primer"))
 rt=g.input("rt","SuperScript IV reverse transcriptase with RNaseOUT (primer/dNTP 55 °C 5 min, then 4→55 °C gradient)",MaterialKind.CATALYST,source=SM)
 copied=g.product("copied","nucleus retaining first-well-indexed cDNA and chromatin",MaterialKind.MOLECULAR,state=molecule("RT-indexed nucleus",(role("RT-well barcode"),role("UMI"),role("first-strand cDNA"),role("accessible chromatin"))))
 g.react("reverse_transcribe","Index RNA by in-situ reverse transcription",StepKind.REVERSE_TRANSCRIPTION,((nucleus,InputRole.TEMPLATE),(rtprimer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(copied,))
 transposome=g.input("atac_transposome","well-indexed TDE1 (Illumina) in Nextera TD buffer, 55 °C 30 min, stopped with EDTA/spermidine; transposon bases not printed (Amini 2014 architecture, inferred)",MaterialKind.ADAPTER,inferred=True,source=AMINI)
 tn5=g.input("atac_tn5","Tn5 transposase",MaterialKind.CATALYST)
 tagged=g.product("tagged","nucleus carrying indexed cDNA and tagmented chromatin",MaterialKind.MOLECULAR,state=copied.state)
 twaste=g.product("tn5_waste","candidate non-library tagmentation products",MaterialKind.MOLECULAR)
 g.react("tagment_atac","Index accessible DNA by in-situ tagmentation",StepKind.TAGMENTATION,((copied,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,twaste)); g.mark_terminal(Terminal(twaste.id, Disposition.UNRESOLVED, 'Candidate non-library end classes remain with the nuclei and enter ATAC PCR. They lack the intended P5/P7 end combination; amplification and physical recovery are unresolved.'))
 facs=g.input("facs_wells","pool, DAPI stain (3 µM), singlet-gated FACS of 25 nuclei per well into 5 µL EB",MaterialKind.CONDITION,source=SM)
 redistributed=g.product("redistributed","uniquely redistributed indexed nuclei",MaterialKind.MOLECULAR,state=tagged.state)
 excluded=g.product("excluded","DAPI-gated doublets and unsorted nuclei",MaterialKind.WASTE)
 g.react("redistribute","Pool nuclei and redistribute them by FACS",StepKind.PARTITION,((tagged,InputRole.TEMPLATE),(facs,InputRole.CONDITION)),(redistributed,excluded),changes_molecule=False); g.discard(excluded,"physically excluded by the sorter gate; barcode collisions between separately sorted nuclei remain in the wells and can only be recognised computationally")
 second_pol=g.input("second_polymerase","mRNA Second Strand Synthesis buffer and enzyme (NEB), 16 °C 180 min; no added primer is printed",MaterialKind.CATALYST,source=SM)
 ds=g.product("double_stranded","nuclei containing double-stranded cDNA and ATAC fragments",MaterialKind.MOLECULAR,state=redistributed.state)
 g.react("second_strand","Make second-strand cDNA before lysis",StepKind.EXTENSION,((redistributed,InputRole.TEMPLATE),(second_pol,InputRole.CATALYST)),(ds,))
 rna=g.product("rna_branch","RNA-dedicated lysate",MaterialKind.MOLECULAR,state=ds.state); atac=g.product("atac_branch","ATAC-dedicated lysate",MaterialKind.MOLECULAR,state=ds.state)
 cleanup=g.input("cleanup","6 µL Zymo DNA binding buffer (5 min), 24 µL AMPure XP, elute 12.5 µL EB",MaterialKind.CONDITION,source=SM)
 cwaste=g.product("cleanup_waste","bead supernatant and washes (unbound material)",MaterialKind.WASTE)
 g.react("split","Stop with DNA binding buffer, AMPure-purify and split each eluate (6 µL) into RNA and ATAC plates",StepKind.SPLIT,((ds,InputRole.TEMPLATE),(cleanup,InputRole.CONDITION)),(rna,atac,cwaste),changes_molecule=False,source=SM)
 g.discard(cwaste,"removed with the bead supernatant and washes")
 rna_transposome=g.input("rna_transposome","unindexed cDNA transposome ('i7 only TDE1', 25 nM, Illumina) in Nextera TD buffer, 55 °C 5 min; adaptor bases not in inspected sources",MaterialKind.ADAPTER,source=SM)
 rna_tn5=g.input("rna_tn5","cDNA Tn5 transposase",MaterialKind.CATALYST)
 tagged_rna=g.product("tagged_rna","tagmented cell-indexed cDNA",MaterialKind.MOLECULAR,state=RNA)
 rwaste=g.product("rna_tn5_waste","candidate non-library cDNA tagmentation products",MaterialKind.MOLECULAR)
 g.react("tagment_rna","Tagment second-strand cDNA",StepKind.TAGMENTATION,((rna,InputRole.TEMPLATE),(rna_transposome,InputRole.STRUCTURAL),(rna_tn5,InputRole.CATALYST)),(tagged_rna,rwaste)); g.mark_terminal(Terminal(rwaste.id, Disposition.UNRESOLVED, 'Candidate RNA-side products lacking the RT-primer end remain in the RNA PCR pool. The end-class calculation does not establish absence of amplification or physical removal.'))
 reviewed=["nucleus","rt_primer","rt","reverse_transcribe","atac_transposome","atac_tn5","tagment_atac","facs_wells","redistribute","second_polymerase","second_strand","cleanup","split","rna_transposome","rna_tn5","tagment_rna"]
 finish_pcr_branch(g,"rna",tagged_rna,RNA,RNA_P,reviewed,label="sci-CAR RNA",read_lengths=READ_LENGTHS["RNA library"],
  step_label="Release Tn5 with SDS/BSA and P7 primer (55 °C 15 min), quench with Tween-20, add indexed P5 primer and PCR 18–22 cycles (66 °C anneal); pool, 0.8× AMPure",
  forward_label="RNA indexed P5 primer P5-[i5]-ACACTCTTTCCCTACACGACGCTCTTCCGATCT (SM p. 17)",
  reverse_label="indexed P7 primer CAAGCAGAAGACGGCATACGAGAT-[i7]-GTCTCGTGGGCTCGG (SM p. 17)",
  polymerase_label="NEBNext High-Fidelity 2X PCR Master Mix",primer_source=SM,
  limitation="The printed primers carry [i5]/[i7] placeholders; extension across the Tn5 gap is not drawn as a separate state.")
 finish_pcr_branch(g,"atac",atac,ATAC,ATAC_P,reviewed,label="sci-CAR ATAC",read_lengths=READ_LENGTHS["ATAC library"],
  step_label="PCR with indexed P5 and P7 primers, 72 °C 3 min then 18–22 cycles (63 °C anneal); pool, 1× AMPure, gel-extract >200 bp",
  forward_label="ATAC indexed P5 primer P5-[i5]-TCGTCGGCAGCGTC (SM p. 18)",
  reverse_label="indexed P7 primer CAAGCAGAAGACGGCATACGAGAT-[i7]-GTCTCGTGGGCTCGG (SM p. 18)",
  polymerase_label="NEBNext High-Fidelity 2X PCR Master Mix",primer_source=SM,
  limitation="The printed primers carry [i5]/[i7] placeholders; the 72 °C gap fill of the Tn5 ends is not drawn as a separate state.")
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
