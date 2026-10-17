"""GAGE-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct,MolecularState,Row,Scene,Segment,Workflow,feature,oligo,revcomp
from multimodal_spatial import HybridizationLigationRound,finish_pcr_branch
from pcr import PCRGeometry,PCRPrimer
from protocol_extensions import proximity_product
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,StepKind
CELL1=feature("gage_cell_round1","cell_barcode","combinatorial",group="gage_cell",part="ligation round 1",whitelist="96-well adapter set")
CELL2=feature("gage_cell_round2","cell_barcode","combinatorial",group="gage_cell",part="ligation round 2",whitelist="96-well adapter set")
UMI_F=feature("gage_rna_umi","umi","random")
def payload(label,rna=False):
 p=[seg("cell barcode round 1","B"*8,"cbc",placeholder=True,feature=CELL1),seg("cell barcode round 2","C"*8,"cbc",placeholder=True,feature=CELL2)]
 if rna:p.append(seg("RT UMI","U"*10,"umi",placeholder=True,feature=UMI_F))
 p.append(seg(label,"X"*42,placeholder=True));return p

RT_HANDLE="AGACAGGACACTGAC"
RNA_LINKER1="GTCAGTGTCCTGTCTCGACTGACTGACCAG"
RNA_R1_LEFT="CAACGAGGTACGACA"; RNA_R1_BC="CTGATCGT"; RNA_R1_RIGHT="CTGGTCAGTCAGTCG"
RNA_LINKER2_RNA="TGUCGUACCUCGTUGCUCAAGCACGUGGAT"
RNA_LINKER2=RNA_LINKER2_RNA.replace("U","T")
RNA_R2_LEFT="ACGACTCACTATAGGGACGACGCTCTTCCGATCT"
RNA_R2_BC="AACGTGAT"; RNA_R2_RIGHT="ATCCACGTGCTTGAG"
RNA_ROUND1=HybridizationLigationRound("RNA round 1",RT_HANDLE,RNA_R1_LEFT,RNA_R1_RIGHT,RNA_LINKER1,CELL1,phosphorylated_5=True)
RNA_ROUND2=HybridizationLigationRound("RNA round 2",RNA_R1_LEFT,RNA_R2_LEFT,RNA_R2_RIGHT,RNA_LINKER2,CELL2)

# Table 1a prints /iBiodT/ as a base between the 15-nt handle and TC G: a biotin-dT residue.
RT_BIOTIN_DT="T"
RT_DT=RT_HANDLE+RT_BIOTIN_DT+"TCG"+"T"*20+"VN"
RT_RANDOM=RT_HANDLE+RT_BIOTIN_DT+"TCG"+"N"*6
POST_TDT_POLYC="AGAGTGAAC"+"C"*12
POST_TDT_P5_T7="AATTTAATACGACTCACTATAGGGA"
# Supplementary Tables 1e/1f/1g, well A1 (AD-1_F_1 / AD-1_R_1).
DNA_AD1_BC="AGTACACG"
DNA_AD1_LINK="TGACTTGTAGTGCAG"
DNA_AD1_F_TNA="TNAA"+DNA_AD1_BC+nx.ME_RC+"GA"
DNA_AD1_F_TA="TAA"+DNA_AD1_BC+nx.ME_RC+"GA"
DNA_AD1_R=DNA_AD1_LINK+revcomp(DNA_AD1_BC)+"T"
NEXTERA_P5_AD_F_TNA="TNA"+nx.ME_RC+"GA"; NEXTERA_P5_AD_F_TA="TA"+nx.ME_RC+"GA"; NEXTERA_P5_AD_R="GTATAAGAGACAG"
DNA_AD2_LO="CAAGTCA"+"TATCAGCG"+"AGATCGG"
DNA_AD2_HI="GGAGTTCAGACGTGTGCTCTTCCGATCT"+"CGCTGATA"
DNA_PRE_F=nx.ADAPTOR_S5
DNA_PRE_R=il.TRUSEQ_READ2[:24]

DNA_I5="TCAATATCGA"; DNA_I7="CAAGGTCA"
RNA_I5="CATGCGAGAAAC"; RNA_I7="CTCTTAGCGG"
DNA_FINAL_F=il.P5+DNA_I5+nx.S5
DNA_FINAL_R=il.P7+DNA_I7+il.TRUSEQ_READ2[:24]
RNA_FINAL_F=il.P5+RNA_I5+il.TRUSEQ_READ1
RNA_FINAL_R=il.P7+RNA_I7+nx.S7

DNA=Construct([seg("P5",il.P5,"p5"),seg("i5",DNA_I5,"cbc",feature=feature("gage_dna_i5","sample_index","whitelist",whitelist="Supplementary Table 1b Nextera-P5 series")),seg("S5",nx.S5,"s5"),seg("left ME",nx.ME,"me"),seg("genomic contact insert ending in the TNA/TA restriction end","X"*42,placeholder=True),seg("AD-1 spacer","A"),seg("cell barcode round 1","B"*8,"cbc",placeholder=True,feature=CELL1),seg("AD-1/AD-2 joint (reverse complement)",revcomp(DNA_AD1_LINK)),seg("cell barcode round 2","C"*8,"cbc",placeholder=True,feature=CELL2),seg("TruSeq Read 2 arm reverse complement",revcomp(il.TRUSEQ_READ2),"r2"),seg("i7 reverse complement",revcomp(DNA_I7),"cbc",feature=feature("gage_dna_i7","sample_index","whitelist",whitelist="Supplementary Table 1b TruSeq-8bp series")),seg("P7 reverse complement",il.P7_RC,"p7")],name="GAGE-seq chromatin-contact library")
DNA_P=(sp.NEXTERA["R1"],sp.TRUSEQ["I1"],sp.NEXTERA["I2"],sp.TRUSEQ["R2"])
RNA=Construct([seg("P5",il.P5,"p5"),seg("i5",RNA_I5,"cbc",feature=feature("gage_rna_i5","sample_index","whitelist",whitelist="Supplementary Table 1b TruSeq-F-L series")),seg("TruSeq Read 1 arm",il.TRUSEQ_READ1,"r1"),seg("Round-2 N10 UMI","U"*10,"umi",placeholder=True,feature=UMI_F),seg("cell barcode round 2","C"*8,"cbc",placeholder=True,feature=CELL2),seg("round-2/round-1 linker arms",RNA_R2_RIGHT+RNA_R1_LEFT),seg("cell barcode round 1","B"*8,"cbc",placeholder=True,feature=CELL1),seg("round-1 arm and RT handle",RNA_R1_RIGHT+RT_HANDLE),seg("biotin-dT",RT_BIOTIN_DT),seg("RT primer spacer","TCG"),seg("cDNA from dT20VN or N6 priming","X"*42,placeholder=True),seg("right ME reverse complement",nx.ME_RC,"me"),seg("S7 reverse complement",nx.S7_RC,"s7"),seg("i7 reverse complement",revcomp(RNA_I7),"cbc",feature=feature("gage_rna_i7","sample_index","whitelist",whitelist="Supplementary Table 1b Nextera-P7 series")),seg("P7 reverse complement",il.P7_RC,"p7")],name="GAGE-seq RNA library")
RNA_P=(sp.TRUSEQ["R1"],sp.NEXTERA["I1"],sp.TRUSEQ["I2"],sp.NEXTERA["R2"])
for _lib,_primers in ((DNA,DNA_P),(RNA,RNA_P)):
 problems=sp.verify(_lib,_primers,required_roles=tuple(p.role for p in _primers))
 if problems: raise ValueError(_lib.name+": "+"; ".join(problems))

DNA_FINAL_GEOMETRY=PCRGeometry.exact(DNA,PCRPrimer("Nextera-P5-2A",DNA_FINAL_F),PCRPrimer("TruSeq-8bp-R1",DNA_FINAL_R,anneal_nt=len(il.TRUSEQ_READ2[:24])),statement="Published Supplementary Table 1b prints the complete representative hybrid Nextera-P5/TruSeq-P7 primer pair.")
RNA_FINAL_GEOMETRY=PCRGeometry.exact(RNA,PCRPrimer("TruSeq-F-L-1",RNA_FINAL_F),PCRPrimer("Nextera-P7-6A",RNA_FINAL_R),statement="Published Supplementary Table 1b prints the complete representative hybrid TruSeq-P5/Nextera-P7 primer pair.")
DNA_PRE_GEOMETRY=PCRGeometry.exact(DNA,PCRPrimer("Nextera-P5-pre-primer",DNA_PRE_F),PCRPrimer("TruSeq-P7-pre-P-S",DNA_PRE_R,anneal_nt=len(DNA_PRE_R)),statement="Published Supplementary Table 1a prints both Hi-C pre-amplification primers.")
TITLE="GAGE-seq — gene expression and genome architecture"
NOTES="01_gage-seq.html"
SOURCE='<a href="https://doi.org/10.1038/s41588-024-01745-3">Zhou et al., <i>Nature Genetics</i> (2024)</a> (protocol first posted as bioRxiv 10.1101/2023.07.20.549578).'
SUMMARY="In fixed cells, RNA is reverse-transcribed and genomic contacts are proximity-ligated. Two rounds of well-specific adapter ligation barcode both analytes before cDNA/DNA separation and library-specific tagmentation."
CAVEAT=("The published Supplementary Table 1 gives every 96-well RNA and Hi-C adapter sequence plus the non-barcoded and terminal PCR primers. The page shows A1 representatives while retaining the two 96-member barcode whitelists; B/C/U/N/X denote varying barcode, UMI, random or biological sequence. "
 "Model corrections in this pass (2026-10-17): the published Methods order is RT before MseI/CviQI digestion, proximity ligation without fill-in, and DdeI re-digestion; the TdT tail is poly(G) (dGTP/ddGTP), annealed by the printed poly(C)12 primer; the Hi-C AD-1 A1 barcode and 3-prime GA, the /iBiodT/ residue of the RT primers and the Hi-C read-2 barcode order were mis-transcribed earlier. "
 "Open: the Methods trim 12 nt from DNA read 1, which matches an AD-1 BC1+A+TNA end rather than a Tn5 end; the model draws the Tn5 end and does not resolve which P5 end read 1 starts from.")
FINAL_LIBRARIES=(("Chromatin-contact library",DNA,DNA_P,"Proximity-ligated DNA is capped with Nextera-P5-AD, preamplified, tagmented and indexed.","Read 2 reads BC2, the 15-nt AD-1/AD-2 joint, BC1 and the restriction end (Methods trim 35 nt); read 1 reads genomic sequence."),("RNA library",RNA,RNA_P,"Biotinylated cDNA is recovered on streptavidin beads, poly(G)-tailed, amplified and tagmented.","Read 1 reads the Round-2 N10 UMI, BC2, the 30-nt linker, BC1 and the 15-nt round-1 arm (Methods positions 1-71)."))
READ_LENGTHS={"Chromatin-contact library":{"Read 1":150,"Read 2":150},"RNA library":{"Read 1":150,"Read 2":150}}

def hic_ad1_scene():
 sc=Scene()
 sc.strand("AD-1_F_1",[Segment("AD-1 F TNA overhang","TNA"),Segment("AD-1 F duplex","A"+DNA_AD1_BC+"CTG"),Segment("AD-1 F Nextera arm",nx.ME_RC[3:]+"GA","me")],label="AD-1_F_1 (Table 1e A1)",mod5="/5Phos/")
 sc.anneal("AD-1_R_1",[Segment("AD-1 R AD-2 joint arm",DNA_AD1_LINK[:12]),Segment("AD-1 R duplex",DNA_AD1_LINK[12:]+revcomp(DNA_AD1_BC)+"T")],to="AD-1_F_1",pair=("AD-1 R duplex","AD-1 F duplex"),unpaired=("AD-1 R AD-2 joint arm",),label="AD-1_R_1 (Table 1g A1)",mod5="/5Phos/")
 return sc
def hic_ad2_scene():
 sc=Scene()
 sc.strand("SeqShort_BC1_LO",[Segment("AD-2 LO joint overhang","CAAGTCA"),Segment("AD-2 LO duplex",DNA_AD2_LO[7:])],label="SeqShort_BC1_LO (Table 1h A1)")
 sc.anneal("SeqShort_BC1_HI",[Segment("AD-2 HI TruSeq arm",DNA_AD2_HI[:21],"r2"),Segment("AD-2 HI duplex",DNA_AD2_HI[21:])],to="SeqShort_BC1_LO",pair=("AD-2 HI duplex","AD-2 LO duplex"),label="SeqShort_BC1_HI (Table 1i A1)")
 return sc
def barcode_geometry(number=None):
 rounds={1:(RNA_ROUND1.adapter_scene(),hic_ad1_scene()),2:(RNA_ROUND2.adapter_scene(),hic_ad2_scene())}
 scenes=rounds[number] if number else (*rounds[1],*rounds[2])
 return AdapterLigationGeometry.known("splinted",*scenes,statement="Published Supplementary Tables 1c-1i print both strands of the A1 RNA splinted adapters and Hi-C AD-1 (5-prime TNA/TA cohesive end, Y-shaped) and AD-2 (CAAGTCA overhang matching the AD-1_R TGACTTG arm) duplexes; the shared constructor verifies the RNA splint arms and the Scene verifies duplex pairing. The Round2 oligo also carries an N10 UMI between its TruSeq arm and BC2, not drawn in the splint panel.")

def oligos():
 for name,seq,mods in (("BiodT-dT(20)VN RT primer",RT_DT,"/5Phos/; base 16 is /iBiodT/"),("BiodT random-hexamer RT primer",RT_RANDOM,"/5Phos/; base 16 is /iBiodT/"),("RNA-ligation linker-1",RNA_LINKER1,""),("RNA-ligation linker-2",RNA_LINKER2_RNA,"RNA"),("RNA Ligation-1-blocker","CGAAGTCAGTGTCCTGTC",""),("Round1_01",RNA_R1_LEFT+RNA_R1_BC+RNA_R1_RIGHT,"/5Phos/"),("Round2_01",RNA_R2_LEFT+"N"*10+RNA_R2_BC+RNA_R2_RIGHT,""),("Post-TdT-PolyC12-S",POST_TDT_POLYC,""),("Post-TdT-P5(T7)-primer",POST_TDT_P5_T7,""),("Hi-C AD-1_F_1 (TNA)",DNA_AD1_F_TNA,"/5Phos/"),("Hi-C AD-1_R_1",DNA_AD1_R,"/5Phos/"),("Hi-C AD2 low A1",DNA_AD2_LO,""),("Hi-C AD2 high A1",DNA_AD2_HI,""),("Nextera-P5-AD-F-TNA",NEXTERA_P5_AD_F_TNA,"/5Phos/"),("Nextera-P5-AD-F-TA",NEXTERA_P5_AD_F_TA,"/5Phos/"),("Nextera-P5-AD-R",NEXTERA_P5_AD_R,""),("Nextera-P5-pre-primer",DNA_PRE_F,""),("TruSeq-P7-pre-P-S",DNA_PRE_R,""),("Nextera-P5-2A",DNA_FINAL_F,""),("TruSeq-8bp-R1",DNA_FINAL_R,""),("TruSeq-F-L-1",RNA_FINAL_F,""),("Nextera-P7-6A",RNA_FINAL_R,"")):
  yield oligo(name,[Segment(name,seq,placeholder=any(x in seq for x in "N"))],mods=mods)

def sections():return [("Reverse-transcribe RNA in lysed fixed cells",[Row(chunks=[("5′-phosphorylated /iBiodT/ poly(dT)VN or random-hexamer RT primer → cDNA (Maxima H Minus)","umi",False)])],"Published Methods: SDS-treated cells are reverse-transcribed before chromatin digestion; the printed RT primers provide the 15-base handle used by RNA ligation round 1."),("Digest, ligate and re-digest chromatin",proximity_product(left="locus A",right="locus B",bridge="MseI/CviQI TA-overhang junction").rows(),"MseI + CviQI digestion, T4 DNA ligase proximity ligation (16 °C overnight, no fill-in step in the Methods), then DdeI re-digestion."),("Ligate two cell indexes",barcode_geometry().rows(),"The A1 RNA splinted adapters and Hi-C AD-1/AD-2 duplexes are shown exactly; the source prints all 96 members of each plate."),("Separate and build libraries",[Row(chunks=[("streptavidin cDNA ──► poly(G) tail ──► TruSeq-P5 / Nextera-P7 RNA library\nsupernatant DNA  ──► Nextera-P5-AD ligation ──► Nextera-P5 / TruSeq-P7 contact library",None,False)])],"Each branch is preamplified, tagmented and completed with its source-printed hybrid primer pair for PE150.")]
def workflow():
 s=sections();w=Workflow(MolecularState(s[0][0],tuple(s[0][1])));w.react(s[1][0],s[1][1],note=s[1][2]);w.ligate_adapter(s[2][0],barcode_geometry(),s[2][1],note=s[2][2]);return w


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input("cells","formaldehyde-crosslinked cells, lysed and SDS-treated",MaterialKind.BIOLOGICAL,state=sections()[0][1],primary=True)
 rt_primer=g.input("rt_primer","5-prime-phosphorylated /iBiodT/ dT20VN and random-hexamer RT primers (Supplementary Table 1a)",MaterialKind.OLIGO)
 rt=g.input("rt","Maxima H Minus reverse transcriptase",MaterialKind.CATALYST)
 copied=g.product("copied","cells containing handle-tagged biotinylated cDNA",MaterialKind.MOLECULAR,state=sections()[0][1])
 g.react("reverse_transcribe","Reverse-transcribe RNA in situ with the biotinylated handle primers",StepKind.REVERSE_TRANSCRIPTION,((cells,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(copied,))
 nuclease=g.input("restriction_nuclease","MseI and CviQI",MaterialKind.CATALYST)
 digested=g.product("digested","MseI/CviQI-digested crosslinked chromatin (5-prime TA overhangs) with cDNA",MaterialKind.MOLECULAR,state=copied.state)
 g.react("digest","Digest crosslinked chromatin in situ (25 °C, 2 h)",StepKind.DIGESTION,((copied,InputRole.TEMPLATE),(nuclease,InputRole.CATALYST)),(digested,))
 contact_ligase=g.input("contact_ligase","T4 DNA ligase",MaterialKind.CATALYST)
 contacts=g.product("contacts","proximity-ligated chromatin contacts with cDNA",MaterialKind.MOLECULAR,state=sections()[1][1])
 g.react("proximity_ligate","Ligate spatially adjacent TA-overhang ends in situ (16 °C overnight)",StepKind.PROXIMITY_LIGATION,((digested,InputRole.TEMPLATE),(contact_ligase,InputRole.CATALYST)),(contacts,))
 dde=g.input("second_nuclease","DdeI",MaterialKind.CATALYST)
 redigested=g.product("redigested","DdeI-redigested contacts with TNA and TA 5-prime overhang ends",MaterialKind.MOLECULAR,state=contacts.state)
 g.react("redigest","Second-round chromatin fragmentation with DdeI (37 °C, 1 h)",StepKind.DIGESTION,((contacts,InputRole.TEMPLATE),(dde,InputRole.CATALYST)),(redigested,))

 barcode_ligase=g.input("barcode_ligase","T4 DNA ligase",MaterialKind.CATALYST)
 current=redigested
 reviewed=["cells","rt_primer","rt","reverse_transcribe","restriction_nuclease","digest","contact_ligase","proximity_ligate","second_nuclease","redigest","barcode_ligase"]
 labels={1:"in-cell round 1: RNA Round1 + linker-1 and Hi-C AD-1 (TNA and TA forward plates with reverse plate)",2:"in-well round 2 after 60 °C SDS lysis: RNA Round2 + linker-2, Hi-C AD-2, non-barcoded Nextera-P5-AD and RNA ligation-1 blocker"}
 for number in (1,2):
  adapter=g.input(f"barcode_{number}",f"{labels[number]}",MaterialKind.ADAPTER)
  indexed=g.product(f"indexed_{number}",f"DNA and cDNA after cell-barcode round {number}",MaterialKind.MOLECULAR,state=sections()[2][1])
  g.react(f"ligate_{number}",f"Ligate cell-barcode round {number}",StepKind.LIGATION,((current,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(barcode_ligase,InputRole.CATALYST)),(indexed,),geometry=barcode_geometry(number))
  reviewed.extend((adapter.id,f"ligate_{number}"));current=indexed

 protease=g.input("proteinase_k","proteinase K with SDS; 60 °C reverse crosslinking",MaterialKind.CATALYST)
 released=g.product("released","deproteinized gDNA and cDNA, isopropanol-precipitated",MaterialKind.MOLECULAR,state=current.state)
 g.react("reverse_crosslink","Proteinase K digestion and reverse crosslinking, then pooled isopropanol precipitation",StepKind.DEPROTEINIZATION,((current,InputRole.TEMPLATE),(protease,InputRole.CATALYST)),(released,))
 streptavidin=g.input("streptavidin","MyOne C1 streptavidin Dynabeads",MaterialKind.SELECTION)
 rna=g.product("rna_branch","bead-bound biotinylated cDNA branch",MaterialKind.MOLECULAR,state=payload("cDNA insert",True))
 dna=g.product("dna_branch","supernatant genomic contact DNA branch",MaterialKind.MOLECULAR,state=payload("proximity-ligated genomic junction"))
 g.react("separate","Separate biotinylated cDNA from genomic contact DNA",StepKind.AFFINITY_CAPTURE,((released,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),(rna,dna),changes_molecule=False)

 tdt=g.input("tdt","terminal transferase with CoCl2 and 0.95 µM dGTP / 0.05 µM ddGTP",MaterialKind.CATALYST)
 rna_tailed=g.product("rna_tailed","poly(G)-tailed, ddG-capped cell-indexed cDNA on beads",MaterialKind.MOLECULAR,state=rna.state)
 g.react("tail_rna","Add a 3-prime poly(G) tail to bead-bound cDNA (37 °C, 20 min)",StepKind.TAILING,((rna,InputRole.TEMPLATE),(tdt,InputRole.CATALYST)),(rna_tailed,))
 anchor=g.input("rna_anchor","printed Post-TdT-PolyC12-S primer",MaterialKind.OLIGO)
 anchor_t7=g.input("rna_anchor_t7","printed Post-TdT-P5(T7) primer",MaterialKind.OLIGO)
 anchor_polymerase=g.input("anchor_polymerase","NEBNext Ultra II Q5 master mix",MaterialKind.CATALYST)
 rna_anchored=g.product("rna_anchored","pre-amplified cDNA (13-18 cycles)",MaterialKind.MOLECULAR,state=rna_tailed.state)
 g.react("extend_rna_anchor","Post-TdT PCR from the poly(G) tail and the Round2 T7 arm",StepKind.PCR,((rna_tailed,InputRole.TEMPLATE),(anchor,InputRole.STRUCTURAL),(anchor_t7,InputRole.STRUCTURAL),(anchor_polymerase,InputRole.CATALYST)),(rna_anchored,),geometry=PCRGeometry.unavailable(RNA,reason="Supplementary Methods name Post-TdT-PolyC12-S (poly(C)12 on the dG tail) and Post-TdT-P5(T7) (on the Round2 T7 arm), but the dG tail length is not defined, so no exact pre-PCR template is drawn."))

 p5ad=g.input("nextera_p5_ad","5-prime-phosphorylated Nextera-P5-AD-F-TNA/TA with Hi-C-AD1-Block (Supplementary Methods)",MaterialKind.ADAPTER)
 p5_ligase=g.input("p5ad_ligase","T4 DNA ligase",MaterialKind.CATALYST)
 dna_capped=g.product("dna_capped","contact DNA with Nextera-P5-AD on remaining restriction ends",MaterialKind.MOLECULAR,state=dna.state)
 g.react("ligate_p5ad","Ligate non-barcoded Nextera-P5-AD to precipitated contact DNA (22 °C, ≥20 h)",StepKind.LIGATION,((dna,InputRole.TEMPLATE),(p5ad,InputRole.STRUCTURAL),(p5_ligase,InputRole.CATALYST)),(dna_capped,),geometry=AdapterLigationGeometry.unavailable("Supplementary Table 1a prints the TNA/TA forward strands and a 13-nt Nextera-P5-AD-R, but the Hi-C-AD1-Block sequence and the annealed duplex used here are not shown in the inspected tables."))
 dna_f=g.input("dna_preamp_forward","printed Nextera-P5-pre-primer",MaterialKind.OLIGO)
 dna_r=g.input("dna_preamp_reverse","printed TruSeq-P7-pre-P-S primer",MaterialKind.OLIGO)
 dna_pol=g.input("dna_preamp_polymerase","NEBNext Ultra II Q5 master mix (9 cycles)",MaterialKind.CATALYST)
 dna_preamp=g.product("dna_preamp","preamplified chromatin-contact products",MaterialKind.MOLECULAR,state=dna.state)
 g.react("preamp_dna","Preamplify the chromatin-contact branch",StepKind.PCR,((dna_capped,InputRole.TEMPLATE),(dna_f,InputRole.STRUCTURAL),(dna_r,InputRole.STRUCTURAL),(dna_pol,InputRole.CATALYST)),(dna_preamp,),geometry=DNA_PRE_GEOMETRY)

 transposome=g.input("finishing_transposome","Illumina TDE1 (Hi-C) / Nextera XT (cDNA) transposome",MaterialKind.ADAPTER)
 tn5=g.input("finishing_tn5","Tn5 transposase",MaterialKind.CATALYST)
 branch_templates=[]
 for prefix,template,final in (("dna",dna_preamp,DNA),("rna",rna_anchored,RNA)):
  tagged=g.product(prefix+"_tagged",prefix.upper()+" branch after Nextera tagmentation",MaterialKind.MOLECULAR,state=final)
  waste=g.product(prefix+"_tn5_waste","non-library "+prefix.upper()+" tagmentation products",MaterialKind.WASTE)
  g.react("tagment_"+prefix,"Tagment the "+prefix.upper()+" branch for library completion (55 °C, 5 min)",StepKind.TAGMENTATION,((template,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,waste))
  g.discard(waste,"does not retain the selected unlike-ended Nextera insert")
  reviewed.append("tagment_"+prefix);branch_templates.append(tagged)

 reviewed.extend(("proteinase_k","reverse_crosslink","streptavidin","separate","tdt","tail_rna","rna_anchor","rna_anchor_t7","anchor_polymerase","extend_rna_anchor","nextera_p5_ad","p5ad_ligase","ligate_p5ad","dna_preamp_forward","dna_preamp_reverse","dna_preamp_polymerase","preamp_dna","finishing_transposome","finishing_tn5"))
 finish_pcr_branch(g,"dna_final",branch_templates[0],DNA,DNA_P,reviewed,label="GAGE-seq chromatin-contact",read_lengths=READ_LENGTHS["Chromatin-contact library"],geometry=DNA_FINAL_GEOMETRY)
 finish_pcr_branch(g,"rna_final",branch_templates[1],RNA,RNA_P,reviewed,label="GAGE-seq RNA",read_lengths=READ_LENGTHS["RNA library"],geometry=RNA_FINAL_GEOMETRY)
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
