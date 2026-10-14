"""GAGE-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct,MolecularState,Row,Segment,Workflow,feature,oligo,revcomp
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

RT_DT=RT_HANDLE+"TCG"+"T"*20+"VN"
RT_RANDOM=RT_HANDLE+"TCG"+"N"*6
POST_TDT_POLYC="AGAGTGAAC"+"C"*12
POST_TDT_P5_T7="AATTTAATACGACTCACTATAGGGA"
DNA_AD1_F_TNA="TNAA"+RNA_R1_BC+nx.ME_RC+"A"
DNA_AD1_F_TA="TAA"+RNA_R1_BC+nx.ME_RC+"A"
DNA_AD1_R="TGACTTGTAGTGCAGCGTGTACTT"
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

DNA=Construct([seg("P5",il.P5,"p5"),seg("i5",DNA_I5,"cbc",feature=feature("gage_dna_i5","sample_index","whitelist",whitelist="Supplementary Table 1b Nextera-P5 series")),seg("S5",nx.S5,"s5"),seg("left ME",nx.ME,"me"),*payload("proximity-ligated genomic junction"),seg("TruSeq Read 2 arm reverse complement",revcomp(il.TRUSEQ_READ2),"r2"),seg("i7 reverse complement",revcomp(DNA_I7),"cbc",feature=feature("gage_dna_i7","sample_index","whitelist",whitelist="Supplementary Table 1b TruSeq-8bp series")),seg("P7 reverse complement",il.P7_RC,"p7")],name="GAGE-seq chromatin-contact library")
DNA_P=(sp.NEXTERA["R1"],sp.TRUSEQ["I1"],sp.NEXTERA["I2"],sp.TRUSEQ["R2"])
RNA=Construct([seg("P5",il.P5,"p5"),seg("i5",RNA_I5,"cbc",feature=feature("gage_rna_i5","sample_index","whitelist",whitelist="Supplementary Table 1b TruSeq-F-L series")),seg("TruSeq Read 1 arm",il.TRUSEQ_READ1,"r1"),seg("RT UMI","U"*10,"umi",placeholder=True,feature=UMI_F),seg("cell barcode round 2","C"*8,"cbc",placeholder=True,feature=CELL2),seg("cell barcode round 1","B"*8,"cbc",placeholder=True,feature=CELL1),seg("cDNA insert","X"*42,placeholder=True),seg("right ME reverse complement",nx.ME_RC,"me"),seg("S7 reverse complement",nx.S7_RC,"s7"),seg("i7 reverse complement",revcomp(RNA_I7),"cbc",feature=feature("gage_rna_i7","sample_index","whitelist",whitelist="Supplementary Table 1b Nextera-P7 series")),seg("P7 reverse complement",il.P7_RC,"p7")],name="GAGE-seq RNA library")
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
CAVEAT="The published Supplementary Table 1 gives every 96-well RNA and Hi-C adapter sequence plus the non-barcoded and terminal PCR primers. The page shows A1 representatives while retaining the two 96-member barcode whitelists; B/C/U/N/X denote varying barcode, UMI, random or biological sequence."
FINAL_LIBRARIES=(("Chromatin-contact library",DNA,DNA_P,"Proximity-ligated DNA is preamplified, tagmented and indexed.","Paired reads traverse the two genomic sides of captured contacts."),("RNA library",RNA,RNA_P,"Biotinylated cDNA is recovered, dC-tailed, amplified and tagmented.","Cell barcode rounds and the RT UMI identify each transcript observation."))
READ_LENGTHS={"Chromatin-contact library":{"Read 1":150,"Read 2":150},"RNA library":{"Read 1":150,"Read 2":150}}
def barcode_geometry():return AdapterLigationGeometry.known("splinted",RNA_ROUND1.adapter_scene(),RNA_ROUND2.adapter_scene(),statement="Published Supplementary Table 1c-d prints both strands of the two representative RNA splinted-ligation rounds; the shared constructor verifies both complementary arms.")

def oligos():
 for name,seq,mods in (("BiodT-dT(20)VN RT primer",RT_DT,"/5Phos/; internal biotin-dT"),("BiodT random-hexamer RT primer",RT_RANDOM,"/5Phos/; internal biotin-dT"),("RNA-ligation linker-1",RNA_LINKER1,""),("RNA-ligation linker-2",RNA_LINKER2_RNA,"RNA"),("Round1_01",RNA_R1_LEFT+RNA_R1_BC+RNA_R1_RIGHT,"/5Phos/"),("Round2_01",RNA_R2_LEFT+"N"*10+RNA_R2_BC+RNA_R2_RIGHT,""),("Post-TdT-PolyC12-S",POST_TDT_POLYC,""),("Post-TdT-P5(T7)-primer",POST_TDT_P5_T7,""),("Hi-C AD-1_F_1 (TNA)",DNA_AD1_F_TNA,"/5Phos/"),("Hi-C AD-1_R_1",DNA_AD1_R,"/5Phos/"),("Hi-C AD2 low A1",DNA_AD2_LO,""),("Hi-C AD2 high A1",DNA_AD2_HI,""),("Nextera-P5-pre-primer",DNA_PRE_F,""),("TruSeq-P7-pre-P-S",DNA_PRE_R,""),("Nextera-P5-2A",DNA_FINAL_F,""),("TruSeq-8bp-R1",DNA_FINAL_R,""),("TruSeq-F-L-1",RNA_FINAL_F,""),("Nextera-P7-6A",RNA_FINAL_R,"")):
  yield oligo(name,[Segment(name,seq,placeholder=any(x in seq for x in "N"))],mods=mods)

def sections():return [("Fix and prepare contacts",proximity_product(left="locus A",right="locus B",bridge="proximity-ligation junction").rows(),"Crosslinked chromatin is digested, filled and ligated in situ."),("Reverse-transcribe RNA",[Row(chunks=[("5′-phosphorylated biotin-poly(dT) or random hexamer → cDNA","umi",False)])],"The printed RT primers provide the exact 15-base handle used by RNA ligation round 1."),("Ligate two cell indexes",barcode_geometry().rows(),"The A1 RNA adapter pair is shown exactly; the source also prints all 96 RNA members and separate 96-member Hi-C duplex-adapter plates."),("Separate and build libraries",[Row(chunks=[("streptavidin cDNA ──► TruSeq-P5 / Nextera-P7 RNA library\nsupernatant DNA  ──► Nextera-P5 / TruSeq-P7 contact library",None,False)])],"Each branch is preamplified, tagmented and completed with its source-printed hybrid primer pair for PE150.")]
def workflow():
 s=sections();w=Workflow(MolecularState(s[0][0],tuple(s[0][1])));w.react(s[1][0],s[1][1],note=s[1][2]);w.ligate_adapter(s[2][0],barcode_geometry(),s[2][1],note=s[2][2]);return w


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input("cells","fixed cells containing chromatin contacts and RNA",MaterialKind.BIOLOGICAL,state=sections()[0][1],primary=True)
 nuclease=g.input("restriction_nuclease","chromatin restriction nuclease",MaterialKind.CATALYST)
 digested=g.product("digested","restriction-digested crosslinked chromatin with retained RNA",MaterialKind.MOLECULAR,state=cells.state)
 g.react("digest","Digest crosslinked chromatin in situ",StepKind.DIGESTION,((cells,InputRole.TEMPLATE),(nuclease,InputRole.CATALYST)),(digested,))
 fill_polymerase=g.input("fill_polymerase","restriction-end fill-in polymerase",MaterialKind.CATALYST)
 filled=g.product("filled","filled chromatin restriction ends with retained RNA",MaterialKind.MOLECULAR,state=digested.state)
 g.react("fill","Fill restriction-fragment ends",StepKind.EXTENSION,((digested,InputRole.TEMPLATE),(fill_polymerase,InputRole.CATALYST)),(filled,))
 contact_ligase=g.input("contact_ligase","proximity-ligation DNA ligase",MaterialKind.CATALYST)
 contacts=g.product("contacts","proximity-ligated chromatin contacts with retained RNA",MaterialKind.MOLECULAR,state=sections()[0][1])
 g.react("proximity_ligate","Ligate spatially adjacent chromatin ends in situ",StepKind.PROXIMITY_LIGATION,((filled,InputRole.TEMPLATE),(contact_ligase,InputRole.CATALYST)),(contacts,))

 rt_primer=g.input("rt_primer","biotinylated UMI-bearing poly(dT)/random-hexamer RT primers",MaterialKind.OLIGO)
 rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
 copied=g.product("copied","cells containing UMI-labelled biotinylated cDNA and proximity-ligated DNA",MaterialKind.MOLECULAR,state=sections()[1][1])
 g.react("reverse_transcribe","Reverse-transcribe RNA with the affinity-tagged UMI primers",StepKind.REVERSE_TRANSCRIPTION,((contacts,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(copied,))

 barcode_ligase=g.input("barcode_ligase","cell-index DNA ligase",MaterialKind.CATALYST)
 current=copied
 reviewed=["cells","restriction_nuclease","digest","fill_polymerase","fill","contact_ligase","proximity_ligate","rt_primer","rt","reverse_transcribe","barcode_ligase"]
 for number in (1,2):
  adapter=g.input(f"barcode_{number}",f"source-printed round-{number} analyte-specific barcode adapter",MaterialKind.ADAPTER)
  indexed=g.product(f"indexed_{number}",f"DNA and cDNA after cell-barcode round {number}",MaterialKind.MOLECULAR,state=sections()[2][1])
  g.react(f"ligate_{number}",f"Ligate cell-barcode round {number}",StepKind.LIGATION,((current,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(barcode_ligase,InputRole.CATALYST)),(indexed,),geometry=barcode_geometry())
  reviewed.extend((adapter.id,f"ligate_{number}"));current=indexed

 streptavidin=g.input("streptavidin","streptavidin affinity matrix",MaterialKind.SELECTION)
 rna=g.product("rna_branch","biotin-captured cDNA branch",MaterialKind.MOLECULAR,state=payload("cDNA insert",True))
 dna=g.product("dna_branch","proximity-ligated DNA flow-through branch",MaterialKind.MOLECULAR,state=payload("proximity-ligated genomic junction"))
 g.react("separate","Separate biotinylated cDNA from genomic contact DNA",StepKind.AFFINITY_CAPTURE,((current,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),(rna,dna),changes_molecule=False)

 tdt=g.input("tdt","terminal deoxynucleotidyl transferase",MaterialKind.CATALYST)
 rna_tailed=g.product("rna_tailed","poly(C)-tailed cell-indexed cDNA",MaterialKind.MOLECULAR,state=rna.state)
 g.react("tail_rna","Add a 3-prime poly(C) tail to cDNA",StepKind.TAILING,((rna,InputRole.TEMPLATE),(tdt,InputRole.CATALYST)),(rna_tailed,))
 anchor=g.input("rna_anchor","printed Post-TdT poly(C)12-S and P5(T7) primers",MaterialKind.OLIGO)
 anchor_polymerase=g.input("anchor_polymerase","anchor-extension polymerase",MaterialKind.CATALYST)
 rna_anchored=g.product("rna_anchored","anchor-extended cDNA",MaterialKind.MOLECULAR,state=rna_tailed.state)
 g.react("extend_rna_anchor","Anneal and extend the cDNA amplification anchor",StepKind.EXTENSION,((rna_tailed,InputRole.TEMPLATE),(anchor,InputRole.STRUCTURAL),(anchor_polymerase,InputRole.CATALYST)),(rna_anchored,))

 dna_f=g.input("dna_preamp_forward","printed Nextera-P5-pre-primer",MaterialKind.OLIGO)
 dna_r=g.input("dna_preamp_reverse","printed TruSeq-P7-pre-P-S primer",MaterialKind.OLIGO)
 dna_pol=g.input("dna_preamp_polymerase","contact-branch preamplification polymerase",MaterialKind.CATALYST)
 dna_preamp=g.product("dna_preamp","preamplified chromatin-contact products",MaterialKind.MOLECULAR,state=dna.state)
 g.react("preamp_dna","Preamplify the chromatin-contact branch",StepKind.PCR,((dna,InputRole.TEMPLATE),(dna_f,InputRole.STRUCTURAL),(dna_r,InputRole.STRUCTURAL),(dna_pol,InputRole.CATALYST)),(dna_preamp,),geometry=DNA_PRE_GEOMETRY)

 transposome=g.input("finishing_transposome","Nextera finishing transposome",MaterialKind.ADAPTER)
 tn5=g.input("finishing_tn5","Tn5 transposase",MaterialKind.CATALYST)
 branch_templates=[]
 for prefix,template,final in (("dna",dna_preamp,DNA),("rna",rna_anchored,RNA)):
  tagged=g.product(prefix+"_tagged",prefix.upper()+" branch after Nextera tagmentation",MaterialKind.MOLECULAR,state=final)
  waste=g.product(prefix+"_tn5_waste","non-library "+prefix.upper()+" tagmentation products",MaterialKind.WASTE)
  g.react("tagment_"+prefix,"Tagment the "+prefix.upper()+" branch for library completion",StepKind.TAGMENTATION,((template,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,waste))
  g.discard(waste,"does not retain the selected unlike-ended Nextera insert")
  reviewed.append("tagment_"+prefix);branch_templates.append(tagged)

 reviewed.extend(("streptavidin","separate","tdt","tail_rna","rna_anchor","anchor_polymerase","extend_rna_anchor","dna_preamp_forward","dna_preamp_reverse","dna_preamp_polymerase","preamp_dna","finishing_transposome","finishing_tn5"))
 finish_pcr_branch(g,"dna_final",branch_templates[0],DNA,DNA_P,reviewed,label="GAGE-seq chromatin-contact",read_lengths=READ_LENGTHS["Chromatin-contact library"],geometry=DNA_FINAL_GEOMETRY)
 finish_pcr_branch(g,"rna_final",branch_templates[1],RNA,RNA_P,reviewed,label="GAGE-seq RNA",read_lengths=READ_LENGTHS["RNA library"],geometry=RNA_FINAL_GEOMETRY)
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
