"""GAGE-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from dataclasses import replace
from protocol_model import molecule,role
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct,MolecularState,Row,Scene,Segment,Workflow,feature,oligo,revcomp
from pcr import PCRGeometry,PCRPrimer
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,StepKind,SequencingPlan,Terminal,Disposition,Highlight
seg=Segment
CELL1=feature("gage_cell_round1","cell_barcode","combinatorial",group="gage_cell",part="ligation round 1",whitelist="96-well adapter set")
CELL2=feature("gage_cell_round2","cell_barcode","combinatorial",group="gage_cell",part="ligation round 2",whitelist="96-well adapter set")
UMI_F=feature("gage_rna_umi","umi","random")
RT_HANDLE="AGACAGGACACTGAC"
RNA_LINKER1="GTCAGTGTCCTGTCTCGACTGACTGACCAG"
RNA_R1_LEFT="CAACGAGGTACGACA"; RNA_R1_BC="CTGATCGT"; RNA_R1_RIGHT="CTGGTCAGTCAGTCG"
RNA_LINKER2_RNA="TGUCGUACCUCGTUGCUCAAGCACGUGGAT"
RNA_LINKER2=RNA_LINKER2_RNA.replace("U","T")
RNA_R2_LEFT="ACGACTCACTATAGGGACGACGCTCTTCCGATCT"
RNA_R2_BC="AACGTGAT"; RNA_R2_RIGHT="ATCCACGTGCTTGAG"

# Table 1a prints /iBiodT/ as a base between the 15-nt handle and TC G: a biotin-dT residue.
RT_BIOTIN_DT="T"
RT_DT=RT_HANDLE+RT_BIOTIN_DT+"TCG"+"T"*20+"VN"
RT_RANDOM=RT_HANDLE+RT_BIOTIN_DT+"TCG"+"N"*6
POST_TDT_POLYC="AGAGTGAA"+"C"*12
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
RNA_I5="CATGCGAGAA"; RNA_I7="CTCTTAGCGG"
DNA_FINAL_F=il.P5+DNA_I5+nx.S5
DNA_FINAL_R=il.P7+DNA_I7+il.TRUSEQ_READ2[:24]
RNA_FINAL_F=il.P5+RNA_I5+il.TRUSEQ_READ1
RNA_FINAL_R=il.P7+RNA_I7+nx.S7

DNA=Construct([seg("P5",il.P5,"p5"),seg("i5",DNA_I5,"cbc",feature=feature("gage_dna_i5","sample_index","whitelist",whitelist="Supplementary Table 1b Nextera-P5 series")),seg("S5",nx.S5,"s5"),seg("left ME",nx.ME,"me"),seg("genomic contact insert ending in the TNA/TA restriction end","X"*42,placeholder=True),seg("AD-1 spacer","A"),seg("cell barcode round 1","B"*8,"cbc",placeholder=True,feature=CELL1),seg("AD-1/AD-2 joint (reverse complement)",revcomp(DNA_AD1_LINK)),seg("cell barcode round 2","C"*8,"cbc",placeholder=True,feature=CELL2),seg("TruSeq Read 2 arm reverse complement",revcomp(il.TRUSEQ_READ2),"r2"),seg("i7 reverse complement",revcomp(DNA_I7),"cbc",feature=feature("gage_dna_i7","sample_index","whitelist",whitelist="Supplementary Table 1b TruSeq-8bp series")),seg("P7 reverse complement",il.P7_RC,"p7")],name="GAGE-seq chromatin-contact library")
DNA=Construct([replace(part,inferred=True) for part in DNA],name=DNA.name+" — conditional endpoint")
DNA_P=(sp.NEXTERA["R1"],sp.TRUSEQ["R2"])
RNA=Construct([seg("P5",il.P5,"p5"),seg("i5",RNA_I5,"cbc",feature=feature("gage_rna_i5","sample_index","whitelist",whitelist="Supplementary Table 1b TruSeq-F-L series")),seg("TruSeq Read 1 arm",il.TRUSEQ_READ1,"r1"),seg("Round-2 N10 UMI","U"*10,"umi",placeholder=True,feature=UMI_F),seg("cell barcode round 2","C"*8,"cbc",placeholder=True,feature=CELL2),seg("round-2/round-1 linker arms",RNA_R2_RIGHT+RNA_R1_LEFT),seg("cell barcode round 1","B"*8,"cbc",placeholder=True,feature=CELL1),seg("round-1 arm and RT handle",RNA_R1_RIGHT+RT_HANDLE),seg("biotin-dT",RT_BIOTIN_DT),seg("RT primer spacer","TCG"),seg("cDNA from dT20VN or N6 priming","X"*42,placeholder=True),seg("right ME reverse complement",nx.ME_RC,"me"),seg("S7 reverse complement",nx.S7_RC,"s7"),seg("i7 reverse complement",revcomp(RNA_I7),"cbc",feature=feature("gage_rna_i7","sample_index","whitelist",whitelist="Supplementary Table 1b Nextera-P7 series")),seg("P7 reverse complement",il.P7_RC,"p7")],name="GAGE-seq RNA library")
RNA_P=(sp.TRUSEQ["R1"],sp.NEXTERA["R2"])
for _lib,_primers in ((DNA,DNA_P),(RNA,RNA_P)):
 problems=sp.verify(_lib,_primers,required_roles=tuple(p.role for p in _primers))
 if problems: raise ValueError(_lib.name+": "+"; ".join(problems))

TITLE="GAGE-seq — gene expression and genome architecture"
NOTES="01_gage-seq.html"
SOURCE='<a href="https://doi.org/10.1038/s41588-024-01745-3">Zhou et al., <i>Nature Genetics</i> (2024)</a> (protocol first posted as bioRxiv 10.1101/2023.07.20.549578).'
SUMMARY="In fixed cells, RNA is reverse-transcribed and genomic contacts are proximity-ligated. Two rounds of well-specific adapter ligation barcode both analytes before cDNA/DNA separation and library-specific tagmentation."
CAVEAT=("The published Supplementary Table 1 gives every 96-well RNA and Hi-C adapter sequence plus the non-barcoded and terminal PCR primers. Index-primer sequences specify physical index lengths; index acquisition and acquired cycle counts remain unassigned. Native first-cycle/PCR and tagmentation geometry is bounded separately from these conditional endpoints. The page shows A1 representatives while retaining the two 96-member barcode whitelists; B/C/U/N/X denote varying barcode, UMI, random or biological sequence. "
 "Model corrections in this pass (2026-10-17): the published Methods order is RT before MseI/CviQI digestion, proximity ligation without fill-in, and DdeI re-digestion; the TdT tail is poly(G) (dGTP/ddGTP), annealed by the printed poly(C)12 primer; the Hi-C AD-1 A1 barcode and 3-prime GA, the /iBiodT/ residue of the RT primers and the Hi-C read-2 barcode order were mis-transcribed earlier. "
 "Open: the Methods trim 12 nt from DNA read 1, which matches an AD-1 BC1+A+TNA end rather than a Tn5 end; the model draws the Tn5 end and does not resolve which P5 end read 1 starts from.")
FINAL_LIBRARIES=(("Chromatin-contact library",DNA,DNA_P,"INFERRED endpoint — proximity-ligated DNA is capped with Nextera-P5-AD, preamplified, tagmented and indexed.","Read 2 reads BC2, the 15-nt AD-1/AD-2 joint, BC1 and the restriction end (Methods trim 35 nt); the drawn read-1 start is conditional because the source 12-nt trim remains unresolved."),("RNA library",RNA,RNA_P,"Biotinylated cDNA is recovered on streptavidin beads, poly(G)-tailed, amplified and tagmented.","Read 1 reads the Round-2 N10 UMI, BC2, the 30-nt linker, BC1 and the 15-nt round-1 arm (Methods positions 1-71)."))
READ_LENGTHS={"Chromatin-contact library":{"Read 1":150,"Read 2":150},"RNA library":{"Read 1":150,"Read 2":150}}


# Published PE150 is positive; omitting index primer geometry is not a claim
# that index reads were not acquired. Keep this distinction in rendered tables.
SHOW_READ_CONTENT=False
READ_SPECS={
 "Chromatin-contact library":(
  sp.read_spec("Read 1",150,"Genomic contact insert; drawn Tn5 start conditional against published12nt trim","include"),
  sp.read_spec("Read 2",150,"BC2, AD-1/AD-2 joint, BC1 and restriction end; source trims35nt before genomic alignment","extract_then_include"),
  sp.read_spec("Index 1 (i7)","unassigned","Physical i7 sample index in printed indexed reverse primer; acquisition and run cycles unassigned","undetermined"),
  sp.read_spec("Index 2 (i5)","unassigned","Physical i5 sample index in printed indexed forward primer; acquisition and run cycles unassigned","undetermined")),
 "RNA library":(
  sp.read_spec("Read 1",150,"Round2 N10 UMI, BC2, linker and BC1; source identifier extraction precedes cDNA alignment","extract_then_include"),
  sp.read_spec("Read 2",150,"RNA-derived cDNA insert","include"),
  sp.read_spec("Index 1 (i7)","unassigned","Physical i7 sample index in printed indexed reverse primer; acquisition and run cycles unassigned","undetermined"),
  sp.read_spec("Index 2 (i5)","unassigned","Physical i5 sample index in printed indexed forward primer; acquisition and run cycles unassigned","undetermined"))}
READ_SPECS_NOTE="Published paired-end150 acquisition is retained. Printed indexed primers establish physical sample-index identities, not whether index reads were acquired or their cycle counts. R1/R2 primer geometry is drawn separately; the DNA Read1 start remains conditional."

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
 rounds={1:(rna_round_scene(1),hic_ad1_scene()),2:(rna_round_scene(2),hic_ad2_scene())}
 scenes=rounds[number] if number else (*rounds[1],*rounds[2])
 return AdapterLigationGeometry.known("splinted",*scenes,statement="Published Supplementary Tables 1c-1i print both strands of the A1 RNA splinted adapters and Hi-C AD-1 (5-prime TNA/TA cohesive end, Y-shaped) and AD-2 (CAAGTCA overhang matching the AD-1_R TGACTTG arm) duplexes; local source-owned declarations verify the RNA splint arms and the Scene verifies duplex pairing. The Round2 panel includes the N10 UMI between its TruSeq arm and BC2.")

def oligos():
 for name,seq,mods in (("BiodT-dT(20)VN RT primer",RT_DT,"/5Phos/; base 16 is /iBiodT/"),("BiodT random-hexamer RT primer",RT_RANDOM,"/5Phos/; base 16 is /iBiodT/"),("RNA-ligation linker-1",RNA_LINKER1,""),("RNA-ligation linker-2",RNA_LINKER2_RNA,"RNA"),("RNA Ligation-1-blocker","CGAAGTCAGTGTCCTGTC",""),("Round1_01",RNA_R1_LEFT+RNA_R1_BC+RNA_R1_RIGHT,"/5Phos/"),("Round2_01",RNA_R2_LEFT+"N"*10+RNA_R2_BC+RNA_R2_RIGHT,""),("Post-TdT-PolyC12-S",POST_TDT_POLYC,""),("Post-TdT-P5(T7)-primer",POST_TDT_P5_T7,""),("Hi-C AD-1_F_1 (TNA)",DNA_AD1_F_TNA,"/5Phos/"),("Hi-C AD-1_R_1",DNA_AD1_R,"/5Phos/"),("Hi-C AD2 low A1",DNA_AD2_LO,""),("Hi-C AD2 high A1",DNA_AD2_HI,""),("Nextera-P5-AD-F-TNA",NEXTERA_P5_AD_F_TNA,"/5Phos/"),("Nextera-P5-AD-F-TA",NEXTERA_P5_AD_F_TA,"/5Phos/"),("Nextera-P5-AD-R",NEXTERA_P5_AD_R,""),("Nextera-P5-pre-primer",DNA_PRE_F,""),("TruSeq-P7-pre-P-S",DNA_PRE_R,""),("Nextera-P5-2A",DNA_FINAL_F,""),("TruSeq-8bp-R1",DNA_FINAL_R,""),("TruSeq-F-L-1",RNA_FINAL_F,""),("Nextera-P7-6A",RNA_FINAL_R,"")):
  yield oligo(name,[Segment(name,seq,placeholder=any(x in seq for x in "N"))],mods=mods)

def rna_round_scene(number):
    previous=RT_HANDLE if number==1 else RNA_R1_LEFT
    left=RNA_R1_LEFT if number==1 else RNA_R2_LEFT
    right=RNA_R1_RIGHT if number==1 else RNA_R2_RIGHT
    bc=RNA_R1_BC if number==1 else RNA_R2_BC
    linker=RNA_LINKER1 if number==1 else RNA_LINKER2
    if linker != revcomp(previous)+revcomp(right):
        raise ValueError("Published RNA splint does not match the two supplied arms")
    sc=Scene()
    sc.strand("previous",[Segment("previous 5′ handle",previous),role("preceding cDNA molecule")],label="preceding cDNA",mod5="/5Phos/")
    sc.anneal("splint",[Segment("previous arm",revcomp(previous)),Segment("new arm",revcomp(right))],to="previous",pair=("previous arm","previous 5′ handle"),label="linker-1 DNA" if number==1 else "linker-2 RNA (U shown as T)")
    parts=[Segment("new 5′ arm",left)]
    if number==2:parts.append(Segment("Round2 N10 UMI","N"*10,"umi",placeholder=True,feature=UMI_F))
    parts += [Segment("barcode",bc,"cbc",feature=CELL1 if number==1 else CELL2),Segment("new 3′ arm",right)]
    sc.anneal("new oligo",parts,to="splint",pair=("new 3′ arm","new arm"),above=True,label="Round1_01" if number==1 else "Round2_01",mod5="/5Phos/" if number==1 else "")
    sc.ligation_end("new oligo","new 3′ arm",kind="splinted",polarity="3-prime",owner="adapter",text="3′OH aligned by published linker")
    sc.ligation_end("previous","previous 5′ handle",kind="splinted",polarity="5-prime",owner="substrate",text="5′phosphate on preceding cDNA/round1 oligo")
    return sc


def rna_parts(rounds=0,tailed=False):
    parts=[]
    if rounds==2:parts += [Segment("Round2 5′ arm",RNA_R2_LEFT),Segment("Round2 N10 UMI","N"*10,"umi",placeholder=True,feature=UMI_F),Segment("Round2 barcode",RNA_R2_BC,"cbc",feature=CELL2),Segment("Round2 3′ arm",RNA_R2_RIGHT)]
    if rounds:parts += [Segment("Round1 5′ arm",RNA_R1_LEFT),Segment("Round1 barcode",RNA_R1_BC,"cbc",feature=CELL1),Segment("Round1 3′ arm",RNA_R1_RIGHT)]
    parts += [Segment("RT handle",RT_HANDLE),Segment("biotin-dT",RT_BIOTIN_DT,note="/iBiodT/"),Segment("RT spacer","TCG"),role("dT20VN or N6 priming residues"),role("single-strand reverse-transcribed RNA")]
    if tailed:parts += [role("poly(G) tail; variable extent"),Segment("terminal ddG","G",note="ddGTP terminator")]
    return parts


def dna_scene(label):
    return molecule(label,(role(label+"; native paired ends and complete strand assembly unassigned"),))


def mixed_scene(label,rounds=0,copied=True):
    sc=Scene()
    sc.strand("DNA population",[role(label+" chromatin DNA; restriction/contact/adapter extent source-scoped")],label="independent chromatin-DNA population")
    sc.strand("RNA population",rna_parts(rounds) if copied else [role("cellular RNA")],label="independent cDNA population" if copied else "independent RNA population")
    return sc


def protocol_graph():
    scope="Published2024 Methods and Supplementary Methods/Table1; not superseded bioRxiv chemistry"
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cells=g.input("cells","formaldehyde-crosslinked cells, lysed and SDS-treated",MaterialKind.BIOLOGICAL,state=mixed_scene("starting",copied=False),primary=True)
    primers=Scene()
    primers.strand("dT primer",[Segment("RT dT",RT_DT,placeholder=True)],label="/5Phos/; base16 /iBiodT/")
    primers.strand("random primer",[Segment("RT random",RT_RANDOM,placeholder=True)],label="independent /5Phos/ random-hexamer primer; base16 /iBiodT/")
    rp=g.input("rt_primer","Table1a biotinylated dT20VN and random-hexamer primer species",MaterialKind.OLIGO,state=primers)
    rt=g.input("rt","Maxima H Minus reverse transcriptase",MaterialKind.CATALYST)
    copied=g.product("copied","first-strand cDNA alongside chromatin DNA",MaterialKind.MOLECULAR,state=mixed_scene("RT"))
    g.react("reverse_transcribe","Reverse-transcribe RNA before chromatin digestion",StepKind.REVERSE_TRANSCRIPTION,((cells,InputRole.TEMPLATE),(rp,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(copied,))
    restriction=g.input("restriction","MseI and CviQI;25°C2h;5′TA ends",MaterialKind.CATALYST)
    cut=g.product("cut","restriction-cut chromatin plus independent cDNA",MaterialKind.MOLECULAR,state=mixed_scene("MseI/CviQI-cut"))
    g.react("digest","Digest chromatin with MseI/CviQI",StepKind.DIGESTION,((copied,InputRole.TEMPLATE),(restriction,InputRole.CATALYST)),(cut,))
    ligase=g.input("ligase","T4 DNA ligase",MaterialKind.CATALYST)
    contact_state=mixed_scene("TA-end proximity-ligated")
    # Contact loci on one genomic molecule; cDNA remains a separate molecular row.
    contact_state=Scene()
    contact_state.strand("contact DNA",[role("genomic locus A"),role("TA-overhang ligation junction; no fill-in"),role("genomic locus B")],label="one proximity-ligated genomic molecule")
    contact_state.strand("cDNA",rna_parts(),label="independent single-strand cDNA")
    contact=g.product("contact","proximity-ligated chromatin plus cDNA",MaterialKind.MOLECULAR,state=contact_state)
    g.react("proximity_ligate","Ligate TA ends16°C overnight, without fill-in",StepKind.PROXIMITY_LIGATION,((cut,InputRole.TEMPLATE),(ligase,InputRole.CATALYST)),(contact,))
    dde=g.input("dde","DdeI;37°C1h;TNA/TA5′ ends",MaterialKind.CATALYST)
    current=g.product("redigested","DdeI-redigested DNA plus independent cDNA",MaterialKind.MOLECULAR,state=mixed_scene("DdeI-cut"))
    g.react("redigest","Re-digest proximity products with DdeI",StepKind.DIGESTION,((contact,InputRole.TEMPLATE),(dde,InputRole.CATALYST)),(current,))
    for number in (1,2):
        adapter=g.input(f"barcode_{number}","published RNA and Hi-C A1 adapter species; full96-well families source-scoped",MaterialKind.ADAPTER,state=tuple((rna_round_scene(number),hic_ad1_scene() if number==1 else hic_ad2_scene())))
        if number==2:
            lysis=g.input("well_lysis_condition","60°C24h in-well SDS lysis before round2; source blocker/NexteraP5AD additions retained as recipe boundary",MaterialKind.CONDITION)
            lysed=g.product("well_lysed","lysed mixed DNA/cDNA populations",MaterialKind.MOLECULAR,state=current.state)
            g.react("well_lysis","Lyse cells in wells before second ligation",StepKind.THERMAL_TREATMENT,((current,InputRole.TEMPLATE),(lysis,InputRole.CONDITION)),(lysed,));current=lysed
        product=g.product(f"indexed_{number}",f"DNA/cDNA populations after round{number}",MaterialKind.MOLECULAR,state=mixed_scene("barcoded",number))
        g.react(f"ligate_{number}",f"Ligate source barcode round{number}",StepKind.LIGATION,((current,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(product,),geometry=barcode_geometry(number));current=product
    pk=g.input("proteinase_k","proteinaseK reverse-crosslink60°C20h; pooled isopropanol precipitation",MaterialKind.CONDITION)
    released=g.product("released","deproteinized precipitated DNA/cDNA populations",MaterialKind.MOLECULAR,state=current.state)
    g.react("release","Reverse crosslinks and recover pooled nucleic acids",StepKind.DEPROTEINIZATION,((current,InputRole.TEMPLATE),(pk,InputRole.CONDITION)),(released,))
    beads=g.input("streptavidin","MyOne C1 streptavidin beads",MaterialKind.SELECTION)
    rna=g.product("rna_branch","biotinylated single-strand cDNA on beads",MaterialKind.MOLECULAR,state=molecule("barcoded ss cDNA",tuple(rna_parts(2))))
    dna=g.product("dna_branch","chromatin-contact DNA supernatant",MaterialKind.MOLECULAR,state=dna_scene("source DNA branch"))
    g.react("separate","Capture biotinylated cDNA; retain genomic supernatant",StepKind.AFFINITY_CAPTURE,((released,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(rna,dna),changes_molecule=False)
    tdt=g.input("tdt","TdT/CoCl2;0.95µMdGTP+0.05µMddGTP;37°C20min",MaterialKind.CATALYST)
    tailed=g.product("rna_tailed","single-strand poly(G)-tailed/ddG-terminated cDNA",MaterialKind.MOLECULAR,state=molecule("actual tailed ss cDNA",tuple(rna_parts(2,True))))
    g.react("tail_rna","Tail cDNA with G/ddG",StepKind.TAILING,((rna,InputRole.TEMPLATE),(tdt,InputRole.CATALYST)),(tailed,))
    tailwash=g.input("tail_wash","remove TdT supernatant; wash beads with EB",MaterialKind.SELECTION)
    washed=g.product("rna_tail_clean","washed tailed cDNA on beads",MaterialKind.MOLECULAR,state=tailed.state)
    waste=g.product("tail_waste","TdT supernatant and bead washes",MaterialKind.WASTE)
    g.react("wash_tail","Remove TdT solution and wash beads",StepKind.CLEANUP,((tailed,InputRole.TEMPLATE),(tailwash,InputRole.SELECTION)),(washed,waste),changes_molecule=False);g.discard(waste,"source physically removes TdT supernatant and washes")
    af=g.input("rna_pre_f","Post-TdT-P5(T7) primer",MaterialKind.OLIGO,state=Construct([Segment("P5T7",POST_TDT_P5_T7)]))
    ar=g.input("rna_pre_r","Post-TdT-PolyC12-S primer",MaterialKind.OLIGO,state=Construct([Segment("polyC12",POST_TDT_POLYC)]))
    q5=g.input("q5","NEBNext UltraII Q5 master mix",MaterialKind.CATALYST)
    rna_pre=g.product("rna_preamp","post-TdT amplified cDNA; full native assembly unassigned",MaterialKind.MOLECULAR,state=dna_scene("duplex amplified RNA-derived cDNA"))
    g.react("amplify_rna_pre","Amplify post-TdT cDNA13–18cycles",StepKind.PCR,((washed,InputRole.TEMPLATE),(af,InputRole.STRUCTURAL),(ar,InputRole.STRUCTURAL),(q5,InputRole.CATALYST)),(rna_pre,),geometry=PCRGeometry.unavailable(washed.state,reason="Actual precursor is single-strand cDNA with variable G tail. The final sequencing library is not this PCR template."))
    p5=g.input("p5_adapter","Nextera-P5-AD F-TNA/TA+R and HiC-AD1-Block",MaterialKind.ADAPTER)
    capped=g.product("dna_capped","P5AD-capped contact DNA; native duplex assembly unassigned",MaterialKind.MOLECULAR,state=dna_scene("actual capped DNA"))
    g.react("p5_ligate","Ligate P5AD22°C≥20h after DNA precipitation",StepKind.LIGATION,((dna,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(capped,),geometry=AdapterLigationGeometry.unavailable("Source prints F-TNA/TA and R components; full blocker-containing duplex assembly remains a model boundary."))
    capped=cleanup(g,"dna_ligation",capped,"80µlAMPure")
    df=g.input("dna_pre_f","Nextera-P5-pre-primer",MaterialKind.OLIGO,state=Construct([Segment("DNApreF",DNA_PRE_F)]))
    dr=g.input("dna_pre_r","TruSeq-P7-pre-P-S",MaterialKind.OLIGO,state=Construct([Segment("DNApreR",DNA_PRE_R)]))
    dna_pre=g.product("dna_preamp","preamplified contact DNA",MaterialKind.MOLECULAR,state=dna_scene("pre-amplified contact DNA"))
    g.react("amplify_dna_pre","Pre-amplify contact DNA9cycles98/60/65°C",StepKind.PCR,((capped,InputRole.TEMPLATE),(df,InputRole.STRUCTURAL),(dr,InputRole.STRUCTURAL),(q5,InputRole.CATALYST)),(dna_pre,),geometry=PCRGeometry.unavailable(capped.state,reason="Native capped-DNA assembly is unassigned. Published primer sequences are explicit inputs, not proof that the final library is the precursor."))
    dna_pre=cleanup(g,"dna_pre",dna_pre,"0.75×AMPure;42µl elution")
    rna_pre=cleanup(g,"rna_pre",rna_pre,"0.8×AMPure;32µl elution")
    for prefix,template,final,primers,fs,rs,recipe in (("dna",dna_pre,DNA,DNA_P,DNA_FINAL_F,DNA_FINAL_R,"40ng;TDE1;7PCRcycles"),("rna",rna_pre,RNA,RNA_P,RNA_FINAL_F,RNA_FINAL_R,"1.5ng;NexteraXT;12PCRcycles")):
        enzyme=g.input(prefix+"_tn5",recipe+";55°C5min;NT stop",MaterialKind.CATALYST)
        adapter=g.input(prefix+"_tn5_adapter","source kit-loaded Tn5 adapter population; native product ends unassigned",MaterialKind.ADAPTER)
        tagged=g.product(prefix+"_tagged","actual pre-index tagmented "+prefix+" population",MaterialKind.MOLECULAR,state=molecule("pre-index "+prefix,(role("tagmented "+prefix+" molecule; native strand/gap assembly unassigned"),)))
        other=g.product(prefix+"_other","other tagmentation end classes",MaterialKind.MOLECULAR,inferred=True)
        g.react("tagment_"+prefix,"Tagment "+prefix+" branch55°C5min; NT stop",StepKind.TAGMENTATION,((template,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(enzyme,InputRole.CATALYST)),(tagged,other))
        g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,"Amplification compatibility does not establish physical removal; recovery through source cleanup remains unassigned."))
        f=g.input(prefix+"_final_f","source indexed forward primer",MaterialKind.OLIGO,state=Construct([Segment("F",fs)]))
        r=g.input(prefix+"_final_r","source indexed reverse primer",MaterialKind.OLIGO,state=Construct([Segment("R",rs)]))
        npm=g.input(prefix+"_npm","NPM;72°C5min then "+("7" if prefix=="dna" else "12")+" cycles",MaterialKind.CATALYST)
        amplified=g.product(prefix+"_amplified","conditional source-primer-derived "+prefix+" endpoint",MaterialKind.MOLECULAR,state=final)
        g.react(prefix+"_index_pcr","INFERRED native assembly — index "+prefix+" branch",StepKind.PCR,((tagged,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(npm,InputRole.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(tagged.state,reason="Actual input is pre-index tagmented material. Gap filling and full native product assembly are not validated by the conditional final endpoint."))
        library=cleanup(g,prefix+"_final",amplified,"0.7×AMPure;40µl elution",library=True)
        g.sequence(library,SequencingPlan(library=final,primers=primers,read_lengths={"Read 1":150,"Read 2":150},limitation="Published PE150. Index acquisition/cycle counts unassigned; DNA read1 Tn5-start representation remains conditional against source12nt trim."))
    g.annotate(Highlight("reverse_transcribe","RNA-template fate unassigned","cDNA is tracked separately from chromatin; RNA-template persistence/degradation is not declared a physical discard."))
    g.declare_reviewed_scope(REVIEWED_TARGETS,scope)
    return g.require_publishable()


def cleanup(g,prefix,template,recipe,library=False):
    selection=g.input(prefix+"_beads",recipe,MaterialKind.SELECTION)
    product=g.product(prefix+"_clean","purified "+template.label,MaterialKind.LIBRARY if library else MaterialKind.MOLECULAR,state=template.state)
    waste=g.product(prefix+"_waste","source bead-cleanup supernatant/washes",MaterialKind.WASTE)
    g.react(prefix+"_cleanup","Purify with "+recipe,StepKind.CLEANUP,((template,InputRole.TEMPLATE),(selection,InputRole.SELECTION)),(product,waste),changes_molecule=False)
    g.discard(waste,"Source physical bead cleanup; exact recovery of other end classes is not inferred.")
    return product


REVIEWED_TARGETS=(
 "cells","rt_primer","rt","reverse_transcribe","restriction","digest","ligase","proximity_ligate","dde","redigest",
 "barcode_1","ligate_1","barcode_2","well_lysis_condition","well_lysis","ligate_2","proteinase_k","release","streptavidin","separate",
 "tdt","tail_rna","tail_wash","wash_tail","rna_pre_f","rna_pre_r","q5","amplify_rna_pre","p5_adapter","p5_ligate",
 "dna_ligation_beads","dna_ligation_cleanup","dna_pre_f","dna_pre_r","amplify_dna_pre","dna_pre_beads","dna_pre_cleanup","rna_pre_beads","rna_pre_cleanup",
 "dna_tn5","dna_tn5_adapter","tagment_dna","dna_final_f","dna_final_r","dna_npm","dna_index_pcr","dna_final_beads","dna_final_cleanup","dna_final_clean",
 "rna_tn5","rna_tn5_adapter","tagment_rna","rna_final_f","rna_final_r","rna_npm","rna_index_pcr","rna_final_beads","rna_final_cleanup","rna_final_clean")
