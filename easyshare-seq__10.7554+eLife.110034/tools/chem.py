"""easySHARE-seq dual-modality, three-part combinatorial indexing.

Authority: the published eLife article (Soltys et al., eLife 2026 version of record, eLife.110034,
PMC13384510) Methods, and its Supplementary file 1 (md5-identical to the bioRxiv
media-6 workbook; sheets "PCR oligos", "BC oligos R1", "BC oligos R2").
"""
import illumina as il
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, MolecularState, Row, Scene, Segment, Workflow, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SourceClaim,
                            SourceDiscrepancy, StepKind, Disposition, Terminal, SequencingPlan, Highlight)
from protocol_model import molecule, role

seg=Segment

TITLE="easySHARE-seq — joint RNA and chromatin accessibility"
NOTES="01_easyshare-seq.html"
SOURCE=('Defining source: <a href="https://doi.org/10.7554/eLife.110034">Soltys et al., <i>eLife</i> (2026, version of record)</a> '
        'and its Supplementary file 1; first posted as <a href="https://doi.org/10.1101/2024.02.26.581705">bioRxiv 2024.02.26.581705</a>.')
SUMMARY=("Fixed nuclei are tagmented with a mixed Tn5-A/Tn5-B2S transposome and reverse-transcribed "
         "with a biotinylated RT primer; both handles receive two 7-nt ligation barcodes. "
         "Streptavidin separates cDNA from ATAC fragments, and an N5XX PCR index completes each library.")
CAVEAT=("Representative Round1_001/Round2_001 are shown. N501(RNA) and N502(ATAC) are distinct illustrative N5XX choices, as the published joint run distinguishes modalities with different Index2 values; actual run mapping remains unassigned. Each well or sub-library "
        "substitutes another published family member. Fragment-end classes and the 72 °C gap fill "
        "are derived from the published oligos, not drawn by the source.")
R1_001="GGTTAGCGTCTCGT"
R1_001_SCAFFOLD="CCGAGCCCACGAGACGCTAACCTCT"
R2_001="CAAGCAGAAGACGGCATACGAGATTCTGGCAAGA"
R2_001_SCAFFOLD="TGCCAGAATCTCGTATGCCGTCTTCTGCTTG"
TN5_A="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
TN5_REVERSE="CTGTCTCTTATACACATCT"
TN5_B2S="GGGCTCGGAGATGTGTATAAGAGACAG"
RT_PRIMER="GGGCTCGGAGATGTGTATAAGAGACAGNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTVN"  # T/iBiodT/T24
TSO="UUUTCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGG"
KRT_B2S="CAGGGCATTTGCACCCCATGCC"   # eLife Key resources table value (see discrepancy)
I7_TRUSEQ_LONG=il.P7
N501=nx.n5xx_primer("AGCGCTAG")
N502=nx.n5xx_primer("GATATCGA")
LINKER="TCT"  # Supplementary file 1: bases 8-10 of Index 1
READ_LENGTHS={"Read 1":150,"Index 1 (i7)":17,"Index 2 (i5)":8,"Read 2":150}

def _bc1(): return seg("round-1 barcode (Index 1 bases 1–7)",revcomp(R1_001[:7]),"cbc",feature=feature("cell_bc1","cell_barcode","whitelist",group="cell_barcode",part="round 1",whitelist="Supplementary file 1 BC oligos R1 (192)"))
def _bc2(): return seg("round-2 barcode (Index 1 bases 11–17)",revcomp(R2_001[24:31]),"cbc",feature=feature("cell_bc2","cell_barcode","whitelist",group="cell_barcode",part="round 2",whitelist="Supplementary file 1 BC oligos R2 (192)"))
def _n5(index="AGCGCTAG"): return seg("N5XX index (round 3 / sub-library)",index,"cbc",feature=feature("cell_bc3","cell_barcode","whitelist",group="cell_barcode",part="round 3",whitelist="Supplementary file 1 Nextera N5XX set"))
def _umi(): return seg("10-nt UMI","N"*10,"umi",placeholder=True,feature=feature("easyshare_umi","umi","random"))
def _barcode_end(): return [seg("ME + S7 (Tn5-B2S/RT handle + BC1 tail)",nx.INDEX1_PRIMER,"r2"),_bc1(),seg("TCT linker",LINKER),_bc2(),seg("P7 reverse complement",il.P7_RC,"p7")]
def _polya(): return [seg("VN anchor (NB on this strand)","NB",placeholder=True),seg("RT-primer T26 tail (A26 on this strand)","A"*26)]

def round1_geometry():
    """Round-1 duplex meets the 5′-phosphorylated GGGCTCGG overhang of a B2S/RT end."""
    sc=Scene()
    sc.strand("round-1 oligo",[seg("Round1_001 barcode bases",R1_001[:7],"cbc",feature=_bc1().feature),seg("Round1_001 S7 tail",R1_001[7:],"r2")],label="/5Phos/Round1_001")
    sc.anneal("round-1 scaffold",[seg("scaffold B2S/RT arm",R1_001_SCAFFOLD[:8]),seg("scaffold barcode arm",R1_001_SCAFFOLD[8:-3]),seg("TCT 3′ overhang",R1_001_SCAFFOLD[-3:])],
              to="round-1 oligo",pair=("scaffold barcode arm","Round1_001 barcode bases"),label="Round1_001 scaffold")
    sc.anneal("B2S end",[seg("B2S/RT 5′ arm",TN5_B2S[:8],"r2"),seg("ME",nx.ME,"me"),seg("insert or UMI-T tail","X"*10,placeholder=True)],
              to="round-1 scaffold",pair=("B2S/RT 5′ arm","scaffold B2S/RT arm"),above=True,label="/5Phos/ Tn5-B2S or RT-primer end",unpaired=("ME","insert or UMI-T tail"))
    sc.ligation_end("round-1 scaffold","scaffold B2S/RT arm",kind="cohesive",polarity="5-prime",owner="adapter",text="8-nt scaffold arm")
    sc.ligation_end("B2S end","B2S/RT 5′ arm",kind="cohesive",polarity="5-prime",owner="substrate",text="5′-phosphorylated GGGCTCGG")
    return AdapterLigationGeometry.known("cohesive",sc,statement="Methods 'Preparing the barcoding oligonucleotides': the round-1 duplex has two single-stranded linkers, one complementary to the 5′ end of the cDNA or transposed molecule (RT primer or Tn5 adapter). Appendix 1 says this initial overhang is 8 bp and that the Read 2 primer site is split between barcode oligo and Tn5/RT primer. Supplementary file 1 prints Round1_001 and its scaffold; Scene checks that the scaffold pairs with both the barcode oligo and the GGGCTCGG arm.")

def round2_geometry():
    sc=Scene()
    sc.strand("round-2 oligo",[seg("Round2_001 P7",R2_001[:24],"p7"),seg("Round2_001 barcode bases",R2_001[24:31],"cbc",feature=_bc2().feature),seg("AGA 3′ overhang",R2_001[31:])],label="Round2_001")
    sc.anneal("round-1 scaffold",[seg("scaffold B2S/RT arm",R1_001_SCAFFOLD[:8]),seg("scaffold barcode arm",R1_001_SCAFFOLD[8:-3]),seg("TCT 3′ overhang",R1_001_SCAFFOLD[-3:])],
              to="round-2 oligo",pair=("TCT 3′ overhang","AGA 3′ overhang"),label="Round1_001 scaffold (in ligated round-1 product)")
    sc.anneal("round-1 oligo",[seg("Round1_001 barcode bases",R1_001[:7],"cbc",feature=_bc1().feature),seg("Round1_001 S7 tail",R1_001[7:],"r2")],
              to="round-1 scaffold",pair=("Round1_001 barcode bases","scaffold barcode arm"),above=True,label="/5Phos/Round1_001")
    sc.anneal("round-2 scaffold",[seg("Round2_001 scaffold",R2_001_SCAFFOLD)],to="round-2 oligo",pair=("Round2_001 scaffold","Round2_001 P7"),label="Round2_001 scaffold")
    sc.overhang("round-2 oligo","AGA 3′ overhang",polarity="3-prime",owner="adapter",text="Round2_001 3′-AGA")
    sc.ligation_end("round-1 scaffold","TCT 3′ overhang",kind="cohesive",polarity="3-prime",owner="substrate",text="round-1 scaffold 3′-TCT")
    return AdapterLigationGeometry.known("cohesive",sc,statement="Methods: the second 3-bp overhang of the round-1 duplex ligates the round-2 barcodes; round-1 oligos are 5′-phosphorylated. Supplementary file 1 prints Round2_001 and its scaffold; Scene checks the AGA/TCT pairing.")

def _library(name,payload_parts,index):
    lib=Construct([seg("P5",il.P5,"p5"),_n5(index),seg("Nextera Read 1 arm (S5 + ME)",nx.READ1_PRIMER,"r1"),*payload_parts,*_barcode_end()],name=name)
    primers=sp.standard_primers(lib)
    problems=sp.verify(lib,primers)
    if problems: raise ValueError(name+": "+"; ".join(problems))
    return lib,primers
RNA,RP=_library("easySHARE-seq RNA library",[seg("cDNA insert","X"*42,placeholder=True),*_polya(),_umi()],"AGCGCTAG")
ATAC,AP=_library("easySHARE-seq ATAC library",[seg("accessible chromatin insert","X"*42,placeholder=True)],"GATATCGA")
ATAC_TEMPLATE=Construct([seg("Tn5-A end (S5 + ME)",nx.READ1_PRIMER,"r1"),seg("accessible chromatin insert","X"*42,placeholder=True),*_barcode_end()],name="released ATAC fragment after 72 °C fill-in")
CDNA_TEMPLATE=Construct([seg("TSO-derived S5 + ME + GGG",TSO[3:].replace("U","T"),"r1"),seg("full-length cDNA","X"*42,placeholder=True),*_polya(),_umi(),*_barcode_end()],name="template-switched barcoded cDNA")
RNA_TEMPLATE=Construct([seg("Tn5-A end (S5 + ME)",nx.READ1_PRIMER,"r1"),seg("cDNA insert","X"*42,placeholder=True),*_polya(),_umi(),*_barcode_end()],name="Tn5-A-tagmented cDNA after 72 °C fill-in")
FINAL_LIBRARIES=(("RNA library",RNA,RP,"Tn5-A-fragmented cDNA: Read 2 reads the UMI and poly(T); Index 1 reads BC1–TCT–BC2.","Published run: Read 1 150, Index 1 17, Index 2 8, Read 2 150 cycles (NovaSeq 6000 S4)."),
                 ("ATAC library",ATAC,AP,"Tn5-A/Tn5-B2S fragment carrying the same three-part cellular identity.","Published run: Read 1 150, Index 1 17, Index 2 8, Read 2 150 cycles (NovaSeq 6000 S4)."))
def _o(name,seq,mod=""): return f"<b>{name}</b> — 5′-{mod}{seq}-3′"
def oligos(): return [
    _o("Tn5-A",TN5_A),
    _o("Tn5-reverse",TN5_REVERSE,"/5Phos/"),
    _o("Tn5-B2S (Supplementary file 1)",TN5_B2S,"/5Phos/"),
    "<b>Reverse-transcription primer</b> — 5′-/5Phos/GGGCTCGGAGATGTGTATAAGAGACAGNNNNNNNNNNT/iBiodT/TTTTTTTTTTTTTTTTTTTTTTTTVN-3′",
    "<b>Template-switch oligo</b> — 5′-rUrUrUTCGTCGGCAGCGTCAGATGTGTATAAGAGACArGrGrG-3′",
    _o("Round1_001 barcode oligo",R1_001,"/5Phos/")+" with scaffold 5′-"+R1_001_SCAFFOLD+"-3′",
    _o("Round2_001 barcode oligo",R2_001)+" with scaffold 5′-"+R2_001_SCAFFOLD+"-3′",
    _o("i7-Tru-Seq-long primer",I7_TRUSEQ_LONG),
    _o("Nextera N501 (representative RNA N5XX)",N501),
    _o("Nextera N502 (distinct representative ATAC N5XX)",N502),
]
def sections(): return [
 ("Tag accessible DNA and copy RNA in fixed nuclei",[Row(chunks=[("chromatin ** Tn5-A / Tn5-B2S     RNA → biotin RT primer","me",False)])],"Both handles end in a 5′-phosphorylated GGGCTCGG arm."),
 ("Ligate round-1 and round-2 barcodes",round1_geometry().rows()+[Row()]+round2_geometry().rows(),"BC1–TCT–BC2 becomes Index 1."),
]

def workflow():
    w=Workflow(MolecularState("RNA and ATAC handles",tuple(sections()[0][1])))
    w.ligate_adapter("Ligate round-1 barcodes",round1_geometry(),[Row(chunks=[("BC1 ** B2S/RT handle","cbc",False)])],note="Round-1 duplexes ligate to the handle ends.")
    w.ligate_adapter("Ligate round-2 barcodes",round2_geometry(),[Row(chunks=[("BC2 ** TCT ** BC1 ** handle","cbc",False)])],note="Round-2 duplexes ligate through the TCT/AGA overhang.")
    return w

def _exact(name,seq): return Construct([seg(name,seq)],name=name)

def population(label,*parts):
    sc=Scene()
    for i,group in enumerate(parts):sc.strand(f"population{i+1}",group,label=f"{label}: population{i+1}")
    return sc


def rt_parts():
    return [Segment("B2SRT handle",TN5_B2S),_umi(),Segment("T/iBiodT/T24","T"*26,note="secondT is internalbiotin-dT"),Segment("VN","VN",placeholder=True),role("first-strand antisense cDNA; variable length")]


def barcode_prefix(rounds):
    first=[Segment("Round1 barcode",R1_001[:7],"cbc",feature=_bc1().feature),Segment("Round1 S7 tail",R1_001[7:])]
    if rounds==1:return first
    if rounds==2:return [Segment("Round2 P7",R2_001[:24]),Segment("Round2 barcode",R2_001[24:31],"cbc",feature=_bc2().feature),Segment("Round2 AGA",R2_001[31:]),*first]
    return []


def carrier(label,rounds=0,copied=False):
    dna=[role("accessible DNA population with Tn5A/B2S endclasses; native pairing unassigned")]
    if rounds:dna=[role(f"round{rounds} barcoded accessibleDNA population; independent fragments")]
    rna=[*barcode_prefix(rounds),*rt_parts()] if copied else [role("poly(A)RNA")]
    return population(label,tuple(dna),tuple(rna))


def source_barcode_duplex(number):
    sc=Scene()
    if number==1:
        sc.strand("barcode",[Segment("R1barcode",R1_001[:7],"cbc",feature=_bc1().feature),Segment("R1tail",R1_001[7:])],label="5′phosphorylated Round1_001")
        sc.anneal("scaffold",[Segment("overhang",R1_001_SCAFFOLD[:8]),Segment("paired",R1_001_SCAFFOLD[8:-3]),Segment("TCT",R1_001_SCAFFOLD[-3:])],to="barcode",pair=("paired","R1barcode"),label="Round1 scaffold")
    else:
        sc.strand("barcode",[Segment("R2P7",R2_001[:24]),Segment("R2barcode",R2_001[24:31],"cbc",feature=_bc2().feature),Segment("AGA","AGA")],label="Round2_001")
        sc.anneal("scaffold",[Segment("R2scaffold",R2_001_SCAFFOLD)],to="barcode",pair=("R2scaffold","R2P7"),label="Round2scaffold")
    return sc


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    nuclei=g.input("nuclei","formaldehyde-fixed nuclei/cells",MaterialKind.BIOLOGICAL,state=carrier("fixedcarrier"),primary=True)
    ta=g.input("tn5_a_oligo","Tn5A",MaterialKind.OLIGO,state=_exact("Tn5A",TN5_A))
    tb=g.input("tn5_b2s_oligo","5′phosphateTn5B2S; publishedSuppfile1 value",MaterialKind.OLIGO,state=_exact("Tn5B2S",TN5_B2S))
    tr=g.input("tn5_reverse","5′phosphateTn5reverse",MaterialKind.OLIGO,state=_exact("Tn5reverse",TN5_REVERSE))
    tn5=g.input("tn5","mixedA/B2S loadedin-houseTn5",MaterialKind.CATALYST)
    tagged=g.product("tagged","tagmentedchromatin and independentRNA",MaterialKind.MOLECULAR,state=carrier("tagmentedcarrier"))
    g.react("tagment","Tagment10000nuclei37°C30min",StepKind.TAGMENTATION,((nuclei,InputRole.TEMPLATE),(ta,InputRole.STRUCTURAL),(tb,InputRole.STRUCTURAL),(tr,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,))
    primer=g.input("rt_primer","5′phosphate/biotin-dT10ntUMI RTprimer",MaterialKind.OLIGO,state=Construct(rt_parts()[:-1],name="sourceRTprimer"))
    rt=g.input("rt","MaximaHMinus withPEG; published52°C/rampRT",MaterialKind.CATALYST)
    copied=g.product("copied","first-strandcDNA beside taggedchromatin",MaterialKind.MOLECULAR,state=carrier("RTcarrier",copied=True))
    g.react("reverse_transcribe","CopyRNA with biotinylatedprimer",StepKind.REVERSE_TRANSCRIPTION,((tagged,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(copied,))
    ligase=g.input("ligase","T4DNA ligase",MaterialKind.CATALYST)
    r1=g.input("round1_barcodes","sourceRound1_001/scaffold duplex",MaterialKind.ADAPTER,state=source_barcode_duplex(1))
    round1=g.product("round1","independent cDNA/ATAC populations afterround1",MaterialKind.MOLECULAR,state=carrier("round1",1,True))
    g.react("ligate_round1","Round1 splitligation25°C30min",StepKind.LIGATION,((copied,InputRole.TEMPLATE),(r1,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(round1,),geometry=round1_geometry())
    pooling=g.input("pool_condition","poolnuclei,pellet500g3min,newligationmix,redistribute",MaterialKind.CONDITION)
    pooled=g.product("pooled","redistributed carrierpopulations",MaterialKind.MOLECULAR,state=round1.state)
    wash=g.product("redistribution_wash","physically aspirated wash/pooling supernatants",MaterialKind.WASTE)
    g.react("redistribute","Pool, wash and redistribute nuclei",StepKind.PARTITION,((round1,InputRole.TEMPLATE),(pooling,InputRole.CONDITION)),(pooled,wash),changes_molecule=False)
    g.discard(wash,"Publishedaspiration step; nuclearcarryover unquantified and unrelated to computational cell-QC exclusion")
    r2=g.input("round2_barcodes","sourceRound2_001/scaffold duplex",MaterialKind.ADAPTER,state=source_barcode_duplex(2))
    indexed=g.product("indexed","round2 barcoded cDNA/ATAC populations",MaterialKind.MOLECULAR,state=carrier("round2",2,True))
    g.react("ligate_round2","Round2 splitligation25°C45min",StepKind.LIGATION,((pooled,InputRole.TEMPLATE),(r2,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(indexed,),geometry=round2_geometry())
    reversal=g.input("reversal_condition","3500cell sublibraries;SDS/NaCl/PK62°C2h,80°C10min,Tweenquench",MaterialKind.CONDITION)
    released=g.product("released","released independent barcodedpopulations",MaterialKind.MOLECULAR,state=indexed.state)
    g.react("reverse_crosslinks","Make sublibraries and reversecrosslinks",StepKind.DEPROTEINIZATION,((indexed,InputRole.TEMPLATE),(reversal,InputRole.CONDITION)),(released,))
    bead=g.input("streptavidin","M280streptavidin1h25°C, B&W washes",MaterialKind.SELECTION)
    bound=g.product("cdna_beads","bead-retained biotinylated single-strandcDNA",MaterialKind.MOLECULAR,state=molecule("barcodedsscDNA",tuple([*barcode_prefix(2),*rt_parts()])))
    atac_sup=g.product("atac_supernatant","ATAC-enriched supernatant; native ends unassigned",MaterialKind.MOLECULAR,state=molecule("actualATACpool",(role("barcodedTn5-fragment population; gaps/ends unassigned"),)))
    g.react("streptavidin_split","Separate bead-bound cDNA and ATACsupernatant",StepKind.AFFINITY_CAPTURE,((released,InputRole.TEMPLATE),(bead,InputRole.SELECTION)),(bound,atac_sup))
    column=g.input("minelute","QiagenMinElute;30µlTris eluate",MaterialKind.SELECTION)
    atac=g.product("atac_branch","purified pre-PCR ATACpopulation",MaterialKind.MOLECULAR,state=atac_sup.state)
    awaste=g.product("atac_cleanup_waste","MinElute flowthrough/washes",MaterialKind.WASTE)
    g.react("atac_cleanup","Purify ATACsupernatant",StepKind.CLEANUP,((atac_sup,InputRole.TEMPLATE),(column,InputRole.SELECTION)),(atac,awaste),changes_molecule=False)
    g.discard(awaste,"Physicalcolumnpurification")
    tso=g.input("tso","rU×3/rG×3 sourceTSO",MaterialKind.OLIGO,state=_exact("sourceTSO",TSO))
    rt2=g.input("ts_rt","MaximaHMinus;25°C30min,42°C90min",MaterialKind.CATALYST)
    switched=g.product("switched","INFERRED template-switchedsscDNA",MaterialKind.MOLECULAR,state=molecule("switchedsscDNA",tuple([*barcode_prefix(2),*rt_parts(),role("copiedTSOend; nativejunction/extent unassigned",inferred=True)])))
    g.react("template_switch","INFERRED nativeend — on-bead templateswitch",StepKind.TEMPLATE_SWITCHING,((bound,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt2,InputRole.CATALYST)),(switched,))
    n5c=g.input("cdna_n5","representativeRNA N501; same index reused at finalRNA PCR",MaterialKind.OLIGO,state=Construct([Segment("P5",il.P5),_n5("AGCGCTAG"),Segment("S5",nx.S5)]))
    i7=g.input("i7_long","source i7-TruSeq-long",MaterialKind.OLIGO,state=_exact("i7long",I7_TRUSEQ_LONG))
    q5c=g.input("cdna_q5","Q5;13–15cycles,66°Canneal/3minextension",MaterialKind.CATALYST)
    full=Construct([Segment("P5",il.P5),_n5(),*list(CDNA_TEMPLATE)],name="amplifiedfull-lengthcDNA; sourceprimer-derivedreference")
    cdna=g.product("cdna_amplified","duplexfull-lengthcDNA",MaterialKind.MOLECULAR,state=full)
    g.react("cdna_pcr","INFERRED endpoint — sourcecDNA PCR13–15cycles",StepKind.PCR,((switched,InputRole.TEMPLATE),(n5c,InputRole.STRUCTURAL),(i7,InputRole.STRUCTURAL),(q5c,InputRole.CATALYST)),(cdna,),geometry=PCRGeometry.unavailable(switched.state,reason="Actualinput is single-strand template-switchedcDNA; full first-cycle pairing/copyextent unresolved, not an expectedduplex used as template."))
    cb=g.input("cdna_cleanup_beads","0.7×lowercutoff;25µlTris",MaterialKind.SELECTION)
    clean=g.product("cdna_clean","purifiedfull-lengthcDNA",MaterialKind.MOLECULAR,state=full)
    cwaste=g.product("cdna_cleanup_waste","0.7×supernatant/washes",MaterialKind.WASTE)
    g.react("cdna_cleanup","Purify sourcecDNA PCR",StepKind.CLEANUP,((cdna,InputRole.TEMPLATE),(cb,InputRole.SELECTION)),(clean,cwaste),changes_molecule=False);g.discard(cwaste,"Sourcebeadpurification")
    only=g.input("tn5_a_only_adapter","Tn5A/reverse duplex",MaterialKind.ADAPTER,state=(Construct([Segment("A",TN5_A)]),Construct([Segment("reverse",TN5_REVERSE)])))
    tna=g.input("tn5_a_only","Tn5A-only;25ngcDNA55°C10min,SDSstop55°C10min",MaterialKind.CATALYST)
    rna_raw=g.product("rna_tagmented","tagmentedRNA-library precursor; nativeends/gaps unassigned",MaterialKind.MOLECULAR,state=molecule("actualRNApreindexpool",(role("tagmented cDNA retaining barcodedRTend; nativegeometry unassigned"),)))
    other=g.product("rna_tn5_other","candidatecDNA fragments withoutbarcodeend",MaterialKind.MOLECULAR,inferred=True)
    g.react("tagment_cdna","Tagment cDNA with A-onlyTn5",StepKind.TAGMENTATION,((clean,InputRole.TEMPLATE),(only,InputRole.STRUCTURAL),(tna,InputRole.CATALYST)),(rna_raw,other))
    g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,"Nonamplification doesnot prove physicalremoval; recovery through cleanup unassigned."))
    rb=g.input("rna_tagment_cleanup_beads","1.3×custombeads;30µl elution",MaterialKind.SELECTION)
    rna=g.product("rna_branch","purifiedactual pre-indexcDNA pool",MaterialKind.MOLECULAR,state=rna_raw.state)
    rwaste=g.product("rna_tagment_cleanup_waste","1.3×supernatant/washes",MaterialKind.WASTE)
    g.react("rna_tagment_cleanup","Clean tagmentedcDNA",StepKind.CLEANUP,((rna_raw,InputRole.TEMPLATE),(rb,InputRole.SELECTION)),(rna,rwaste),changes_molecule=False);g.discard(rwaste,"Sourcephysicalcleanup")
    n5a=g.input("atac_n5","representativeATAC N502; different fromRNA N501",MaterialKind.OLIGO,state=Construct([Segment("P5",il.P5),_n5("GATATCGA"),Segment("S5",nx.S5)]))
    rpol=g.input("rna_q5","Q5;72°C6min,5cycles66°Canneal",MaterialKind.CATALYST)
    apol=g.input("atac_q5","Q5;72°C6min,15–17cycles66°Canneal",MaterialKind.CATALYST)
    for key,template,primer,pol,final,primers,ratios in (("rna",rna,n5c,rpol,RNA,RP,"0.5×upper/0.8×lower"),("atac",atac,n5a,apol,ATAC,AP,"0.55×upper/1.4×lower")):
        raw=g.product(key+"_pcr_product","indexed "+key+" PCRproduct",MaterialKind.MOLECULAR,state=final)
        g.react(key+"_pcr","INFERRED gap-fill geometry — source "+key+" PCR",StepKind.PCR,((template,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(i7,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(raw,),geometry=PCRGeometry.unavailable(template.state,reason="Actualnative pre-PCR tagmentedpopulation;72°Cincubation is interpreted as gapfill but fullnative first-cycle geometry is not reconstructed. Sourceprimer-derived finishedlibrary is separate."))
        selection=g.input(key+"_final_beads",ratios+" custombeads,25µlTris",MaterialKind.SELECTION)
        library=g.product(key+"_library","purified "+key+" sequencinglibrary",MaterialKind.LIBRARY,state=final)
        waste=g.product(key+"_final_waste","sizeexcluded bead/supernatant/washfractions",MaterialKind.WASTE)
        g.react(key+"_final_cleanup","Source double-sided "+key+" librarycleanup",StepKind.CLEANUP,((raw,InputRole.TEMPLATE),(selection,InputRole.SELECTION)),(library,waste),changes_molecule=False);g.discard(waste,"Sourcephysicalupper/lowercutoffs, not absence of exponentialamplification")
        g.sequence(library,SequencingPlan(library=final,primers=primers,acquired_roles=("Read 1","Index 1 (i7)","Index 2 (i5)","Read 2"),read_lengths=READ_LENGTHS))
    g.annotate(Highlight("reverse_transcribe","RNA-template disposition unassigned","InitialRT adds cDNA beside chromatin. RNA pairing/recovery afterRT, heat and washes is not quantitatively reconstructed or treated as a physicaldiscard."))
    g.annotate(SourceDiscrepancy(
        "tn5_b2s_oligo", "Tn5-B2S sequence differs between the eLife Key resources table and Supplementary file 1",
        (SourceClaim("eLife Methods, Key resources table (eLife v1 XML line 257; PMC13384510 XML lines 286–289)", KRT_B2S),
         SourceClaim("eLife Supplementary file 1, sheet 'PCR oligos' row Tn5-B2S", "/5Phos/"+TN5_B2S)),
        "The model uses the Supplementary file 1 sequence.",
        "The table value has no Tn5 mosaic end and does not match the RT primer's shared handle or the round-1 scaffold arm; the workbook value does both.",
        "The origin and resolution of the table value remain unassigned within the inspected published-source scope."))
    g.declare_reviewed_scope(("nuclei","tn5_a_oligo","tn5_b2s_oligo","tn5_reverse","tn5","tagment","rt_primer","rt","reverse_transcribe","ligase","round1_barcodes","ligate_round1","pool_condition","redistribute","round2_barcodes","ligate_round2","reversal_condition","reverse_crosslinks","streptavidin","streptavidin_split","minelute","atac_cleanup","tso","ts_rt","template_switch","cdna_n5","i7_long","cdna_q5","cdna_pcr","cdna_cleanup_beads","cdna_cleanup","tn5_a_only_adapter","tn5_a_only","tagment_cdna","rna_tagment_cleanup_beads","rna_tagment_cleanup","atac_n5","rna_q5","atac_q5","rna_pcr","atac_pcr","rna_final_beads","rna_final_cleanup","rna_library","atac_final_beads","atac_final_cleanup","atac_library"),SOURCE)
    return g.require_publishable()
