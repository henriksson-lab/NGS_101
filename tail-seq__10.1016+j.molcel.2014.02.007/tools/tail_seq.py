"""Original TAIL-seq molecular model (Chang et al. 2014, Molecular Cell)."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct, Scene, Segment, feature, oligo, revcomp
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule,role

TITLE="TAIL-seq — poly(A)-tail length and terminal-base sequencing"
NOTES="01_tail-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.molcel.2014.02.007">Chang et al., <i>Molecular Cell</i> (2014)</a>.'
SUMMARY="A biotinylated 3-prime adapter preserves native RNA termini through partial RNase T1 digestion and capture; paired reads identify the transcript and traverse its complete poly(A) tail."
SI="Chang et al. 2014 Supplemental Information (mmc1), Extended Experimental Procedures: 'Construction of TAIL-seq library' and 'Oligonucleotide sequences for TAIL-seq'"
ARTICLE="Chang et al. 2014, Molecular Cell 53:1044-1052, Results (TAIL-seq features 1-8) and Experimental Procedures"
CAVEAT=("Reaction order, reagents, every adapter and primer base and the 51 &times; 251 HiSeq 2500 run come from the article and its Supplemental Experimental Procedures (mmc1). "
 "The 3′ adapter is a 5′-adenylated DNA oligo with 15 degenerate bases after the delimiter CTGAC (read as GTCAG in read 2), two internal biotin-dT and a 3′ dideoxy-C; the 5′ adapter is RNA. "
 "The ligases and the RT and PCR conditions are not stated, and no sequencing primer is printed: the Read 1 and Read 2 primer sites drawn below are derived from the 5′-adapter and RT-primer sequences. "
 "The reverse PCR primer's base opposite the terminal C of the adapter RTP-binding arm (before the modified TTC suffix) is T, so the amplified library carries A at that position. Native poly(A) and terminal additions have variable lengths; displayed placeholder widths are not read-cycle counts. The source does not say whether the 6-nt index was read.")

# mmc1 'Oligonucleotide sequences for TAIL-seq' (text lines 1349-1430), verbatim bases.
# Reverse primer 2 represents the set: primer 1 wraps across the page-14 number in the
# extracted text (verified by eye), so the source guard cannot match it.
ADAPTER_5P_RNA="GUUCAGAGUUCUACAGUCCGACGAUC"
DELIMITER_TOP="CTGAC"                      # read as 'GTCAG' in read-2 direction (mmc1 Sequence processing)
DEGENERATE_LEN=15
RTP_SITE="TGGAATTCTCGGGTGCCAAGGC"         # 3′-adapter 3′ arm, before /iBiodT//iBiodT//3ddC/
ADAPTER_3P_CORE=DELIMITER_TOP+"N"*DEGENERATE_LEN+RTP_SITE
# Both /iBiodT/ tokens and /3ddC/ are nucleotides, not zero-length decorations.
ADAPTER_3P=ADAPTER_3P_CORE+"TTC"

def three_prime_adapter_segments():
    return [Segment("3′-adapter delimiter",DELIMITER_TOP,"r2"),
            Segment("3′-adapter degenerate bases","N"*DEGENERATE_LEN,"umi",placeholder=True,feature=DEGENERATE),
            Segment("RTP site",RTP_SITE,"r2"),*modified_three_prime_tail()]

def modified_three_prime_tail():
    return [Segment("first internal biotin-dT","T",note="/iBiodT/"),
            Segment("second internal biotin-dT","T",note="/iBiodT/"),
            Segment("blocked terminal dideoxy-C","C",note="/3ddC/: no 3′-OH")]
RT_PRIMER="GCCTTGGCACCCGAGAATTCCA"
RP1_TAIL="GTTCAGAGTTCTACAGTCCGA"
PCR_FORWARD=il.P5+RP1_TAIL
RPI_TAIL="GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
RPI_INDEXES=("CGTGAT","ACATCG","GCCTAA","TGGTCA","CACTGT","ATTGGC","GATCTG","TCAAGT","TACAAG")
def pcr_reverse(n:int=1)->str:
    """mmc1 'PCR primer (reverse) n', n = 1-9."""
    return il.P7+RPI_INDEXES[n-1]+RPI_TAIL

if ADAPTER_5P_RNA.replace("U","T")!=il.TRUSEQ_SMALL_RNA_RA5:
    raise ValueError("5′ adapter no longer matches the RA5 constant used for the drawn arm")
if revcomp(RT_PRIMER)!=RTP_SITE:
    raise ValueError("RT primer is not the complement of the 3′-adapter arm")
if not ADAPTER_5P_RNA.replace("U","T").startswith(RP1_TAIL):
    raise ValueError("forward PCR primer does not end on the 5′-adapter sequence")
if not RT_PRIMER.endswith(RPI_TAIL[-21:]):
    raise ValueError("reverse PCR primer does not overlap the RT-primer site")

DEGENERATE=feature("tail_seq_3p_degenerate","umi","random",note="3′-adapter degenerate bases (deduplication with read-1 cycles 21-35; mmc1 Sequence processing)")
INDEX=feature("tail_seq_rpi_index","sample_index","fixed",note="6-nt index of reverse PCR primers 1-9 (mmc1); acquisition not stated")

# PCR installs RP1 on the 5′-adapter side and RPI on the 3′-adapter side, so the library is
# P5 | RA5 (DNA copy) | insert | CTGAC | N15 | revcomp(RPI tail) | revcomp(index) | P7′.
RPI_ARM=revcomp(RPI_TAIL)
FINAL_LIBRARY=Construct([
    Segment("P5","" + il.P5,"p5"),
    Segment("5′-adapter copy",il.TRUSEQ_SMALL_RNA_RA5,"r1"),
    Segment("gene-proximal RNA fragment","X"*34,placeholder=True),
    Segment("native poly(A) tail (variable length)","A"*24,placeholder=True),
    Segment("optional native terminal additions (variable length)","X"*6,placeholder=True),
    Segment("3′-adapter delimiter",DELIMITER_TOP,"r2"),
    Segment("3′-adapter degenerate bases","N"*DEGENERATE_LEN,"umi",placeholder=True,feature=DEGENERATE),
    Segment("RTP site / reverse-primer arm",RPI_ARM,"r2"),
    Segment("index (reverse complement)","I"*6,"cbc",placeholder=True,feature=INDEX),
    Segment("P7 reverse complement",il.P7_RC,"p7"),
],name="TAIL-seq library")
SEQ_PRIMERS=(
    sp.custom("Read 1","Read 1 site (5′-adapter sequence)",il.TRUSEQ_SMALL_RNA_RA5,
              SI,"Derived: no sequencing primer is printed; read 1 must start at the insert (51-nt read 1 maps the transcript)."),
    sp.custom("Read 2","Read 2 site (RT-primer sequence)",RT_PRIMER,
              SI,"Derived: no sequencing primer is printed; read-2 cycles 1-15 are the degenerate bases (mmc1). The primer's 5′ G faces the PCR-installed A."),
)
READ_LENGTHS={"Read 1":51,"Read 2":251}
FINAL_CAPTION="Read 1 identifies the gene-proximal fragment; read 2 reads the 15 degenerate bases, GTCAG, then terminal additions and the poly(A) tract from the native 3-prime end."
SEQUENCING_INTRO=("HiSeq 2500, 51 &times; 251 paired-end with 20% PhiX and 1% spike-ins (mmc1). After the 15 degenerate bases and the "
                  "5-nt delimiter, 231 read-2 cycles cover the insert from its 3′ end (article: '231 nt from the 3′ end'). "
                  "Tail length is called from raw fluorescence because base calling is unreliable through long T stretches.")
SHOW_READ_CONTENT=True


def oligos():
    return [
        f"<b>3′ adapter</b> — <code>5′-/5App/{ADAPTER_3P_CORE}/iBiodT//iBiodT//3ddC/-3′</code>",
        oligo("5′ adapter (RNA)",[Segment("5′ adapter",ADAPTER_5P_RNA,"r1")]),
        oligo("RT primer",[Segment("RTP",RT_PRIMER,"r2")]),
        oligo("PCR primer (forward)",[Segment("P5",il.P5,"p5"),Segment("5′-adapter arm",RP1_TAIL,"r1")]),
        oligo("PCR primer (reverse) 2",[Segment("P7",il.P7,"p7"),Segment("index 2",RPI_INDEXES[1],"cbc"),Segment("RTP arm",RPI_TAIL,"r2")]),
    ]


def _three_prime_ligation(substrate:list[Segment])->AdapterLigationGeometry:
    sc=Scene()
    sc.strand("rna",substrate,label="rRNA-depleted RNA (native 3′-OH)")
    sc.strand("adapter",three_prime_adapter_segments(),label="mmc1 3′ adapter",mod5="App",mod3="ddC")
    sc.ligation_end("rna",substrate[-1].name,kind="single-strand",polarity="3-prime",owner="substrate")
    sc.ligation_end("adapter","3′-adapter delimiter",kind="single-strand",polarity="5-prime",owner="adapter")
    return AdapterLigationGeometry.known("single-strand",sc,statement="mmc1 prints the 3′ adapter as a 5′-adenylated (/5App/) DNA oligo blocked by 3′ dideoxy-C, ligated to native RNA 3′ ends before fragmentation; the ligase is not named.")


def _five_prime_ligation(substrate:list[Segment])->AdapterLigationGeometry:
    sc=Scene()
    sc.strand("adapter",[Segment("5′ adapter",ADAPTER_5P_RNA)],label="mmc1 5′ RNA adapter")
    sc.strand("rna",substrate,label="PNK-phosphorylated captured fragment",mod5="P")
    sc.ligation_end("adapter","5′ adapter",kind="single-strand",polarity="3-prime",owner="adapter")
    sc.ligation_end("rna","gene fragment",kind="single-strand",polarity="5-prime",owner="substrate")
    return AdapterLigationGeometry.known("single-strand",sc,statement="mmc1 prints the 5′ adapter as RNA, ligated after PNK phosphorylation and 500-1000-nt gel purification; the ligase and the adapter's 3′-end chemistry are not stated.")


def _rna_state(segments, *, label="single-stranded RNA with DNA 3′ adapter", mod5=""):
    sc=Scene()
    sc.strand("rna",list(segments),label=label,mod5=mod5,mod3="ddC")
    return sc


def protocol_graph():
    authority=SOURCE+" Reaction order, reagents and oligos: "+SI+"; run: "+ARTICLE+"."
    graph=ProtocolGraph(TITLE,authority=authority)
    total=graph.input("total_rna","total RNA (TRIzol, RNeasy MinElute; >200 nt fraction)",MaterialKind.BIOLOGICAL,state=molecule("total RNA",(role("rRNA"),role("other RNA with native 3-prime ends"))),primary=True,source=ARTICLE)
    ribozero=graph.input("ribozero","Ribo-Zero rRNA removal kit (Epicentre MRZH11124)",MaterialKind.SELECTION,source=SI)
    rna=graph.product("rna","rRNA-depleted RNA with native 3-prime ends",MaterialKind.MOLECULAR,state=molecule("rRNA-depleted RNA",(role("gene body"),role("native poly(A) and terminal additions"))))
    rrna=graph.product("rrna","removed rRNA",MaterialKind.WASTE)
    graph.react("deplete","Deplete rRNA",StepKind.SELECTION,((total,InputRole.TEMPLATE),(ribozero,InputRole.SELECTION)),(rna,rrna))
    graph.discard(rrna,"removed by Ribo-Zero before 3-prime-adapter ligation")
    adapter3=graph.input("three_prime_adapter","3′ adapter /5App/CTGAC-N15-TGGAATTCTCGGGTGCCAAGGC/iBiodT//iBiodT//3ddC/",MaterialKind.ADAPTER,state=_rna_state(three_prime_adapter_segments(),label="single-stranded DNA 3′ adapter",mod5="App"),source=SI)
    ligase3=graph.input("three_prime_ligase","3′-adapter ligase (not named in the source)",MaterialKind.CATALYST)
    native=[Segment("gene RNA","X"*20,placeholder=True),Segment("native poly(A) tail (variable length)","A"*12,placeholder=True),Segment("native 3′ end (variable length)","X"*3,placeholder=True)]
    tagged_con=Construct([*native,*three_prime_adapter_segments()],name="biotin-adapter-tagged RNA")
    tagged=graph.product("tagged","biotin-adapter-tagged native RNA termini",MaterialKind.MOLECULAR,state=_rna_state(tagged_con))
    graph.react("ligate_three_prime","Ligate the biotinylated 3-prime adapter",StepKind.LIGATION,((rna,InputRole.TEMPLATE),(adapter3,InputRole.STRUCTURAL),(ligase3,InputRole.CATALYST)),(tagged,),geometry=_three_prime_ligation(native),source=SI)
    rnase=graph.input("rnase_t1","RNase T1 at low concentration (Ambion AM2283; cleaves after G, sparing poly(A))",MaterialKind.CATALYST,source=SI)
    fragments=graph.product("fragments","partially digested adapter-bearing terminal RNA fragments",MaterialKind.MOLECULAR,state=tagged.state)
    graph.react("digest","Partially digest RNA with RNase T1",StepKind.DIGESTION,((tagged,InputRole.TEMPLATE),(rnase,InputRole.CATALYST)),(fragments,),source=SI)
    beads=graph.input("streptavidin","streptavidin beads (Invitrogen 11206D)",MaterialKind.SELECTION,source=SI)
    captured=graph.product("captured","captured native-terminal RNA fragments",MaterialKind.MOLECULAR,state=fragments.state)
    rejected=graph.product("rejected","non-terminal RNA fragments",MaterialKind.WASTE)
    graph.react("capture","Retain biotin-linked native terminal fragments",StepKind.AFFINITY_CAPTURE,((fragments,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(captured,rejected),source=SI)
    graph.discard(rejected,"lacks the biotinylated native 3-prime-end adapter")
    repair=graph.input("end_repair","polynucleotide kinase (Takara 2021B)",MaterialKind.CATALYST,source=SI)
    phosphorylated=graph.product("phosphorylated","5-prime-phosphorylated captured RNA",MaterialKind.MOLECULAR,state=captured.state)
    graph.react("repair","Phosphorylate the selected RNA with PNK",StepKind.PHOSPHORYLATION,((captured,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(phosphorylated,),source=SI)
    gel=graph.input("gel_selection","gel purification, 500-1000 nt",MaterialKind.SELECTION,source=SI)
    repaired=graph.product("repaired","500-1000-nt captured RNA",MaterialKind.MOLECULAR,state=captured.state)
    offsize=graph.product("offsize","RNA outside 500-1000 nt",MaterialKind.WASTE)
    graph.react("size_select","Gel-purify 500-1000 nt",StepKind.SELECTION,((phosphorylated,InputRole.TEMPLATE),(gel,InputRole.SELECTION)),(repaired,offsize),source=SI)
    graph.discard(offsize,"outside the 500-1000-nt gel window")
    adapter5=graph.input("five_prime_adapter","5′ RNA adapter",MaterialKind.ADAPTER,state=ADAPTER_5P_RNA,source=SI)
    ligase5=graph.input("five_prime_ligase","5′-adapter ligase (not named in the source)",MaterialKind.CATALYST)
    selected=[Segment("gene fragment","X"*20,placeholder=True),Segment("native poly(A) tail (variable length)","A"*12,placeholder=True),Segment("native 3′ end (variable length)","X"*3,placeholder=True),Segment("3′ adapter","X"*8,placeholder=True)]
    dual_con=Construct([Segment("5′ adapter",ADAPTER_5P_RNA,"r1"),*tagged_con],name="two-adapter terminal RNA")
    dual=graph.product("dual_adapter","two-adapter terminal RNA",MaterialKind.MOLECULAR,state=_rna_state(dual_con))
    graph.react("ligate_five_prime","Ligate the 5-prime RNA adapter",StepKind.LIGATION,((repaired,InputRole.TEMPLATE),(adapter5,InputRole.STRUCTURAL),(ligase5,InputRole.CATALYST)),(dual,),geometry=_five_prime_ligation(selected),source=SI)
    rt_primer=graph.input("rt_primer","RT primer",MaterialKind.OLIGO,state=RT_PRIMER,source=SI)
    rt=graph.input("reverse_transcriptase","reverse transcriptase (Invitrogen 18080085; conditions not stated)",MaterialKind.CATALYST,source=SI)
    cdna=graph.product("cdna","TAIL-seq cDNA",MaterialKind.MOLECULAR,state=molecule("TAIL-seq reverse-transcription product",(role("cDNA primed by the printed RT primer"),role("RNA-template recovery unresolved"))))
    graph.react("reverse_transcribe","Reverse-transcribe from the 3′-adapter RT primer",StepKind.REVERSE_TRANSCRIPTION,((dual,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,),source=SI)
    fwd=graph.input("pcr_forward","PCR primer (forward)",MaterialKind.OLIGO,state=PCR_FORWARD,source=SI)
    rev=graph.input("pcr_reverse","PCR primer (reverse) 2, a representative of the nine indexed primers (per-sample member not stated)",MaterialKind.OLIGO,state=pcr_reverse(2),source=SI)
    pol=graph.input("pcr_polymerase","Phusion DNA polymerase (Thermo F-530L; cycles not stated)",MaterialKind.CATALYST,source=SI)
    library=graph.product("library",FINAL_LIBRARY.name+" (PCR products purified again)",MaterialKind.LIBRARY,state=FINAL_LIBRARY,source=SI)
    graph.react("library_pcr","Amplify with the forward and indexed reverse primers",StepKind.PCR,((cdna,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),(rev,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(cdna.state,reason="The cDNA intermediate is not drawn base-for-base; the library is assembled from the printed adapter and primer sequences."),source=SI)
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=SEQ_PRIMERS,acquired_roles=("Read 1","Read 2"),read_lengths=READ_LENGTHS,read_uses={"Read 1":"include","Read 2":"extract_then_include"}))
    graph.declare_reviewed_scope(["total_rna","ribozero","deplete","gel_selection","size_select","three_prime_adapter","three_prime_ligase","ligate_three_prime","rnase_t1","digest","streptavidin","capture","end_repair","repair","five_prime_adapter","five_prime_ligase","ligate_five_prime","rt_primer","reverse_transcriptase","reverse_transcribe","pcr_forward","pcr_reverse","pcr_polymerase","library_pcr","library"],authority)
    return graph.require_publishable()
