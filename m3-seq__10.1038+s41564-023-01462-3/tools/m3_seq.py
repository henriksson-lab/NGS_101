"""M3-seq molecular model from the paper's Supplementary Table 3."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct,MolecularState,Row,Scene,Segment,Workflow,feature,revcomp
from pcr import PCRGeometry,PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE="M3-seq — multiplexed microbial single-cell RNA-seq"
NOTES="01_m3-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41564-023-01462-3">McNulty et al., <i>Nature Microbiology</i> (2023)</a>.'
SUMMARY="Random in-situ RT supplies a plate barcode and UMI; 10x ATAC droplets ligate a second barcode, followed by Tn5/T7 amplification and post hoc RNase-H rRNA depletion."
CAVEAT="The published RT, bridge, amplification, Tn5 and indexing oligos are shown exactly. The 10x ATAC vendor distributes the 16-base barcode whitelist; BC2 remains a role token because the barcode acquired by a particular cell varies by experiment."

UMI_FEATURE=feature("m3seq_umi","umi","random")
BC1_FEATURE=feature("m3seq_bc1","cell_barcode","combinatorial",group="m3seq_cell",part="RT plate")
BC2_FEATURE=feature("m3seq_bc2","cell_barcode","combinatorial",group="m3seq_cell",part="10x droplet")
I7_FEATURE=feature("m3seq_i7","sample_index","fixed")
BC1_EXAMPLE="AAGTGATTAGCAA"
I7_PRIMER_INDEX="TGCCTCTT"
I7_EXAMPLE=revcomp(I7_PRIMER_INDEX)
BRIDGE="CGTCGTGTAGGGAAAGAGTGTGACGCTGCCGACGA"
P5_PRIMER="AATGATACGGCGACCACCGAGA"
SHORT_P5="ACGGCGACCACCGAGA"
SMART_PRIMER=rt.SMART_HANDLE
S3_RANDOMER=rt.SMART_HANDLE+"GAATNNNGGNNN"
T7_S7="TAATACGACTCACTATAGGG"+nx.S7+"AG"

def seg(name,top,tag=None,**kw):
 kw.setdefault("feature",{"UMI":UMI_FEATURE,"BC1":BC1_FEATURE,
                           "BC2":BC2_FEATURE,"i7":I7_FEATURE}.get(name))
 return Segment(name,top,tag,**kw)

def rt_primer():
 return [seg("TruSeq Read 1",il.TRUSEQ_READ1,"r1"),
         seg("UMI","N"*8,"umi",placeholder=True),
         seg("BC1",BC1_EXAMPLE,"cbc"),
         seg("random hexamer","N"*6,placeholder=True)]

def bead_oligo():
 return [seg("P5",il.P5,"p5"),seg("BC2","C"*16,"cbc",placeholder=True,inferred=True),
         seg("s5",nx.S5,"s5")]

def bridge_oligo():
 return [seg("RT-primer bridge",BRIDGE[:21]),seg("bead bridge",BRIDGE[21:])]

def rt_scene():
 sc=Scene();sc.strand("RNA",[seg("bacterial RNA","X"*34,placeholder=True)],label="fixed bacterial RNA")
 sc.anneal("RT primer",rt_primer(),to="RNA",pair=("random hexamer","bacterial RNA"),
           label="published Barcode 1 example",mod5="p",
           unpaired=("TruSeq Read 1","UMI","BC1","random hexamer"))
 sc.arrow("RT primer","in-situ reverse transcription");sc.labels("RT primer")
 return sc

def ligated_first_strand(insert_nt=32):
 return Construct([*bead_oligo(),*rt_primer(),seg("cDNA","X"*insert_nt,placeholder=True)],
                  "M3-seq dual-barcoded first strand")

def ligation_ready_scene():
 sc=Scene();sc.strand("bead oligo",bead_oligo(),label="10x scATAC gel-bead oligo")
 sc.strand("RT product",[*rt_primer(),seg("cDNA","X"*24,placeholder=True)],
           label="phosphorylated round-one RT product",mod5="p")
 sc.strand("bridge",bridge_oligo(),label="oDS025 bridge",mod3="ddC")
 sc.ligation_end("bead oligo","s5",kind="splinted",polarity="3-prime",owner="adapter",
                 text="bead s5 3′-OH at the nick")
 sc.ligation_end("RT product","TruSeq Read 1",kind="splinted",polarity="5-prime",owner="substrate",
                 text="RT-primer 5′ phosphate at the nick")
 return sc

def ligation_scene():
 joined=[*bead_oligo(),*rt_primer(),seg("cDNA","X"*24,placeholder=True)]
 sc=Scene();sc.strand("product",joined,label="BC2–BC1 ligation product")
 sc.anneal("bridge",bridge_oligo(),to="product",pair=("RT-primer bridge","TruSeq Read 1"),
           label="oDS025; 3′-ddC",above=True)
 if revcomp(BRIDGE)!=nx.S5+il.TRUSEQ_READ1[:21]:
  raise ValueError("oDS025 no longer bridges s5 directly to the RT-primer handle")
 sc.junction("product","s5","TruSeq Read 1")
 return sc

def amplified_cdna(insert_nt=32):
 return Construct([*bead_oligo(),*rt_primer(),seg("cDNA","X"*insert_nt,placeholder=True),
                   seg("S3 randomer handle reverse complement",revcomp(SMART_PRIMER),"tso")],
                  "M3-seq full-length amplified cDNA")

def tagmented_template(insert_nt=32):
 return Construct([*bead_oligo(),*rt_primer(),seg("cDNA","X"*insert_nt,placeholder=True),
                   seg("mosaic end reverse complement",nx.ME_RC,"me"),
                   seg("s7 reverse complement",nx.S7_RC,"s7")],
                  "M3-seq selected Read-2-only Tn5 fragment")

def final_library(insert_nt=32):
 lib=Construct([*bead_oligo(),*rt_primer(),seg("cDNA","X"*insert_nt,placeholder=True),
                seg("mosaic end reverse complement",nx.ME_RC,"me"),
                seg("s7 reverse complement",nx.S7_RC,"s7"),
                seg("i7",I7_EXAMPLE,"cbc"),seg("P7 reverse complement",il.P7_RC,"p7")],
               "M3-seq RNA library")
 problems=sp.verify(lib,SEQ_PRIMERS)
 if problems: raise ValueError("invalid M3-seq final construct: "+"; ".join(problems))
 return lib

INDEX2_FROM_P5=sp.custom("Index 2 (i5)","flow-cell P5 oligo (forward-strand workflow)",
                         il.P5,"Illumina forward-strand indexed sequencing",
                         "The 30-cycle i5 read traverses s5 and the 16-base droplet barcode.")
SEQ_PRIMERS=(sp.TRUSEQ["R1"],sp.NEXTERA["I1"],INDEX2_FROM_P5,sp.NEXTERA["R2"])
FINAL_LIBRARY=final_library()
FINAL_CAPTION="The final cDNA library retains the round-one barcode, round-two barcode and UMI through T7 amplification and depletion."
SEQUENCING_INTRO="Read cycles recover the combinatorial cell identity, UMI and bacterial transcript insert."
READ_LENGTHS={"Read 1":26,"Index 1 (i7)":8,"Index 2 (i5)":30,"Read 2":74}
READ_USES={"Read 1":"extract","Index 1 (i7)":"demultiplex",
           "Index 2 (i5)":"extract","Read 2":"include"}
SHOW_READ_CONTENT=True
INITIAL_ROWS=(Row(chunks=[("fixed bacterial cells containing total RNA",None,False)]),)
INITIAL_NAME="Fixed bacterial cells"

def oligos():
 return (
  f"Barcode 1 RT primer — 5′-/5Phos/{il.TRUSEQ_READ1}{'N'*8}{BC1_EXAMPLE}{'N'*6}-3′",
  f"oDS025 bridge — 5′-{BRIDGE}/3ddC/-3′",
  f"oBW140 S³ randomer — 5′-{S3_RANDOMER}-3′",
  f"oBW170 P5 primer — 5′-{P5_PRIMER}-3′",
  f"oBW171 short P5 primer — 5′-{SHORT_P5}-3′",
  f"oDS029 Tn5 mosaic end — 5′-/5Phos/{nx.ME_RC}-3′",
  f"oDS030 Read-2-only transposon strand — 5′-{nx.S7}{nx.ME}-3′",
  f"oBW168 T7–s7 primer — 5′-{T7_S7}-3′",
  f"oBW176 representative indexed P7–s7 primer — 5′-{il.P7}{I7_PRIMER_INDEX}{nx.S7}-3′",
 )

def sections():
 wf=Workflow(MolecularState("Fixed bacteria",(Row(chunks=[("fixed bacterial RNA",None,False)]),)))
 wf.react("Index fixed bacteria by random-primed in-situ RT",rt_scene().rows(),note="The published 5′-phosphorylated RT primer carries an 8-base UMI, a 13-base BC1 and a 3′ random hexamer.")
 wf.ligate_adapter("Ligate the second index in 10x ATAC droplets",AdapterLigationGeometry.known("splinted",ligation_ready_scene(),ligation_scene(),statement="Published oDS025 bridges the 10x s5 arm to the phosphorylated TruSeq Read-1 end of the RT product."),ligation_scene().rows(),note="BC2 is ligated to cell-associated cDNA; BC1+BC2 defines a cell.")
 wf.react("Lyse cells and synthesize the second cDNA strand",[Row(chunks=[("single-stranded indexed cDNA → Klenow/random-primer dsDNA",None,False)])],note="Bulk second-strand synthesis makes the library amplifiable.")
 wf.react("Tagment, PCR-add a T7 promoter and transcribe",[Row(chunks=[("dsDNA → Tn5 fragment → T7-promoter PCR → amplified RNA","me",False)])],note="T7 transcription linearly amplifies indexed library molecules.")
 wf.react("Deplete rRNA library molecules",[Row(chunks=[("RNA library + rRNA DNA probes → RNase H cleavage",None,False)])],note="Post hoc depletion acts on amplified RNA rather than intact cells.")
 wf.react("Reverse-transcribe depleted RNA and add sequencing adapters",[Row(chunks=[("depleted RNA → cDNA → final adapter PCR",None,False)])],note="A final RT and PCR return the selected molecules to an Illumina-compatible cDNA library.")
 return [(r.action,r.after.rows,r.note) for r in wf.reactions]


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input("fixed_cells","fixed permeabilized bacterial cells",MaterialKind.BIOLOGICAL,
               state=INITIAL_ROWS,primary=True)
 rt_primer=g.input("indexed_rt_primer","BC1/UMI random RT primer",MaterialKind.OLIGO)
 rt=g.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
 cdna=g.product("bc1_cdna","BC1- and UMI-tagged bacterial cDNA",MaterialKind.MOLECULAR,
                state=tuple(rt_scene().rows()))
 g.react("in_situ_rt","Index fixed bacteria by random-primed in-situ RT",
         StepKind.REVERSE_TRANSCRIPTION,
         ((cells,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),
         (cdna,))
 droplet=g.input("droplet_condition","10x ATAC droplet partitioning condition",
                 MaterialKind.CONDITION)
 occupied=g.product("occupied_droplets","droplets containing BC1-tagged cells and beads",
                    MaterialKind.MOLECULAR)
 empty=g.product("excluded_droplets","empty or unusable droplets",MaterialKind.WASTE)
 g.react("partition","Partition indexed cells into 10x ATAC droplets",StepKind.PARTITION,
         ((cdna,InputRole.TEMPLATE),(droplet,InputRole.CONDITION)),(occupied,empty),
         changes_molecule=False)
 g.discard(empty,"empty and unusable droplets do not yield retained molecules")
 bc2=g.input("bc2_adapter","commercial 10x ATAC BC2 adapter",MaterialKind.ADAPTER,
             inferred=True,source="10x CG000209 v1.1 guide and Cell Ranger ATAC barcode whitelist; experiment-specific barcode assignment")
 ligase=g.input("dna_ligase","DNA ligase",MaterialKind.CATALYST)
 dual=g.product("dual_barcoded_cdna","BC1/BC2-tagged bacterial cDNA",MaterialKind.MOLECULAR,
                state=ligated_first_strand())
 g.react("ligate_bc2","Ligate the second cell index in droplets",StepKind.LIGATION,
         ((occupied,InputRole.TEMPLATE),(bc2,InputRole.STRUCTURAL),
          (ligase,InputRole.CATALYST)),(dual,),
         geometry=AdapterLigationGeometry.known(
             "splinted",ligation_ready_scene(),ligation_scene(),
             statement="Supplementary Table 3 prints oDS025; its two arms are the reverse complement of s5 plus the first 21 bases of the phosphorylated RT primer."))
 random_primer=g.input("random_second_strand_primer","random second-strand primer",
                       MaterialKind.OLIGO)
 klenow=g.input("klenow","Klenow DNA polymerase",MaterialKind.CATALYST)
 ds=g.product("double_stranded_cdna","double-stranded dual-indexed cDNA",MaterialKind.MOLECULAR)
 g.react("second_strand","Synthesize the second cDNA strand with S3 randomer oBW140",StepKind.EXTENSION,
         ((dual,InputRole.TEMPLATE),(random_primer,InputRole.STRUCTURAL),
          (klenow,InputRole.CATALYST)),(ds,))
 tn5_adapter=g.input("tn5_adapter","Tn5 mosaic-end adapter",MaterialKind.ADAPTER)
 tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
 selected=tagmented_template()
 tagged=g.product("tagmented_cdna","Tn5-tagged double-stranded cDNA",MaterialKind.MOLECULAR,
                  state=selected)
 tag_waste=g.product("tagment_waste","non-amplifiable tagmentation products",MaterialKind.WASTE)
 g.react("tagment","Tagment double-stranded indexed cDNA",StepKind.TAGMENTATION,
         ((ds,InputRole.TEMPLATE),(tn5_adapter,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),
         (tagged,tag_waste))
 g.discard(tag_waste,"products lacking the required amplifiable transfer-end combination")
 t7_fwd=g.input("t7_forward","T7-promoter PCR primer",MaterialKind.OLIGO)
 t7_rev=g.input("t7_reverse","T7-amplification reverse primer",MaterialKind.OLIGO)
 pol=g.input("t7_pcr_polymerase","T7-promoter PCR polymerase",MaterialKind.CATALYST)
 t7_template=g.product("t7_template","T7-promoter-bearing indexed DNA",MaterialKind.MOLECULAR,
                       state=selected)
 g.react("t7_pcr","PCR-install the T7 promoter",StepKind.PCR,
         ((tagged,InputRole.TEMPLATE),(t7_fwd,InputRole.STRUCTURAL),
          (t7_rev,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(t7_template,),
         geometry=PCRGeometry.exact(
             selected,PCRPrimer("oBW170 P5",P5_PRIMER),
             PCRPrimer("oBW168 T7-s7",T7_S7,anneal_nt=len(nx.S7)),
             statement="The Methods and Supplementary Table 3 specify oBW170 and oBW168 for this PCR."))
 t7=g.input("t7_polymerase","T7 RNA polymerase",MaterialKind.CATALYST)
 amplified_rna=g.product("amplified_rna","linearly amplified indexed RNA",MaterialKind.MOLECULAR)
 g.react("transcribe","Linearly amplify the indexed library by T7 transcription",
         StepKind.TRANSCRIPTION,
         ((t7_template,InputRole.TEMPLATE),(t7,InputRole.CATALYST)),(amplified_rna,))
 rrna_probe=g.input("rrna_probe","rRNA depletion DNA probes",MaterialKind.OLIGO)
 hybrid=g.product("rrna_hybrid","probe-hybridized amplified RNA",MaterialKind.MOLECULAR)
 g.react("hybridize_rrna","Hybridize depletion probes to rRNA-derived molecules",
         StepKind.HYBRIDIZATION,
         ((amplified_rna,InputRole.TEMPLATE),(rrna_probe,InputRole.STRUCTURAL)),(hybrid,))
 rnaseh=g.input("rnase_h","RNase H",MaterialKind.CATALYST)
 depleted=g.product("depleted_rna","rRNA-depleted amplified RNA",MaterialKind.MOLECULAR)
 rrna_waste=g.product("rrna_waste","RNase-H-cleaved rRNA library molecules",MaterialKind.WASTE)
 g.react("deplete_rrna","Cleave rRNA-derived amplified molecules",StepKind.DIGESTION,
         ((hybrid,InputRole.TEMPLATE),(rnaseh,InputRole.CATALYST)),(depleted,rrna_waste))
 g.discard(rrna_waste,"cleaved rRNA-derived library molecules are not reverse-transcribed")
 final_rt_primer=g.input("final_rt_primer","final reverse-transcription primer",MaterialKind.OLIGO)
 final_cdna=g.product("final_cdna","depleted second-generation cDNA",MaterialKind.MOLECULAR)
 g.react("final_rt","Reverse-transcribe the depleted amplified RNA",
         StepKind.REVERSE_TRANSCRIPTION,
         ((depleted,InputRole.TEMPLATE),(final_rt_primer,InputRole.STRUCTURAL),
          (rt,InputRole.CATALYST)),(final_cdna,))
 fwd=g.input("library_forward","final library PCR forward primer",MaterialKind.OLIGO)
 rev=g.input("library_reverse","final library PCR reverse primer",MaterialKind.OLIGO)
 final_pol=g.input("library_polymerase","final library PCR polymerase",MaterialKind.CATALYST)
 library=g.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 g.react("library_pcr","Add sequencing adapters to the depleted cDNA",StepKind.PCR,
         ((final_cdna,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
          (rev,InputRole.STRUCTURAL),(final_pol,InputRole.CATALYST)),(library,),
         geometry=PCRGeometry.exact(
             FINAL_LIBRARY,PCRPrimer("oBW170 P5",P5_PRIMER),
             PCRPrimer("oBW176 indexed P7",il.P7+revcomp(I7_EXAMPLE)+nx.S7,
                       anneal_nt=len(nx.S7)),
             statement="The final PCR uses oBW170 and one printed P7-index-s7 primer; oBW176 is shown as the representative index."))
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),
                                   acquired_roles=tuple(READ_LENGTHS),read_lengths=READ_LENGTHS,
                                   read_uses=READ_USES))
 g.declare_reviewed_scope(("fixed_cells","indexed_rt_primer","reverse_transcriptase",
   "in_situ_rt","droplet_condition","partition","bc2_adapter","dna_ligase","ligate_bc2",
   "random_second_strand_primer","klenow","second_strand","tn5_adapter","tn5","tagment",
   "t7_forward","t7_reverse","t7_pcr_polymerase","t7_pcr","t7_polymerase","transcribe",
   "rrna_probe","hybridize_rrna","rnase_h","deplete_rrna","final_rt_primer","final_rt",
   "library_forward","library_reverse","library_polymerase","library_pcr","library"),SOURCE)
 return g.require_publishable()
