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
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind, Terminal)

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
# oBW168 (Table S3-Standard oligos): T7 promoter TAATACGACTCACTATAG, one G, s7, then AG
# (the first two mosaic-end bases). An earlier model wrote GGG+s7 and so had one extra G.
T7_S7="TAATACGACTCACTATAGG"+nx.S7+"AG"
SMART_PCR=rt.SMART_HANDLE  # oDS028 AAGCAGTGGTATCAACGCAGAGT
TN5_OLIGOS=(nx.ME_RC,nx.S7+nx.ME)  # oDS029 (/5Phos/) and oDS030

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
  f"oDS028 S³-handle PCR primer — 5′-{SMART_PCR}-3′",
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
 wf.react("Lyse cells, strip RNA and synthesize the second cDNA strand",[Row(chunks=[("cDNA:RNA hybrid → RNase H → single-stranded cDNA → Klenow exo− + oBW140 dsDNA",None,False)])],note="Silane cleanup, RNase H and 1.8× SPRI precede bulk second-strand synthesis.")
 wf.react("Amplify the full-length GEX library",[Row(chunks=[("dsDNA → Q5 PCR with oBW170 + oDS028, 14 cycles",None,False)])],note="Two 1.2× SPRI cleanups remove primer dimers.")
 wf.react("Tagment, PCR-add a T7 promoter and transcribe",[Row(chunks=[("dsDNA → Tn5 fragment → T7-promoter PCR → amplified RNA","me",False)])],note="T7 transcription linearly amplifies indexed library molecules.")
 wf.react("Deplete rRNA library molecules",[Row(chunks=[("RNA library + rRNA DNA probes → RNase H cleavage",None,False)])],note="Post hoc depletion acts on amplified RNA rather than intact cells.")
 wf.react("Reverse-transcribe depleted RNA and add sequencing adapters",[Row(chunks=[("depleted RNA → cDNA → final adapter PCR",None,False)])],note="The final RT is primed from the P5 end (oBW170 + oBW171), so only full-length constructs are copied; a 9-cycle oBW170 + P7-index PCR follows.")
 return [(r.action,r.after.rows,r.note) for r in wf.reactions]


def _spri(g,ident,template,ratio,label,state=None,kind=MaterialKind.MOLECULAR):
 beads=g.input(ident+"_beads",ratio+" SPRI",MaterialKind.SELECTION)
 kept=g.product(ident+"_kept",label,kind,state=state)
 waste=g.product(ident+"_waste",ratio+" SPRI supernatant",MaterialKind.WASTE)
 g.react(ident,"Purify with "+ratio+" SPRI",StepKind.CLEANUP,
         ((template,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(kept,waste))
 g.discard(waste,"unbound material removed by the bead cleanup")
 return kept


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input("fixed_cells","fixed, ethanol- and lysozyme-permeabilized bacterial cells",MaterialKind.BIOLOGICAL,
               state=INITIAL_ROWS,primary=True)
 rt_oligo=g.input("indexed_rt_primer","/5Phos/ BC1/UMI random RT primer (Barcode 1 shown, Table S3-RT oligos)",MaterialKind.OLIGO,
                   state="".join(x.top for x in rt_primer()))
 rt=g.input("reverse_transcriptase","Maxima H Minus reverse transcriptase (7.5% PEG 8000)",MaterialKind.CATALYST)
 cdna=g.product("bc1_cdna","BC1- and UMI-tagged bacterial cDNA:RNA hybrids in cells",MaterialKind.MOLECULAR,
                state=tuple(rt_scene().rows()))
 g.react("in_situ_rt","Index fixed bacteria by random-primed in-situ RT (stepped 8–50 °C)",
         StepKind.REVERSE_TRANSCRIPTION,
         ((cells,InputRole.TEMPLATE),(rt_oligo,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),
         (cdna,))
 droplet=g.input("droplet_condition","Chromium scATAC v1.1 Chip H with Ampligase buffer, reducing agent B and oDS025",
                 MaterialKind.CONDITION)
 occupied=g.product("occupied_droplets","droplets containing BC1-tagged cells and beads",
                    MaterialKind.MOLECULAR)
 empty=g.product("excluded_droplets","empty or unusable droplets",MaterialKind.WASTE)
 g.react("partition","Pool, then partition indexed cells into 10x ATAC droplets",StepKind.PARTITION,
         ((cdna,InputRole.TEMPLATE),(droplet,InputRole.CONDITION)),(occupied,empty),
         changes_molecule=False)
 g.discard(empty,"empty and unusable droplets do not yield retained molecules")
 bc2=g.input("bc2_adapter","commercial 10x ATAC BC2 gel-bead oligo",MaterialKind.ADAPTER,
             inferred=True,source="10x CG000209 v1.1 guide and Cell Ranger ATAC barcode whitelist; experiment-specific barcode assignment")
 bridge=g.input("bridge_oligo","oDS025 bridge (3′-ddC)",MaterialKind.OLIGO,state=BRIDGE)
 ligase=g.input("dna_ligase","Ampligase thermostable DNA ligase (2.3 U/µl)",MaterialKind.CATALYST)
 dual=g.product("dual_barcoded_cdna","BC1/BC2-tagged bacterial cDNA",MaterialKind.MOLECULAR,
                state=ligated_first_strand())
 g.react("ligate_bc2","Ligate the second cell index in droplets: 12 × (98 °C 30 s, 59 °C 2 min)",StepKind.LIGATION,
         ((occupied,InputRole.TEMPLATE),(bc2,InputRole.STRUCTURAL),(bridge,InputRole.STRUCTURAL),
          (ligase,InputRole.CATALYST)),(dual,),
         geometry=AdapterLigationGeometry.known(
             "splinted",ligation_ready_scene(),ligation_scene(),
             statement="Supplementary Table 3 prints oDS025; its two arms are the reverse complement of s5 plus the first 21 bases of the phosphorylated RT primer."))
 lysis=g.input("lysis","recovery agent; Lysis-T + proteinase K, 55 °C 1 h",MaterialKind.CONDITION)
 lysed=g.product("lysed_hybrids","released DNA:RNA hybrid library",MaterialKind.MOLECULAR,state=dual.state)
 g.react("break_and_lyse","Break emulsions, reverse-crosslink and lyse cells",StepKind.DEPROTEINIZATION,
         ((dual,InputRole.TEMPLATE),(lysis,InputRole.CONDITION)),(lysed,))
 silane=g.input("silane","Dynabeads MyOne Silane in 10x cleanup buffer; 2 × 80% ethanol",MaterialKind.SELECTION)
 hybrids=g.product("isolated_hybrids","silane-purified DNA:RNA hybrid library",MaterialKind.MOLECULAR,state=dual.state)
 silane_waste=g.product("silane_waste","silane bead supernatant and washes",MaterialKind.WASTE)
 g.react("silane_cleanup","Isolate hybrid libraries on silane beads",StepKind.CLEANUP,
         ((lysed,InputRole.TEMPLATE),(silane,InputRole.SELECTION)),(hybrids,silane_waste))
 g.discard(silane_waste,"unbound material removed by the silane cleanup")
 rnaseh_strip=g.input("rnase_h_strip","RNase H (NEB M0297), 37 °C 30 min",MaterialKind.CATALYST)
 single=g.product("single_stranded_cdna","single-stranded dual-barcoded cDNA",MaterialKind.MOLECULAR,state=dual.state)
 rna_frag=g.product("hybrid_rna_fragments","RNase-H-degraded template RNA",MaterialKind.WASTE)
 g.react("strip_rna","Strip template RNA with RNase H",StepKind.DIGESTION,
         ((hybrids,InputRole.TEMPLATE),(rnaseh_strip,InputRole.CATALYST)),(single,rna_frag))
 g.discard(rna_frag,"degraded RNA is removed by the following 1.8× SPRI")
 single=_spri(g,"spri_after_rnaseh",single,"1.8×","purified single-stranded cDNA",dual.state)
 random_primer=g.input("random_second_strand_primer","S³ randomer oBW140",MaterialKind.OLIGO,state=S3_RANDOMER)
 klenow=g.input("klenow","Klenow fragment (3′→5′ exo−), PEG 8000",MaterialKind.CATALYST)
 ds=g.product("double_stranded_cdna","double-stranded dual-indexed cDNA",MaterialKind.MOLECULAR,
              state=amplified_cdna())
 g.react("second_strand","Synthesize the second cDNA strand with S3 randomer oBW140, 37 °C 60 min",StepKind.EXTENSION,
         ((single,InputRole.TEMPLATE),(random_primer,InputRole.STRUCTURAL),
          (klenow,InputRole.CATALYST)),(ds,))
 ds=_spri(g,"spri_after_second_strand",ds,"1.8×","purified double-stranded cDNA",amplified_cdna())
 gex_p5=g.input("gex_forward","oBW170 P5 primer",MaterialKind.OLIGO,state=P5_PRIMER)
 gex_smart=g.input("gex_reverse","oDS028 S³-handle primer",MaterialKind.OLIGO,state=SMART_PCR)
 q5=g.input("q5","Q5 High-Fidelity 2× master mix",MaterialKind.CATALYST)
 gex=g.product("gex_library","full-length GEX library",MaterialKind.MOLECULAR,state=amplified_cdna())
 g.react("gex_pcr","Amplify the full-length library with oBW170 + oDS028 (14 cycles)",StepKind.PCR,
         ((ds,InputRole.TEMPLATE),(gex_p5,InputRole.STRUCTURAL),(gex_smart,InputRole.STRUCTURAL),
          (q5,InputRole.CATALYST)),(gex,),
         geometry=PCRGeometry.exact(amplified_cdna(),PCRPrimer("oBW170 P5",P5_PRIMER),
             PCRPrimer("oDS028",SMART_PCR),
             statement="Methods, second-strand cDNA synthesis: Q5 PCR with oDS028 and oBW170, 14 cycles."))
 gex=_spri(g,"spri_gex_1",gex,"1.2×","once-purified GEX library",amplified_cdna())
 gex=_spri(g,"spri_gex_2",gex,"1.2×","GEX library (primer dimers removed)",amplified_cdna())
 tn5_adapter=g.input("tn5_adapter","Read 2-only transposome adapter oDS029 (/5Phos/) : oDS030",MaterialKind.ADAPTER,
                     state=TN5_OLIGOS)
 tn5=g.input("tn5","EZ-Tn5 transposase (Lucigen) in TAPS/MgCl2 buffer, 55 °C 10 min",MaterialKind.CATALYST)
 selected=tagmented_template()
 tagged=g.product("tagmented_cdna","Tn5-tagged GEX fragments with P5/BC end and Read 2 end",MaterialKind.MOLECULAR,
                  state=selected)
 tag_other=g.product("tagment_other","internal fragments with Read 2 adapters at both ends",MaterialKind.MOLECULAR)
 g.react("tagment","Tagment 10 ng GEX library with Read 2-only Tn5; SDS inactivation",StepKind.TAGMENTATION,
         ((gex,InputRole.TEMPLATE),(tn5_adapter,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),
         (tagged,tag_other))
 g.mark_terminal(Terminal(tag_other.id,Disposition.UNRESOLVED,
   "lack the P5 end so oBW170 cannot extend them; carried into the PCR rather than physically removed"))
 tagged=_spri(g,"spri_after_tn5",tagged,"1.2×","fragmented GEX library",selected)
 t7_fwd=g.input("t7_forward","oBW170 P5 primer",MaterialKind.OLIGO,state=P5_PRIMER)
 t7_rev=g.input("t7_reverse","oBW168 T7–Nextera primer",MaterialKind.OLIGO,state=T7_S7)
 t7_template=g.product("t7_template","T7-promoter-bearing indexed DNA",MaterialKind.MOLECULAR,
                       state=selected)
 g.react("t7_pcr","PCR-install the T7 promoter (72 °C 3 min gap fill, 9 cycles)",StepKind.PCR,
         ((tagged,InputRole.TEMPLATE),(t7_fwd,InputRole.STRUCTURAL),
          (t7_rev,InputRole.STRUCTURAL),(q5,InputRole.CATALYST)),(t7_template,),
         geometry=PCRGeometry.exact(
             selected,PCRPrimer("oBW170 P5",P5_PRIMER),
             PCRPrimer("oBW168 T7-s7",T7_S7,anneal_nt=len(nx.S7)+2),
             statement="The Methods and Supplementary Table 3 specify oBW170 and oBW168 for this PCR."))
 t7_template=_spri(g,"spri_after_t7_pcr",t7_template,"1.2×","purified T7 template",selected)
 t7=g.input("t7_polymerase","T7 RNA polymerase (Thermo EP0112), 37 °C 2 h",MaterialKind.CATALYST)
 transcribed=g.product("transcribed_mix","IVT reaction",MaterialKind.MOLECULAR)
 g.react("transcribe","Linearly amplify 100 ng by T7 transcription",
         StepKind.TRANSCRIPTION,
         ((t7_template,InputRole.TEMPLATE),(t7,InputRole.CATALYST)),(transcribed,))
 dnase=g.input("dnase","DNase I (NEB M0303)",MaterialKind.CATALYST)
 amplified_rna=g.product("amplified_rna","linearly amplified indexed RNA",MaterialKind.MOLECULAR)
 dna_waste=g.product("ivt_template_digest","DNase-I-digested DNA template",MaterialKind.WASTE)
 g.react("dnase_ivt","Digest the DNA template with DNase I",StepKind.DIGESTION,
         ((transcribed,InputRole.TEMPLATE),(dnase,InputRole.CATALYST)),(amplified_rna,dna_waste))
 g.discard(dna_waste,"digested DNA removed by the following 2× SPRI")
 amplified_rna=_spri(g,"spri_after_ivt",amplified_rna,"2×","in vitro transcribed GEX library")
 rrna_probe=g.input("rrna_probe","rRNA depletion DNA probes (3 µg; Table S3-rRNA depletion oligos)",MaterialKind.OLIGO)
 hybrid=g.product("rrna_hybrid","probe-hybridized amplified RNA",MaterialKind.MOLECULAR)
 g.react("hybridize_rrna","Hybridize probes to 500 ng RNA (95 °C, ramp to 25 °C)",
         StepKind.HYBRIDIZATION,
         ((amplified_rna,InputRole.TEMPLATE),(rrna_probe,InputRole.STRUCTURAL)),(hybrid,))
 rnaseh=g.input("rnase_h","thermostable RNase H + RNase H, 50 °C 45 min",MaterialKind.CATALYST)
 depleted=g.product("depleted_rna","rRNA-depleted amplified RNA with probes",MaterialKind.MOLECULAR)
 rrna_frag=g.product("rrna_fragments","RNase-H-cleaved rRNA-derived library RNA",MaterialKind.MOLECULAR)
 g.react("deplete_rrna","Cleave rRNA–probe hybrids with RNase H",StepKind.DIGESTION,
         ((hybrid,InputRole.TEMPLATE),(rnaseh,InputRole.CATALYST)),(depleted,rrna_frag))
 g.mark_terminal(Terminal(rrna_frag.id,Disposition.UNRESOLVED,
   "cleaved pieces are not shown as separated; only those keeping both library ends can be RT-primed and PCR-amplified"))
 probe_digest=g.product("probe_digest","DNase-I-digested probes",MaterialKind.WASTE)
 clean_rna=g.product("probe_free_rna","rRNA-depleted RNA library",MaterialKind.MOLECULAR)
 g.react("dnase_probes","Degrade DNA probes with DNase I, 37 °C 45 min",StepKind.DIGESTION,
         ((depleted,InputRole.TEMPLATE),(dnase,InputRole.CATALYST)),(clean_rna,probe_digest))
 g.discard(probe_digest,"digested probes removed by the following 2× SPRI")
 clean_rna=_spri(g,"spri_after_depletion",clean_rna,"2×","purified rRNA-depleted RNA library")
 final_rt_primer=g.input("final_rt_primer","oBW170 P5 primer (end-specific RT primer)",MaterialKind.OLIGO,state=P5_PRIMER)
 final_rt_short=g.input("final_rt_short_primer","oBW171 short P5 primer",MaterialKind.OLIGO,state=SHORT_P5)
 final_cdna=g.product("final_cdna","depleted second-generation cDNA",MaterialKind.MOLECULAR)
 g.react("final_rt","Reverse-transcribe full-length depleted RNA from its P5 end",
         StepKind.REVERSE_TRANSCRIPTION,
         ((clean_rna,InputRole.TEMPLATE),(final_rt_primer,InputRole.STRUCTURAL),
          (final_rt_short,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(final_cdna,))
 final_cdna=_spri(g,"spri_after_final_rt",final_cdna,"1.2×","purified final cDNA")
 fwd=g.input("library_forward","oBW170 P5 primer",MaterialKind.OLIGO,state=P5_PRIMER)
 rev=g.input("library_reverse","unique P7 index primer (oBW176 shown)",MaterialKind.OLIGO,
             state=il.P7+I7_PRIMER_INDEX+nx.S7)
 library=g.product("library_pcr_product",FINAL_LIBRARY.name+" (PCR product)",MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
 g.react("library_pcr","Index the depleted cDNA by PCR (9 cycles)",StepKind.PCR,
         ((final_cdna,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
          (rev,InputRole.STRUCTURAL),(q5,InputRole.CATALYST)),(library,),
         geometry=PCRGeometry.exact(
             FINAL_LIBRARY,PCRPrimer("oBW170 P5",P5_PRIMER),
             PCRPrimer("oBW176 indexed P7",il.P7+revcomp(I7_EXAMPLE)+nx.S7,
                       anneal_nt=len(nx.S7)),
             statement="The final PCR uses oBW170 and one printed P7-index-s7 primer; oBW176 is shown as the representative index."))
 library=_spri(g,"spri_final_1",library,"0.8×","once-purified final library",FINAL_LIBRARY)
 library=_spri(g,"spri_final_2",library,"0.8×",FINAL_LIBRARY.name,FINAL_LIBRARY,kind=MaterialKind.LIBRARY)
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),
                                   acquired_roles=tuple(READ_LENGTHS),read_lengths=READ_LENGTHS,
                                   read_uses=READ_USES))
 g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
 return g.require_publishable()
