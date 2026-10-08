"""10x Chromium 3' v3 Feature Barcode chemistry for direct sgRNA capture."""
from __future__ import annotations
import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,revcomp

CS1 = "TTGCTAGGACCGGCCTTAAAGC"
CS2 = "CCTTAGCCGCTAATAGGTGAGC"
CS1_RC = revcomp(CS1)
CS2_RC = revcomp(CS2)
CR1 = "GTTTAAGAGCTAAGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
PARTIAL_R1N = nx.READ1_PRIMER[-22:]
TSO_DNA = rt.SMART_HANDLE + "ACAT"
FEATURE_CDNA_FORWARD = nx.READ1_PRIMER[-27:]
FEATURE_CDNA_REVERSE = rt.SMART_HANDLE[:22]
FEATURE_SI_FORWARD = il.P5 + nx.READ1_PRIMER
FEATURE_SI_REVERSE = il.TRUSEQ_READ2 + FEATURE_CDNA_REVERSE
def seg(name,top,tag=None,**kw): return Segment(name=name,top=top,tag=tag,**kw)
def bead_primer(): return [seg("Partial Read 1N",PARTIAL_R1N,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*12,"umi",placeholder=True),seg("Capture Sequence 1",CS1,"tso")]
def tso(): return [seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("ACAT","ACAT","tso"),seg("rGrGrG","GGG","tso")]
def sgrna():
    return [seg("protospacer", "N" * 20, "cbc", placeholder=True),
            seg("CR1 scaffold", CR1, "insert"),
            seg("Capture Sequence 1 rc", CS1_RC, "tso"),
            seg("Pol III terminator", "U" * 7, "insert")]
def copied_sgrna(): return [seg("CR1 complement",revcomp(CR1),"insert"),seg("protospacer complement","n"*20,"cbc",placeholder=True)]
def feature_cdna_primers(): return ([seg("Read 1N completion",FEATURE_CDNA_FORWARD,"r1")],[seg("partial TSO",FEATURE_CDNA_REVERSE,"tso")])
def feature_si_primers(): return ([seg("P5 + Read 1N",FEATURE_SI_FORWARD,"p5+r1")],[seg("Read 2 + partial TSO",FEATURE_SI_REVERSE,"r2+tso")])
def capture_scene():
    sc=Scene(); sc.strand("sgRNA",sgrna(),label="capture-sequence sgRNA")
    sc.anneal("bead primer",bead_primer(),to="sgRNA",pair=("Capture Sequence 1","Capture Sequence 1 rc"),label="gel-bead primer",mod5="bead",unpaired=("UMI",))
    sc.mark("sgRNA","protospacer","Feature Barcode / guide identity")
    sc.mark("sgRNA","Capture Sequence 1 rc","engineered capture sequence")
    sc.arrow("bead primer","reverse transcription"); return sc
def first_strand(): return Construct([*bead_primer(),*copied_sgrna(),seg("CCC",rt.UNTEMPLATED_TAIL,"tso")],name="captured sgRNA first strand")
def switch_scene():
    first=first_strand(); sc=Scene(); sc.strand("first",list(first),label="first strand")
    sc.anneal("TSO",tso(),to="first",pair=("rGrGrG","CCC"),label="TSO",above=True); sc.arrow("first","RT switches template"); return sc
def tso_complement(): return rt.UNTEMPLATED_TAIL+revcomp(TSO_DNA)
def feature_cdna(): return Construct([*bead_primer(),*copied_sgrna(),seg("TSO complement",tso_complement(),"tso")],name="CRISPR feature cDNA")
SEQ_PRIMERS=(sp.NEXTERA["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
def final_library():
    lib=Construct([seg("P5",il.P5,"p5"),seg("Nextera Read 1",nx.READ1_PRIMER,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*12,"umi",placeholder=True),seg("Capture Sequence 1",CS1,"tso"),*copied_sgrna(),seg("TSO complement",tso_complement(),"tso"),seg("Read 2 leading A","A","r2"),seg("Index 1 / Read 2 arm",il.INDEX1_PRIMER,"r2"),seg("i7 reverse complement","I"*8,"cbc",placeholder=True),seg("P7 reverse complement",il.P7_RC,"p7")],name="10x 3' v3 CRISPR Screening library")
    problems=sp.verify(lib,SEQ_PRIMERS,required_roles=("Read 1","Index 1 (i7)","Read 2"))
    if problems: raise ValueError("invalid CRISPR Screening library: "+"; ".join(problems))
    return lib
def read_layout(): return [("Read 1, 1–16","cell barcode"),("Read 1, 17–28","UMI"),("Read 2, 1–30","TSO-derived sequence"),("Read 2, 31–50","protospacer / guide identity"),("Read 2, 51 onward","CR1 scaffold")]
