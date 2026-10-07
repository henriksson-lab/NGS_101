"""Molecular construct model for Seq-Well (Gierahn et al., 2017)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, revcomp


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


LINKER = "T" * 7
SPACER = "GCCTGTCCGCGG"
BEAD_TEMPLATE = LINKER + rt.SMART_HANDLE + "AC" + "J" * 12 + "N" * 8 + "T" * 30
TSO = rt.SMART_HANDLE + "GAAT" + "GGG"
SMART_PCR = rt.SMART_HANDLE
P5_HYBRID = il.P5 + SPACER + rt.SMART_HANDLE + "AC"
CUSTOM_R1 = SPACER + rt.SMART_HANDLE + "AC"


def bead_primer() -> list[Segment]:
    return [seg("bead linker", LINKER), seg("SMART handle", rt.SMART_HANDLE, "tso"),
            seg("AC", "AC"), seg("cell barcode", "B" * 12, "cbc", placeholder=True),
            seg("UMI", "U" * 8, "umi", placeholder=True), seg("dT30", "T" * 30)]


def tso() -> list[Segment]:
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GAAT", "GAAT"),
            seg("rGrGrG", "GGG", "tso")]


def p5_hybrid() -> list[Segment]:
    return [seg("P5", il.P5, "p5"), seg("spacer", SPACER, "r1"),
            seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("AC", "AC")]


def capture_scene() -> Scene:
    mrna=[seg("RNA body","X"*28,placeholder=True),seg("poly(A)","A"*30)]
    sc=Scene(); sc.strand("mRNA",mrna,label="captured mRNA")
    sc.anneal("bead primer",bead_primer(),to="mRNA",pair=("dT30","poly(A)"),label="bead primer",mod5="bead")
    sc.arrow("bead primer","reverse transcriptase"); return sc


def first_strand() -> Construct:
    return Construct([*bead_primer(),seg("antisense cDNA","X"*28,placeholder=True),
                      seg("CCC","CCC","tso")],name="Seq-Well first-strand cDNA")


def template_switch_scene() -> Scene:
    first=first_strand(); sc=Scene(); sc.strand("first strand",list(first),label="first strand",mod5="bead")
    sc.anneal("TSO",tso(),to="first strand",pair=("rGrGrG","CCC"),label="TSO",above=True)
    sc.arrow("first strand","RT switches template"); return sc


def wta_product() -> Construct:
    return Construct([*first_strand(),seg("copied GAAT",revcomp("GAAT")),
                      seg("copied SMART handle",revcomp(rt.SMART_HANDLE),"tso")],
                     name="Seq-Well template-switched cDNA")


def selected_tagmented_strand(insert_nt: int=34) -> Construct:
    return Construct([
        seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("AC","AC"),
        seg("cell barcode","B"*12,"cbc",placeholder=True),seg("UMI","U"*8,"umi",placeholder=True),
        seg("poly(T)","T"*30),seg("antisense cDNA","X"*insert_nt,placeholder=True),
        seg("mosaic end reverse complement",nx.ME_RC,"me"),seg("s7 reverse complement",nx.S7_RC,"s7")],
        name="Seq-Well bead-end tagmented strand")


READ1=sp.custom("Read 1","Seq-Well custom Read 1",CUSTOM_R1,
                "Gierahn et al. 2017, Supplementary Table 1")
SEQ_PRIMERS=(READ1,sp.NEXTERA["I1"],sp.NEXTERA["R2"])


def final_library(insert_nt: int=34) -> Construct:
    selected=selected_tagmented_strand(insert_nt)
    lib=Construct([seg("P5",il.P5,"p5"),seg("spacer",SPACER,"r1"),*list(selected),
                   seg("i7 reverse complement","I"*8,"cbc",placeholder=True),
                   seg("P7 reverse complement",il.P7_RC,"p7")],name="Seq-Well sequencing library")
    for p in SEQ_PRIMERS:
        if sp.locate(lib,p) is None: raise ValueError(f"{p.name} does not land on {lib.name}")
    return lib


def primer_landings(lib: Construct|None=None):
    lib=lib or final_library(); out=[]
    for p in SEQ_PRIMERS:
        h=sp.locate(lib,p)
        if h is None: raise ValueError(f"{p.name} does not land on {lib.name}")
        out.append((p,h))
    return out


def read1_layout() -> list[tuple[str,str]]:
    lib=final_library(); bc=lib.get("cell barcode"); umi=lib.get("UMI")
    return [(f"1&ndash;{len(bc)}","cell barcode"),(f"{len(bc)+1}&ndash;{len(bc)+len(umi)}","UMI")]
