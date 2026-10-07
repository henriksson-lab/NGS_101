"""Molecular model for Drop-seq, with unavailable oligo-table regions inferred by default."""
from __future__ import annotations

import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, revcomp


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def inferred(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    """Unavailable Drop-seq Table S6 sequence: always rendered as inferred."""
    return Segment(name=name, top=top, tag=tag, inferred=True, **kw)


LINKER="T"*7
SPACER="GCCTGTCCGCGG"


def bead_primer(batch: str="B") -> list[Segment]:
    if batch not in ("A","B"): raise ValueError("Drop-seq bead batch must be A or B")
    constant="ACGT" if batch=="A" else "AC"
    return [inferred("bead linker",LINKER),inferred("SMART handle",rt.SMART_HANDLE,"tso"),
            inferred(f"batch {batch} constant",constant),
            inferred("cell barcode","B"*12,"cbc",placeholder=True),
            inferred("UMI","U"*8,"umi",placeholder=True),inferred("dT30","T"*30)]


def tso() -> list[Segment]:
    return [inferred("SMART handle",rt.SMART_HANDLE,"tso"),inferred("GAAT","GAAT"),
            inferred("rGrGrG","GGG","tso")]


def p5_hybrid(batch: str="B") -> list[Segment]:
    constant="ACGT" if batch=="A" else "AC"
    return [seg("P5",il.P5,"p5"),inferred("spacer",SPACER,"r1"),
            inferred("SMART handle",rt.SMART_HANDLE,"tso"),inferred(f"batch {batch} constant",constant)]


def custom_r1(batch: str="B") -> str:
    constant="ACGT" if batch=="A" else "AC"
    return SPACER+rt.SMART_HANDLE+constant


def capture_scene(batch: str="B") -> Scene:
    mrna=[seg("RNA body","X"*28,placeholder=True),seg("poly(A)","A"*30)]
    sc=Scene(); sc.strand("mRNA",mrna,label="captured mRNA")
    sc.anneal("bead primer",bead_primer(batch),to="mRNA",pair=("dT30","poly(A)"),label=f"batch {batch} bead",mod5="bead")
    sc.arrow("bead primer","reverse transcriptase"); return sc


def first_strand(batch: str="B") -> Construct:
    return Construct([*bead_primer(batch),seg("antisense cDNA","X"*28,placeholder=True),
                      inferred("CCC","CCC","tso")],name=f"Drop-seq batch {batch} first strand")


def template_switch_scene(batch: str="B") -> Scene:
    first=first_strand(batch); sc=Scene(); sc.strand("first strand",list(first),label="first strand",mod5="bead")
    sc.anneal("TSO",tso(),to="first strand",pair=("rGrGrG","CCC"),label="TSO",above=True)
    sc.arrow("first strand","RT switches template"); return sc


def selected_tagmented_strand(batch: str="B",insert_nt: int=34) -> Construct:
    constant="ACGT" if batch=="A" else "AC"
    return Construct([inferred("SMART handle",rt.SMART_HANDLE,"tso"),
        inferred(f"batch {batch} constant",constant),inferred("cell barcode","B"*12,"cbc",placeholder=True),
        inferred("UMI","U"*8,"umi",placeholder=True),inferred("poly(T)","T"*30),
        seg("antisense cDNA","X"*insert_nt,placeholder=True),seg("mosaic end reverse complement",nx.ME_RC,"me"),
        seg("s7 reverse complement",nx.S7_RC,"s7")],name=f"Drop-seq batch {batch} bead-end fragment")


def seq_primers(batch: str="B"):
    return (sp.custom("Read 1",f"Drop-seq batch {batch} custom Read 1",custom_r1(batch),
                      "secondary reconstruction; Drop-seq Table S6 unavailable"),
            sp.NEXTERA["I1"],sp.NEXTERA["R2"])


def final_library(batch: str="B",insert_nt: int=34) -> Construct:
    selected=selected_tagmented_strand(batch,insert_nt)
    lib=Construct([seg("P5",il.P5,"p5"),inferred("spacer",SPACER,"r1"),*list(selected),
                   seg("i7 reverse complement","I"*8,"cbc",placeholder=True),
                   seg("P7 reverse complement",il.P7_RC,"p7")],name=f"Drop-seq batch {batch} library")
    for p in seq_primers(batch):
        if sp.locate(lib,p) is None: raise ValueError(f"{p.name} does not land on {lib.name}")
    return lib


def primer_landings(batch: str="B"):
    lib=final_library(batch); out=[]
    for p in seq_primers(batch):
        h=sp.locate(lib,p)
        if h is None: raise ValueError(f"{p.name} does not land on {lib.name}")
        out.append((p,h))
    return out


def read1_layout(batch: str="B") -> list[tuple[str,str]]:
    lib=final_library(batch); bc=lib.get("cell barcode"); umi=lib.get("UMI")
    return [(f"1&ndash;{len(bc)}","cell barcode"),(f"{len(bc)+1}&ndash;{len(bc)+len(umi)}","UMI")]
