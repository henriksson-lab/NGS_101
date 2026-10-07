"""Molecular construct model for sci-RNA-seq3 (Cao et al., 2019)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, revcomp


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


HANDLE = il.TRUSEQ_READ1[-18:]
RT_BARCODE_1 = "TCCTACCAGT"
HAIRPIN_BARCODE_1 = "ACAATCAAGT"
I5_1 = "CTCCATCGAG"
I7_1 = "CCGAATCCGA"
RT1 = "CAGAGC" + "N" * 8 + RT_BARCODE_1 + "T" * 30
HAIRPIN1 = "GCTCTG" + HAIRPIN_BARCODE_1 + "U" + HANDLE + revcomp(HAIRPIN_BARCODE_1)
P5_1 = il.P5 + I5_1 + il.TRUSEQ_READ1
P7_1 = il.P7 + I7_1 + nx.S7


def rt_primer() -> list[Segment]:
    return [seg("ligation site", "CAGAGC"), seg("UMI", "U" * 8, "umi", placeholder=True),
            seg("RT barcode", RT_BARCODE_1, "cbc"), seg("dT30", "T" * 30)]


def hairpin_oligo() -> list[Segment]:
    """Representative adaptor; the final arm is derived, never transcribed twice."""
    return [seg("6-nt overhang", "GCTCTG"),
            seg("ligation barcode", HAIRPIN_BARCODE_1, "cbc"),
            seg("dU", "U", "w1", placeholder=True), seg("Read 1 handle", HANDLE, "r1"),
            seg("barcode reverse complement", revcomp(HAIRPIN_BARCODE_1), "cbc")]


def p5_primer() -> list[Segment]:
    return [seg("P5", il.P5, "p5"), seg("i5", I5_1, "cbc"),
            seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1")]


def p7_primer() -> list[Segment]:
    return [seg("P7", il.P7, "p7"), seg("i7", I7_1, "cbc"), seg("s7", nx.S7, "s7")]


def rt_scene() -> Scene:
    mrna=[seg("RNA body","X"*28,placeholder=True),seg("poly(A)","A"*30)]
    sc=Scene(); sc.strand("mRNA",mrna,label="mRNA")
    sc.anneal("RT primer",rt_primer(),to="mRNA",pair=("dT30","poly(A)"),label="barcoded RT primer",mod5="p")
    sc.arrow("RT primer","reverse transcriptase"); return sc


def first_strand() -> Construct:
    return Construct([*rt_primer(),seg("antisense cDNA","X"*34,placeholder=True)],name="sci-RNA-seq3 first strand")


def ligation_scene() -> Scene:
    target=first_strand(); hp=hairpin_oligo()
    sc=Scene(); sc.strand("first strand",list(target),label="first strand",mod5="p")
    sc.anneal("hairpin",hp,to="first strand",pair=("6-nt overhang","ligation site"),
              label="hairpin adaptor",unpaired=("ligation barcode","dU","Read 1 handle","barcode reverse complement"))
    sc.mark("hairpin","6-nt overhang","splints the ligation nick"); return sc


def ligated_product() -> Construct:
    return Construct([*hairpin_oligo(),*list(first_strand())],name="hairpin-ligated sci-RNA-seq3 cDNA")


def user_cleaved_tagmented(insert_nt: int=34) -> Construct:
    """PCR-bearing strand after N7-only tagmentation and USER opening."""
    return Construct([
        seg("Read 1 handle",HANDLE,"r1"),
        seg("ligation barcode, read orientation",revcomp(HAIRPIN_BARCODE_1),"cbc"),
        seg("ligation site","CAGAGC"),seg("UMI","U"*8,"umi",placeholder=True),
        seg("RT barcode",RT_BARCODE_1,"cbc"),seg("poly(T)","T"*30),
        seg("antisense cDNA","X"*insert_nt,placeholder=True),
        seg("mosaic end reverse complement",nx.ME_RC,"me"),
        seg("s7 reverse complement",nx.S7_RC,"s7")],
        name="USER-opened, N7-tagmented sci-RNA-seq3 strand")


SEQ_PRIMERS=(sp.TRUSEQ["R1"],sp.NEXTERA["I1"],sp.TRUSEQ["I2"],sp.NEXTERA["R2"])


def final_library(insert_nt: int=34) -> Construct:
    pre=user_cleaved_tagmented(insert_nt)
    lib=Construct([
        seg("P5",il.P5,"p5"),seg("i5",I5_1,"cbc"),seg("TruSeq Read 1",il.TRUSEQ_READ1,"r1"),
        *list(pre)[1:],seg("i7 reverse complement",revcomp(I7_1),"cbc"),
        seg("P7 reverse complement",il.P7_RC,"p7")],name="sci-RNA-seq3 sequencing library")
    problems=sp.verify(lib,SEQ_PRIMERS)
    if problems: raise ValueError("invalid final sci-RNA-seq3 construct: "+"; ".join(problems))
    return lib


def primer_landings(lib: Construct|None=None):
    lib=lib or final_library(); out=[]
    for p in SEQ_PRIMERS:
        h=sp.locate(lib,p)
        if h is None: raise ValueError(f"{p.name} does not land on {lib.name}")
        out.append((p,h))
    return out


def read1_layout() -> list[tuple[str,str]]:
    lib=final_library(); names=[("ligation barcode, read orientation","ligation barcode"),("ligation site","CAGAGC"),("UMI","UMI"),("RT barcode","RT barcode")]
    rows=[]; cycle=1
    for name,label in names:
        end=cycle+len(lib.get(name))-1; rows.append((f"{cycle}&ndash;{end}",label)); cycle=end+1
    return rows
