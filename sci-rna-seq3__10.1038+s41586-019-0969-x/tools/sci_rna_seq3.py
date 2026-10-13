"""Molecular construct model for sci-RNA-seq3 (Cao et al., 2019)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from topology import HairpinTopology

UMI_FEATURE = feature("umi", "umi", "random")
LIGATION_BARCODE = feature("cell_barcode_ligation", "cell_barcode", "combinatorial",
                           group="cell_barcode", part="ligation round")
RT_BARCODE = feature("cell_barcode_rt", "cell_barcode", "combinatorial",
                     group="cell_barcode", part="RT round")
I5_FEATURE = feature("sample_index_i5", "sample_index", "fixed")
I7_FEATURE = feature("sample_index_i7", "sample_index", "fixed")
FEATURES = {"UMI": UMI_FEATURE, "ligation barcode": LIGATION_BARCODE,
            "ligation barcode, read orientation": LIGATION_BARCODE,
            "barcode reverse complement": LIGATION_BARCODE, "RT barcode": RT_BARCODE,
            "i5": I5_FEATURE, "i7": I7_FEATURE,
            "i7 reverse complement": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


HANDLE = il.TRUSEQ_READ1[-18:]
RT_BARCODE_1 = "TCCTACCAGT"
HAIRPIN_BARCODE_1 = "ACAATCAAGT"
I5_1 = "CTCCATCGAG"
I7_1 = "CCGAATCCGA"
RT1 = "CAGAGC" + "N" * 8 + RT_BARCODE_1 + "T" * 30
HAIRPIN1 = "GCTCTG" + HAIRPIN_BARCODE_1 + "U" + HANDLE + revcomp(HAIRPIN_BARCODE_1)
HAIRPIN_TOPOLOGY = HairpinTopology(
    "sci-RNA-seq3 barcoded dU hairpin", HAIRPIN1,
    first_stem_start=len("GCTCTG"),
    second_stem_start=len(HAIRPIN1) - len(HAIRPIN_BARCODE_1),
    stem_length=len(HAIRPIN_BARCODE_1), colour="var(--c-cbc)")
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
    sc.ligation_end("first strand","ligation site",kind="splinted",polarity="5-prime",
                    owner="substrate",text="5′-phosphate at the ligation nick")
    sc.ligation_end("hairpin","barcode reverse complement",kind="splinted",
                    polarity="3-prime",owner="adapter",text="hairpin 3′-OH at the ligation nick")
    sc.mark("hairpin","6-nt overhang","splints the ligation nick"); return sc


def ligated_scene() -> Scene:
    product = ligated_product()
    sc = Scene(); sc.strand("product", list(product), label="ligated strand")
    sc.junction("product", "barcode reverse complement", "ligation site")
    return sc


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
READ_LENGTHS={"Read 1":34,"Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":52}
READ_USES={"Read 1":"extract"}


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
    spans=sp.feature_spans(final_library(),SEQ_PRIMERS,{"Read 1":34})
    by_id={x.feature.id:(x.cycles,x.feature.label) for x in spans}
    return [by_id["cell_barcode_ligation"],("11–16","CAGAGC"),
            by_id["umi"],by_id["cell_barcode_rt"]]


def protocol_graph() -> ProtocolGraph:
    source="Defining source: https://doi.org/10.1038/s41586-019-0969-x"
    graph=ProtocolGraph("sci-RNA-seq3",authority=source)
    nuclei=graph.input("fixed_nuclei","fixed permeabilized nuclei",MaterialKind.BIOLOGICAL,
                       primary=True)
    rt_oligo=graph.input("rt_primer","5-prime-phosphorylated barcoded RT primer",
                         MaterialKind.OLIGO,state=tuple(rt_primer()))
    rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    first=graph.product("first_strand",first_strand().name,MaterialKind.MOLECULAR,
                        state=first_strand())
    graph.react("reverse_transcribe","Install the RT barcode and UMI in situ",
                StepKind.REVERSE_TRANSCRIPTION,
                ((nuclei,InputRole.TEMPLATE),(rt_oligo,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(first,))
    hairpin=graph.input("hairpin_adapter","barcoded dU hairpin adapter",
                        MaterialKind.ADAPTER,state=tuple(hairpin_oligo()))
    ligase=graph.input("dna_ligase","DNA ligase",MaterialKind.CATALYST)
    ligated=graph.product("ligated_cdna",ligated_product().name,MaterialKind.MOLECULAR,
                          state=ligated_product())
    ligation_geometry=AdapterLigationGeometry.known(
        "splinted",ligation_scene(),
        statement="The defining supplement prints the complementary six-base splint, phosphorylated RT-primer end and hairpin adapter.")
    graph.react("hairpin_ligation","Add the second well barcode through the hairpin splint",
                StepKind.LIGATION,
                ((first,InputRole.TEMPLATE),(hairpin,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(ligated,),geometry=ligation_geometry)
    second_primer=graph.input("second_strand_primer","random second-strand primer",
                              MaterialKind.OLIGO)
    second_pol=graph.input("second_strand_polymerase","second-strand polymerase",
                           MaterialKind.CATALYST)
    duplex=graph.product("double_stranded_cdna","double-stranded hairpin-ligated cDNA",
                         MaterialKind.MOLECULAR,state=ligated_product())
    graph.react("second_strand","Synthesize the second cDNA strand",StepKind.EXTENSION,
                ((ligated,InputRole.TEMPLATE),(second_primer,InputRole.STRUCTURAL),
                 (second_pol,InputRole.CATALYST)),(duplex,))
    transposome=graph.input("n7_transposome","N7-only transposome",MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("tagmented_cdna","N7-tagmented hairpin cDNA",MaterialKind.MOLECULAR)
    tag_waste=graph.product("tagment_waste","tagmentation products lacking the selected ends",
                            MaterialKind.WASTE)
    graph.react("tagment","Tagment double-stranded cDNA with N7-only Tn5",
                StepKind.TAGMENTATION,
                ((duplex,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,tag_waste))
    graph.discard(tag_waste,"fragments without the barcode-bearing hairpin end are not PCR-selected")
    user=graph.input("user_enzyme","USER enzyme",MaterialKind.CATALYST)
    opened=graph.product("opened_cdna",user_cleaved_tagmented().name,
                         MaterialKind.MOLECULAR,state=user_cleaved_tagmented())
    released=graph.product("released_hairpin_arm","released dU-containing hairpin arm",
                           MaterialKind.WASTE)
    graph.react("user_digest","Open the hairpin at dU with USER",StepKind.DIGESTION,
                ((tagged,InputRole.TEMPLATE),(user,InputRole.CATALYST)),(opened,released))
    graph.discard(released,"the cleaved outer hairpin arm is not part of the final library")
    p5=graph.input("p5_primer","indexed P5 PCR primer",MaterialKind.OLIGO,
                   state=tuple(p5_primer()))
    p7=graph.input("p7_primer","indexed P7 PCR primer",MaterialKind.OLIGO,
                   state=tuple(p7_primer()))
    polymerase=graph.input("pcr_polymerase","PCR polymerase",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final)
    geometry=PCRGeometry.exact(
        user_cleaved_tagmented(),
        PCRPrimer("indexed P5 PCR primer",P5_1,anneal_nt=len(HANDLE)),
        PCRPrimer("indexed P7 PCR primer",P7_1,anneal_nt=len(nx.S7)),
        statement="The defining supplement prints both indexed PCR primers and their Read-1-handle and s7 annealing sites.")
    graph.react("library_pcr","Add the PCR-well i5/i7 pair",StepKind.PCR,
                ((opened,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
                 (p7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=geometry)
    graph.sequence(library,SequencingPlan(
        library=final,primers=SEQ_PRIMERS,acquired_roles=tuple(READ_LENGTHS),
        read_lengths=READ_LENGTHS,
        read_uses={"Read 1":"extract","Index 1 (i7)":"demultiplex",
                   "Index 2 (i5)":"demultiplex","Read 2":"include"}))
    graph.declare_reviewed_scope((
        "fixed_nuclei","rt_primer","reverse_transcriptase","reverse_transcribe",
        "hairpin_adapter","dna_ligase","hairpin_ligation","second_strand_primer",
        "second_strand_polymerase","second_strand","n7_transposome","tn5","tagment",
        "user_enzyme","user_digest","p5_primer","p7_primer","pcr_polymerase",
        "library_pcr","library"),source)
    return graph.require_publishable()
