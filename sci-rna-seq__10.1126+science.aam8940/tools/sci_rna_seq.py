"""Molecular construct model for sci-RNA-seq (Cao et al., 2017)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

UMI_FEATURE = feature("umi", "umi", "random")
RT_BARCODE = feature("cell_barcode_rt", "cell_barcode", "combinatorial",
                     group="cell_barcode", part="RT round")
I5_FEATURE = feature("sample_index_i5", "sample_index", "unknown")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")
FEATURES = {"UMI": UMI_FEATURE, "RT barcode": RT_BARCODE, "i5": I5_FEATURE,
            "i7": I7_FEATURE, "i7 reverse complement": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


HANDLE = il.TRUSEQ_READ1[-18:]
RT_TEMPLATE = HANDLE + "N" * 8 + "N" * 10 + "T" * 30 + "VN"
P5_TEMPLATE = il.P5 + "N" * 10 + il.TRUSEQ_READ1
P7_TEMPLATE = il.P7 + "N" * 10 + nx.S7


def rt_primer() -> list[Segment]:
    return [seg("Read 1 handle", HANDLE, "r1"),
            seg("UMI", "U" * 8, "umi", placeholder=True),
            seg("RT barcode", "B" * 10, "cbc", placeholder=True),
            seg("dT30", "T" * 30), seg("V", "V", placeholder=True),
            seg("N", "N", placeholder=True)]


def p5_primer() -> list[Segment]:
    return [seg("P5", il.P5, "p5"), seg("i5", "I" * 10, "cbc", placeholder=True),
            seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1")]


def p7_primer() -> list[Segment]:
    return [seg("P7", il.P7, "p7"), seg("i7", "J" * 10, "cbc", placeholder=True),
            seg("s7", nx.S7, "s7")]


def rt_scene() -> Scene:
    mrna = [seg("RNA body", "X" * 28, placeholder=True), seg("poly(A)", "A" * 30)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA")
    sc.anneal("RT primer", rt_primer(), to="mRNA", pair=("dT30", "poly(A)"),
              label="barcoded RT primer", unpaired=("V", "N"))
    sc.arrow("RT primer", "reverse transcriptase")
    return sc


def first_strand() -> Construct:
    return Construct([*rt_primer(), seg("antisense cDNA", "X" * 34, placeholder=True)],
                     name="sci-RNA-seq first strand")


def second_strand_scene() -> Scene:
    first = first_strand()
    sc = Scene.duplex(list(first), label="cDNA")
    sc.labels("top")
    return sc


def selected_tagmented_strand(insert_nt: int = 34) -> Construct:
    """RT-primer end through the nearest s7 insertion after gap fill."""
    return Construct([
        seg("Read 1 handle", HANDLE, "r1"),
        seg("UMI", "U" * 8, "umi", placeholder=True),
        seg("RT barcode", "B" * 10, "cbc", placeholder=True),
        seg("poly(T)", "T" * 30), seg("V", "V", placeholder=True),
        seg("antisense cDNA", "X" * insert_nt, placeholder=True),
        seg("mosaic end reverse complement", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
    ], name="selected sci-RNA-seq tagmented strand")


SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.NEXTERA["I1"],
               sp.TRUSEQ["I2"], sp.NEXTERA["R2"])


def final_library(insert_nt: int = 34) -> Construct:
    selected = selected_tagmented_strand(insert_nt)
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("i5", "I" * 10, "cbc", placeholder=True),
        seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
        *list(selected)[1:],
        seg("i7 reverse complement", "J" * 10, "cbc", placeholder=True),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="sci-RNA-seq sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid final sci-RNA-seq construct: " + "; ".join(problems))
    return lib


def primer_landings(lib: Construct | None = None):
    lib = lib or final_library()
    out = []
    for primer in SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        if hit is None:
            raise ValueError(f"{primer.name} does not land on {lib.name}")
        out.append((primer, hit))
    return out


def read1_layout() -> list[tuple[str, str]]:
    return [(x.cycles, x.feature.label) for x in
            sp.feature_spans(final_library(), SEQ_PRIMERS, {"Read 1": 18})]


def protocol_graph() -> ProtocolGraph:
    source="Defining source: https://doi.org/10.1126/science.aam8940"
    graph=ProtocolGraph("sci-RNA-seq",authority=source)
    cells=graph.input("fixed_cells","fixed permeabilized cells",MaterialKind.BIOLOGICAL,
                      primary=True)
    primer=graph.input("rt_primer","barcoded anchored oligo-dT RT primer",
                       MaterialKind.OLIGO,state=tuple(rt_primer()))
    rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    first=graph.product("first_strand",first_strand().name,MaterialKind.MOLECULAR,
                        state=first_strand())
    graph.react("reverse_transcribe","Install the RT-well barcode and UMI in fixed cells",
                StepKind.REVERSE_TRANSCRIPTION,
                ((cells,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(first,))
    sorter=graph.input("facs","FACS well assignment",MaterialKind.SELECTION)
    sorted_cdna=graph.product("sorted_cdna","first-strand cDNA retained in sorted cells",
                              MaterialKind.MOLECULAR,state=first_strand())
    rejected=graph.product("sorting_waste","cells not retained for library construction",
                           MaterialKind.WASTE)
    graph.react("sort","Pool cells and sort them into PCR wells",StepKind.SAMPLING,
                ((first,InputRole.TEMPLATE),(sorter,InputRole.SELECTION)),
                (sorted_cdna,rejected),changes_molecule=False)
    graph.discard(rejected,"cells outside the selected PCR wells do not enter the library")
    random_primer=graph.input("second_strand_primer","random second-strand primer",
                              MaterialKind.OLIGO)
    second_pol=graph.input("second_strand_polymerase","second-strand polymerase",
                           MaterialKind.CATALYST)
    duplex=graph.product("double_stranded_cdna","double-stranded barcoded cDNA",
                         MaterialKind.MOLECULAR,state=second_strand_scene())
    graph.react("second_strand","Synthesize the second cDNA strand",StepKind.EXTENSION,
                ((sorted_cdna,InputRole.TEMPLATE),(random_primer,InputRole.STRUCTURAL),
                 (second_pol,InputRole.CATALYST)),(duplex,))
    carrier=graph.input("carrier_dna","carrier DNA",MaterialKind.MOLECULAR)
    transposome=graph.input("nextera_transposome","Nextera transposome",MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("selected_tagmented",selected_tagmented_strand().name,
                         MaterialKind.MOLECULAR,state=selected_tagmented_strand())
    tag_waste=graph.product("tagment_waste","carrier and fragments lacking the RT-handle/s7 pair",
                            MaterialKind.WASTE)
    graph.react("tagment","Tagment cDNA and retain the RT-primer end nearest an s7 insertion",
                StepKind.TAGMENTATION,
                ((duplex,InputRole.TEMPLATE),(carrier,InputRole.SPIKE_IN),
                 (transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),
                (tagged,tag_waste))
    graph.discard(tag_waste,"these products lack the two sites selected by library PCR")
    p5=graph.input("p5_primer","indexed P5 PCR primer",MaterialKind.OLIGO,
                   state=tuple(p5_primer()))
    p7=graph.input("p7_primer","indexed P7 PCR primer",MaterialKind.OLIGO,
                   state=tuple(p7_primer()))
    polymerase=graph.input("pcr_polymerase","PCR polymerase",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final)
    geometry=PCRGeometry.exact(
        selected_tagmented_strand(),
        PCRPrimer("indexed P5 PCR primer",P5_TEMPLATE,anneal_nt=len(HANDLE)),
        PCRPrimer("indexed P7 PCR primer",P7_TEMPLATE,anneal_nt=len(nx.S7)),
        statement="The defining supplement prints both indexed PCR primers and their RT-handle and s7 annealing sites.")
    graph.react("library_pcr","Select the barcode-bearing fragment and add PCR-well indexes",
                StepKind.PCR,
                ((tagged,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
                 (p7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=geometry)
    lengths={"Read 1":18,"Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":52}
    graph.sequence(library,SequencingPlan(
        library=final,primers=SEQ_PRIMERS,acquired_roles=tuple(lengths),
        read_lengths=lengths,
        read_uses={"Read 1":"extract","Index 1 (i7)":"demultiplex",
                   "Index 2 (i5)":"demultiplex","Read 2":"include"}))
    graph.declare_reviewed_scope((
        "fixed_cells","rt_primer","reverse_transcriptase","reverse_transcribe",
        "facs","sort","second_strand_primer","second_strand_polymerase","second_strand",
        "carrier_dna","nextera_transposome","tn5","tagment","p5_primer","p7_primer",
        "pcr_polymerase","library_pcr","library"),source)
    return graph.require_publishable()
