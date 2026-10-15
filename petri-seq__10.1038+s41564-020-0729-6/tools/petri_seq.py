"""Molecular construct model for PETRI-seq (Blattman et al., 2020)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (AnnotationSlot, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)

UMI_FEATURE = feature("umi", "umi", "random")
I5_FEATURE = feature("sample_index_i5", "sample_index", "fixed")
I7_FEATURE = feature("sample_index_i7", "sample_index", "fixed")


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    if name == "UMI": kw.setdefault("feature", UMI_FEATURE)
    elif name == "i5 index": kw.setdefault("feature", I5_FEATURE)
    elif name == "i7 index read": kw.setdefault("feature", I7_FEATURE)
    return Segment(name=name, top=top, tag=tag, **kw)


R1_PREFIX = "GCCAGA"
R2_PREFIX = "GCTTCGC"
R2_SUFFIX = "CCTCCTAC"
R3_PREFIX = "AGAATACACGACGCTCTTCCGATCT"
R3_SUFFIX = "GGTCCTTG"
LINKER_R2 = "STCTGGCGTAGGAGGW"
LINKER_R3 = "GCGAAGCCAAGGACCW"
# Members 1 and 2 of the three 96-well oligo sets in Supplementary Table 3.
R1_MEMBERS = tuple(R1_PREFIX + bc + "N" * 6 for bc in ("CAGAGAA", "CCAGGAA"))
R2_MEMBERS = tuple(R2_PREFIX + bc + R2_SUFFIX for bc in ("GACCTTA", "AGAGTGA"))
R3_MEMBERS = tuple(R3_PREFIX + "N" * 7 + bc + R3_SUFFIX
                   for bc in ("GTGTGAA", "TTGGTGA"))
I5_OLIGO_INDEX = il.NEBNEXT_I5_SET1["i501"]
I7_OLIGO_INDEX = "TAAGGCGA"
I7_INDEX_READ = "TCGCCTTA"


def barcode(name: str) -> Segment:
    round_name = name.replace("barcode ", "round ")
    return seg(name, "B" * 7, "cbc", placeholder=True,
               feature=feature(name.replace(" ", "_"), "cell_barcode", "combinatorial",
                               group="cell_barcode", part=round_name))


def round1_primer() -> list[Segment]:
    return [seg("round-1 constant", R1_PREFIX), barcode("barcode 1"),
            seg("random hexamer", "N" * 6, placeholder=True)]


def round2_oligo() -> list[Segment]:
    return [seg("round-2 constant", R2_PREFIX), barcode("barcode 2"),
            seg("round-2 junction", R2_SUFFIX)]


def round3_oligo() -> list[Segment]:
    return [seg("round-3 library handle", R3_PREFIX, "r1"),
            seg("UMI", "N" * 7, "umi", placeholder=True), barcode("barcode 3"),
            seg("round-3 junction", R3_SUFFIX)]


def after_round1(insert_nt: int = 28) -> Construct:
    return Construct([*round1_primer(), seg("cDNA", "X" * insert_nt, placeholder=True)],
                     name="PETRI-seq after RT barcode")


def after_round2(insert_nt: int = 28) -> Construct:
    return Construct([*round2_oligo(), *after_round1(insert_nt)],
                     name="PETRI-seq after second barcode")


def after_round3(insert_nt: int = 28) -> Construct:
    return Construct([*round3_oligo(), *after_round2(insert_nt)],
                     name="PETRI-seq after third barcode")


def splint_geometry(round_no: int) -> AdapterLigationGeometry:
    linker = LINKER_R2 if round_no == 2 else LINKER_R3
    old_head = R1_PREFIX + "C" if round_no == 2 else R2_PREFIX
    new_tail = "A" + (R2_SUFFIX if round_no == 2 else R3_SUFFIX)
    old_name = f"round-{round_no - 1} phosphorylated head"
    new_name = f"round-{round_no} incoming 3-prime end"
    sc = Scene()
    sc.strand("splint", [seg("old-side splint", linker[:len(old_head)], placeholder=True),
                         seg("new-side splint", linker[len(old_head):], placeholder=True)],
              label=f"round-{round_no} linker")
    sc.anneal("old barcode oligo", [seg(old_name, old_head)],
              to="splint", pair=(old_name, "old-side splint"),
              label="retained 5′-phosphorylated oligo", mod5="p")
    sc.anneal("incoming barcode oligo", [seg(new_name, new_tail)],
              to="splint", pair=(new_name, "new-side splint"),
              label="incoming 3′-OH oligo")
    sc.ligation_end("old barcode oligo", old_name,
                    kind="splinted", polarity="5-prime", owner="adapter",
                    text="5′-phosphorylated nick partner; linker pairs to both source oligos")
    return AdapterLigationGeometry.known(
        "splinted", sc,
        statement=f"The PETRI-seq supplement prints the 16-base round-{round_no} linker; its two eight-base halves align the barcode oligos across the nick.")


def workflow() -> Workflow:
    r1, r2, r3 = after_round1(), after_round2(), after_round3()
    wf = Workflow(MolecularState("Round-1 barcoded cDNA", tuple(Scene.duplex(list(r1), label="after RT").rows())))
    wf.ligate_adapter("Ligate barcode round 2", splint_geometry(2),
                      Scene.duplex(list(r2), label="after round 2").rows(),
                      note="The splint spans the round-2/round-1 nick.")
    wf.ligate_adapter("Ligate UMI and barcode round 3", splint_geometry(3),
                      Scene.duplex(list(r3), label="after round 3").rows(),
                      note="The third oligo supplies the UMI, barcode 3 and future Read-1 handle.")
    wf.react("Second-strand synthesis, tagmentation and selective PCR",
             final_scene().rows(), note="NEB i50x and Nextera N7xx select the two required ends.")
    return wf


SEQ_PRIMERS = tuple([sp.TRUSEQ["R1"], sp.NEXTERA["I1"],
                     sp.TRUSEQ["I2"], sp.NEXTERA["R2"]])


def final_library(insert_nt: int = 28) -> Construct:
    """Selected round-3-handle/s7 product after mixed NEB-i5/Nextera-i7 PCR."""
    lib = Construct([
        seg("P5", il.P5, "p5"), seg("i5 index", I5_OLIGO_INDEX, "cbc"),
        seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
        seg("UMI", "N" * 7, "umi", placeholder=True), barcode("barcode 3"),
        seg("round-3 junction", R3_SUFFIX), seg("round-2 constant", R2_PREFIX),
        barcode("barcode 2"), seg("round-2 junction", R2_SUFFIX),
        seg("round-1 constant", R1_PREFIX), barcode("barcode 1"),
        seg("random-hexamer-derived bases", "N" * 6, placeholder=True),
        seg("cDNA", "X" * insert_nt, placeholder=True),
        seg("mosaic end reverse complement", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 index read", I7_INDEX_READ, "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="PETRI-seq sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid final PETRI-seq construct: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")


def read1_layout() -> list[tuple[str, str]]:
    spans=sp.feature_spans(final_library(),SEQ_PRIMERS,{"Read 1":58})
    by_start={x.cycle_start:(x.cycles,x.feature.label+(
              f" ({x.feature.part})" if x.feature.part else "")) for x in spans}
    return [by_start[1],by_start[8], ("15–29", "round-3/round-2 junction"),
            by_start[30], ("37–50", "round-2/round-1 junction"), by_start[51],
            ("58", "first random-hexamer-derived base")]


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1038/s41564-020-0729-6"
    graph = ProtocolGraph("PETRI-seq", authority=source)
    cells = graph.input("cells", "fixed permeabilized bacterial cells",
                        MaterialKind.BIOLOGICAL, primary=True)
    r1 = graph.input("round1_primers", "round-1 barcoded random-hexamer RT primers",
                     MaterialKind.OLIGO, state=tuple(R1_MEMBERS), source=source + "; Supplementary Table 3")
    rt = graph.input("reverse_transcriptase", "reverse transcriptase",
                     MaterialKind.CATALYST)
    well1 = graph.product("round1_well_a", "representative round-1 well A cells",
                          MaterialKind.BIOLOGICAL)
    well2 = graph.product("round1_well_b", "representative round-1 well B cells",
                          MaterialKind.BIOLOGICAL)
    graph.react("round1_split", "Split fixed cells across round-1 RT wells",
                StepKind.SPLIT, ((cells, InputRole.TEMPLATE),), (well1, well2))
    cdna1 = graph.product("round1_cdna_a", "round-1-barcoded cDNA from well A",
                          MaterialKind.MOLECULAR, state=after_round1())
    cdna2 = graph.product("round1_cdna_b", "round-1-barcoded cDNA from well B",
                          MaterialKind.MOLECULAR, state=after_round1())
    graph.react("round1_rt_a", "Reverse-transcribe with round-1 barcode set A",
                StepKind.REVERSE_TRANSCRIPTION,
                ((well1, InputRole.TEMPLATE), (r1, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna1,))
    graph.react("round1_rt_b", "Reverse-transcribe with round-1 barcode set B",
                StepKind.REVERSE_TRANSCRIPTION,
                ((well2, InputRole.TEMPLATE), (r1, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna2,))
    pooled1 = graph.product("round1_pool", "pooled round-1-barcoded cells",
                            MaterialKind.MOLECULAR)
    graph.react("round1_pooling", "Pool round-1 wells", StepKind.POOL,
                ((cdna1, InputRole.TEMPLATE), (cdna2, InputRole.TEMPLATE)), (pooled1,))

    r2 = graph.input("round2_oligos", "round-2 barcode oligos",
                     MaterialKind.ADAPTER, state=tuple(R2_MEMBERS), source=source + "; Supplementary Table 3")
    linker2 = graph.input("round2_linker", "round-2 splint linker",
                          MaterialKind.OLIGO, state=LINKER_R2, source=source + "; Supplementary Table 2")
    ligase = graph.input("ligase", "DNA ligase", MaterialKind.CATALYST)
    r2a = graph.product("round2_well_a", "representative round-2 well A cells",
                        MaterialKind.MOLECULAR)
    r2b = graph.product("round2_well_b", "representative round-2 well B cells",
                        MaterialKind.MOLECULAR)
    graph.react("round2_split", "Redistribute cells across round-2 ligation wells",
                StepKind.SPLIT, ((pooled1, InputRole.TEMPLATE),), (r2a, r2b))
    lig2a = graph.product("round2_product_a", "round-2-ligated cells from well A",
                          MaterialKind.MOLECULAR, state=after_round2())
    lig2b = graph.product("round2_product_b", "round-2-ligated cells from well B",
                          MaterialKind.MOLECULAR, state=after_round2())
    for suffix, substrate, output in (("a", r2a, lig2a), ("b", r2b, lig2b)):
        graph.react("round2_ligation_" + suffix,
                    "Splint-ligate the round-2 barcode in representative well " + suffix.upper(),
                    StepKind.LIGATION,
                    ((substrate, InputRole.TEMPLATE), (r2, InputRole.STRUCTURAL),
                     (linker2, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)),
                    (output,), geometry=splint_geometry(2))
    pooled2 = graph.product("round2_pool", "pooled two-barcode cells",
                            MaterialKind.MOLECULAR)
    graph.react("round2_pooling", "Pool round-2 wells", StepKind.POOL,
                ((lig2a, InputRole.TEMPLATE), (lig2b, InputRole.TEMPLATE)), (pooled2,))

    r3 = graph.input("round3_oligos", "round-3 UMI/barcode oligos",
                     MaterialKind.ADAPTER, state=tuple(R3_MEMBERS), source=source + "; Supplementary Table 3")
    linker3 = graph.input("round3_linker", "round-3 splint linker",
                          MaterialKind.OLIGO, state=LINKER_R3, source=source + "; Supplementary Table 2")
    r3a = graph.product("round3_well_a", "representative round-3 well A cells",
                        MaterialKind.MOLECULAR)
    r3b = graph.product("round3_well_b", "representative round-3 well B cells",
                        MaterialKind.MOLECULAR)
    graph.react("round3_split", "Redistribute cells across round-3 ligation wells",
                StepKind.SPLIT, ((pooled2, InputRole.TEMPLATE),), (r3a, r3b))
    lig3a = graph.product("round3_product_a", "round-3-ligated cells from well A",
                          MaterialKind.MOLECULAR, state=after_round3())
    lig3b = graph.product("round3_product_b", "round-3-ligated cells from well B",
                          MaterialKind.MOLECULAR, state=after_round3())
    for suffix, substrate, output in (("a", r3a, lig3a), ("b", r3b, lig3b)):
        graph.react("round3_ligation_" + suffix,
                    "Splint-ligate the round-3 UMI/barcode in representative well " + suffix.upper(),
                    StepKind.LIGATION,
                    ((substrate, InputRole.TEMPLATE), (r3, InputRole.STRUCTURAL),
                     (linker3, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)),
                    (output,), geometry=splint_geometry(3))
    pooled3 = graph.product("round3_pool", "pooled three-barcode cDNA",
                            MaterialKind.MOLECULAR)
    graph.react("round3_pooling", "Pool round-3 wells", StepKind.POOL,
                ((lig3a, InputRole.TEMPLATE), (lig3b, InputRole.TEMPLATE)), (pooled3,))
    second_pol = graph.input("second_strand_polymerase", "second-strand polymerase",
                             MaterialKind.CATALYST)
    duplex = graph.product("double_stranded_cdna", "double-stranded barcoded cDNA",
                           MaterialKind.MOLECULAR)
    graph.react("second_strand", "Synthesize the second cDNA strand",
                StepKind.EXTENSION,
                ((pooled3, InputRole.TEMPLATE), (second_pol, InputRole.CATALYST)),
                (duplex,))
    transposome = graph.input("transposome", "Nextera transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagmented_cdna", "tagmented barcoded cDNA",
                           MaterialKind.MOLECULAR)
    waste = graph.product("tagmentation_waste", "non-library tagmentation products",
                          MaterialKind.WASTE)
    graph.react("tagment", "Tagment double-stranded cDNA", StepKind.TAGMENTATION,
                ((duplex, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "fragments lacking the round-3 handle and usable Tn5 end are not amplified")
    i5 = graph.input("i5_primer", "NEB i50x selective PCR primer", MaterialKind.OLIGO)
    i7 = graph.input("i7_primer", "Nextera N7xx selective PCR primer", MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    final = final_library()
    library = graph.product("library", final.name, MaterialKind.LIBRARY, state=final)
    graph.react("selective_pcr", "Select the round-3-handle/Tn5 fragments and complete adapters",
                StepKind.PCR,
                ((tagged, InputRole.TEMPLATE), (i5, InputRole.STRUCTURAL),
                 (i7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The tagmented pre-PCR fragment is variable; the final sequencing construct is exact."))
    graph.sequence(library, SequencingPlan(
        library=final, primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
        read_lengths={"Read 1": 58, "Index 1 (i7)": 8,
                      "Index 2 (i5)": 8, "Read 2": 17},
        read_uses={"Read 1": "extract", "Index 1 (i7)": "demultiplex",
                   "Index 2 (i5)": "demultiplex", "Read 2": "include"}))
    graph.declare_reviewed_scope((
        "cells", "round1_primers", "reverse_transcriptase", "round1_split",
        "round1_rt_a", "round1_rt_b", "round1_pooling", "round2_oligos",
        "round2_linker", "ligase", "round2_split", "round2_ligation_a",
        "round2_ligation_b", "round2_pooling", "round3_oligos", "round3_linker",
        "round3_split", "round3_ligation_a", "round3_ligation_b", "round3_pooling",
        "second_strand_polymerase", "second_strand", "transposome", "tn5",
        "tagment", "i5_primer", "i7_primer", "pcr_polymerase", "selective_pcr",
        "library"), source)
    graph.annotate(SourceDiscrepancy(
        target="round1_rt_a", title="SB5/SB6 qPCR-pair table labels",
        claims=(
            SourceClaim("Blattman et al. Supplementary Table 2, SB5/SB6 descriptions",
                        "Each oligo's pairing description names itself ('with SB5' / 'with SB6')."),
            SourceClaim("Blattman et al. Methods, qPCR quality-control pair",
                        "SB5 and SB6 are used together as the qPCR pair.")),
        model_choice="The schematic's sequenced construct does not include this quality-control PCR; SB5 and SB6 are treated as a pair in the source notes.",
        evidence="The Methods explicitly name the pair used for genomic DNA and random-hexamer cDNA checks.",
        uncertainty="The supplementary table has no formal corrected version in the sources reviewed.",
        slot=AnnotationSlot.INTRO))
    return graph.require_publishable()
