"""Molecular construct model for Microwell-seq (Han et al., 2018)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)

UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "fixed")


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    if name == "UMI":
        kw.setdefault("feature", UMI_FEATURE)
    elif name == "i7 index read":
        kw.setdefault("feature", I7_FEATURE)
    return Segment(name=name, top=top, tag=tag, **kw)


BEAD_HEAD = "TTTAGGGATAACAGGGTAAT"
BEAD_END = "ACGT"
LINK1 = "CGACTCACTACAGGG"
LINK2 = "TCGGTGACACGATCG"
P5_SPACER = "GCCTGTCCGCGG"
I7_OLIGO_INDEX = "TAAGGCGA"
I7_INDEX_READ = revcomp(I7_OLIGO_INDEX)
TSO = rt.SMART_HANDLE + "GAAT" + "GGG"
P5_PRIMER = il.P5 + P5_SPACER + rt.SMART_HANDLE + "AC"
CUSTOM_R1 = P5_SPACER + rt.SMART_HANDLE + BEAD_END


def barcode(name: str) -> Segment:
    part = name.replace("barcode ", "round ")
    return seg(name, "B" * 6, "cbc", placeholder=True,
               feature=feature(name.replace(" ", "_"), "cell_barcode", "combinatorial",
                               group="cell_barcode", part=part))


def bead_oligo() -> Construct:
    """Published seqA bead architecture after three extension rounds."""
    return Construct([
        seg("bead head", BEAD_HEAD), seg("SMART handle", rt.SMART_HANDLE, "tso"),
        seg("seqA constant", BEAD_END), barcode("barcode 1"), seg("linker 1", LINK1),
        barcode("barcode 2"), seg("linker 2", LINK2), barcode("barcode 3"),
        seg("UMI", "N" * 6, "umi", placeholder=True), seg("poly(T)", "T" * 30),
    ], name="Microwell-seq seqA bead oligo")


def tso_oligo() -> list[Segment]:
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GAAT", "GAAT"),
            seg("rGrG + modified G", "rGrG+G", placeholder=True)]


def rt_scene(insert_nt: int = 28) -> Scene:
    mrna = [seg("transcript", "X" * insert_nt, placeholder=True), seg("poly(A)", "A" * 30)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA", mod3="poly(A)")
    sc.anneal("bead", list(bead_oligo()), to="mRNA", pair=("poly(T)", "poly(A)"),
              label="barcoded bead primer")
    sc.arrow("bead", "reverse transcription and template switching")
    return sc


def inferred_s7_end() -> list[Segment]:
    """Same-lab/primer-supported, but not printed for the defining 2018 protocol."""
    return [seg("mosaic end reverse complement", nx.ME_RC, "me", inferred=True),
            seg("s7 reverse complement", nx.S7_RC, "s7", inferred=True),
            seg("i7 index read", I7_INDEX_READ, "cbc", inferred=True),
            seg("P7 reverse complement", il.P7_RC, "p7", inferred=True)]


SEQ_PRIMERS = (
    sp.custom("Read 1", "Microwell-seq custom Read 1", CUSTOM_R1,
              "Han et al. 2018, Supplementary Table S2"),
    sp.NEXTERA["I1"], sp.NEXTERA["R2"],
)
RUN_ROLES = ("Read 1", "Index 1 (i7)", "Read 2")


def final_library(insert_nt: int = 28) -> Construct:
    bead = list(bead_oligo())
    lib = Construct([
        seg("P5", il.P5, "p5"), seg("custom spacer", P5_SPACER),
        *bead[1:], seg("cDNA", "X" * insert_nt, placeholder=True), *inferred_s7_end(),
    ], name="Microwell-seq sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS, required_roles=RUN_ROLES)
    if problems:
        raise ValueError("invalid final Microwell-seq construct: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")


def read1_layout() -> list[tuple[str, str]]:
    spans = sp.feature_spans(final_library(), SEQ_PRIMERS, {"Read 1": 54})
    by_start = {x.cycle_start: (x.cycles, x.feature.label +
                (f" ({x.feature.part})" if x.feature.part else "")) for x in spans}
    return [by_start[1], ("7–21", "linker 1"), by_start[22],
            ("28–42", "linker 2"), by_start[43], by_start[49]]


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1016/j.cell.2018.02.001"
    inferred_source = source + "; Tn5 adapter inferred from printed PCR primers and later same-lab description"
    graph = ProtocolGraph("Microwell-seq", authority=source)
    bead_seed = graph.input("bead_seed", "seqA bead oligo seed", MaterialKind.MOLECULAR,
                            primary=True)
    barcode_oligos = graph.input("barcode_oligos", "three split-pool barcode-extension oligo sets",
                                 MaterialKind.OLIGO)
    extension_polymerase = graph.input("extension_polymerase", "bead-extension polymerase",
                                      MaterialKind.CATALYST)
    bead = graph.product("barcoded_bead", bead_oligo().name, MaterialKind.MOLECULAR,
                         state=bead_oligo())
    graph.react("bead_assembly", "Build the bead barcode by three split-pool extension rounds",
                StepKind.EXTENSION,
                ((bead_seed, InputRole.TEMPLATE), (barcode_oligos, InputRole.STRUCTURAL),
                 (extension_polymerase, InputRole.CATALYST)), (bead,))
    rna = graph.input("rna", "microwell cellular RNA", MaterialKind.BIOLOGICAL, primary=True)
    reverse_transcriptase = graph.input("reverse_transcriptase",
                                       "template-switching reverse transcriptase",
                                       MaterialKind.CATALYST)
    first_strand = graph.product("first_strand", "bead-primed first-strand cDNA",
                                 MaterialKind.MOLECULAR, state=rt_scene())
    graph.react("reverse_transcription", "Capture microwell RNA and reverse-transcribe",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (bead, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (first_strand,))
    tso = graph.input("tso", "published template-switching oligo", MaterialKind.OLIGO)
    switched = graph.product("switched_cdna", "SMART-handle-flanked cDNA",
                             MaterialKind.MOLECULAR, state=rt_scene())
    graph.react("template_switch", "Template-switch at the cDNA 5′ end",
                StepKind.TEMPLATE_SWITCHING,
                ((first_strand, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (switched,))
    smart_primer = graph.input("smart_primer", "single TSO-PCR primer", MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "WTA and library PCR polymerase",
                             MaterialKind.CATALYST)
    wta = graph.product("wta_cdna", "whole-transcriptome-amplified cDNA",
                        MaterialKind.MOLECULAR)
    graph.react("wta_pcr", "Amplify with the single SMART-handle primer", StepKind.PCR,
                ((switched, InputRole.TEMPLATE), (smart_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (wta,),
                geometry=PCRGeometry.unavailable(
                    final_library(), reason="The WTA substrate is symmetric before the final tagmentation shell.",
                    symmetric_species=True))
    tn5_adapter = graph.input("tn5_adapter", "symmetric s7 Tn5 adapter", MaterialKind.ADAPTER,
                              inferred=True, source=inferred_source)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged_cdna", "s7-tagmented WTA fragments", MaterialKind.MOLECULAR,
                           inferred=True, source=inferred_source)
    tag_waste = graph.product("tagmentation_waste", "non-library tagmentation products",
                              MaterialKind.WASTE)
    graph.react("tagmentation", "Tagment WTA with the inferred symmetric s7 adapter",
                StepKind.TAGMENTATION,
                ((wta, InputRole.TEMPLATE), (tn5_adapter, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tag_waste), source=inferred_source)
    graph.discard(tag_waste, "fragments without usable primer sites do not amplify")
    p5 = graph.input("p5_primer", "published bead-end-selective P5 primer", MaterialKind.OLIGO)
    p7 = graph.input("p7_primer", "published N701 indexed P7 primer", MaterialKind.OLIGO)
    final = final_library()
    library = graph.product("library", final.name, MaterialKind.LIBRARY, state=final,
                            inferred=True, source=inferred_source)
    graph.react("library_pcr", "Select the bead end and PCR-complete the library", StepKind.PCR,
                ((tagged, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The final construct is reconstructed from exact PCR primers around an unpublished Tn5 adapter."),
                source=inferred_source)
    graph.sequence(library, SequencingPlan(
        library=final, primers=SEQ_PRIMERS, acquired_roles=RUN_ROLES,
        read_lengths={"Read 1": 54, "Index 1 (i7)": 8},
        read_uses={"Read 1": "extract", "Index 1 (i7)": "demultiplex",
                   "Read 2": "include"},
        limitation="The s7/ME and P7-side construct regions are inferred from the source boundary."))
    graph.declare_reviewed_scope((
        "bead_seed", "barcode_oligos", "extension_polymerase", "bead_assembly",
        "rna", "reverse_transcriptase", "reverse_transcription", "tso",
        "template_switch", "smart_primer", "pcr_polymerase", "wta_pcr",
        "tn5_adapter", "tn5", "tagmentation", "p5_primer", "p7_primer",
        "library_pcr", "library"), source)
    return graph.require_publishable()
