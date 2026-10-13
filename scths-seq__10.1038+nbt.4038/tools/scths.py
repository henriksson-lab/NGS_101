"""Molecular construct model for scTHS-seq (Lake et al., 2018)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


T7_LEADER = "AATTAATACGACTCACTATA"
CONNECTOR = "GGGAGATCCACGCGC"
R5_EXAMPLE = "TCTAAT"
I5_OLIGO_INDEX = "CTCTCTAT"
I7_OLIGO_INDEX = "TCGCCTTA"
I7_INDEX_READ = "TAAGGCGA"


def r5_transposon() -> list[Segment]:
    return [seg("T7 leader/promoter", T7_LEADER, "t7"),
            seg("IVT/i5 connector", CONNECTOR, "t7"),
            seg("round-1 r5 barcode", "B" * 6, "cbc", placeholder=True),
            seg("s5", nx.S5, "s5"), seg("mosaic end", nx.ME, "me")]


def s7_transposon() -> list[Segment]:
    return [seg("s7", nx.S7, "s7"), seg("mosaic end", nx.ME, "me")]


def first_tagmented_scene(insert_nt: int = 28) -> Scene:
    """One end of a genomic fragment tagged by the T7/r5/s5 homodimer."""
    top = [*r5_transposon(), seg("accessible genomic DNA", "X" * insert_nt,
                                placeholder=True)]
    sc = Scene()
    sc.strand("top", top, label="round-1 tagged DNA")
    sc.anneal("ME bottom", [seg("mosaic end'", nx.ME_RC, "me")], to="top",
              pair=("mosaic end'", "mosaic end"), label="", mod5="p")
    sc.mark("top", "mosaic end", "9-nt genomic gap on the opposite strand")
    return sc


def ivt_template(insert_nt: int = 28) -> Construct:
    return Construct([seg("IVT/i5 connector", CONNECTOR, "t7"),
                      seg("round-1 r5 barcode", "B" * 6, "cbc", placeholder=True),
                      seg("s5", nx.S5, "s5"), seg("mosaic end", nx.ME, "me"),
                      seg("genomic DNA", "X" * insert_nt, placeholder=True)],
                     name="scTHS-seq re-made dsDNA")


INDEX2_P5 = sp.custom(
    "Index 2 (i5)", "scTHS nXTv2 i5 index-read primer", il.P5,
    "Lake et al. 2018 Supplementary Table S12",
    "Reads i5, the 15-nt connector, the r5 barcode and three s5 bases.")
SEQ_PRIMERS = (sp.NEXTERA["R1"], sp.NEXTERA["I1"], INDEX2_P5)
RUN_ROLES = ("Read 1", "Index 1 (i7)", "Index 2 (i5)")


def final_library(insert_nt: int = 28) -> Construct:
    lib = Construct([
        seg("P5", il.P5, "p5"), seg("i5 index", I5_OLIGO_INDEX, "cbc",
                                        feature=feature("cell_i5", "cell_barcode", "combinatorial", group="cell_id", part="i5")),
        seg("IVT/i5 connector", CONNECTOR, "t7"),
        seg("round-1 r5 barcode", "B" * 6, "cbc", placeholder=True,
            feature=feature("cell_r5", "cell_barcode", "combinatorial", group="cell_id", part="round 1 r5")),
        seg("s5", nx.S5, "s5"), seg("mosaic end", nx.ME, "me"),
        seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
        seg("second mosaic end reverse complement", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 index read", I7_INDEX_READ, "cbc",
            feature=feature("cell_i7", "cell_barcode", "combinatorial", group="cell_id", part="i7")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="scTHS-seq sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS, required_roles=RUN_ROLES)
    if problems:
        raise ValueError("invalid final scTHS-seq construct: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1038/nbt.4038"
    graph = ProtocolGraph("scTHS-seq", authority=source)
    nuclei = graph.input("nuclei", "permeabilized nuclei in first-round wells",
                         MaterialKind.BIOLOGICAL, primary=True)
    r5_transposome = graph.input("r5_transposome", "barcoded T7/r5/s5 transposome",
                                 MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    first_tagged = graph.product("first_tagged", "round-1 r5-tagged accessible DNA",
                                 MaterialKind.MOLECULAR, state=first_tagmented_scene())
    first_waste = graph.product("first_tagmentation_waste",
                                "non-library first-round tagmentation products",
                                MaterialKind.WASTE)
    graph.react("first_tagmentation",
                "Tagment nuclei with a barcoded T7/r5/s5 homodimer",
                StepKind.TAGMENTATION,
                ((nuclei, InputRole.TEMPLATE), (r5_transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (first_tagged, first_waste))
    graph.discard(first_waste, "molecules without usable first-round insertions do not amplify")
    sort_rule = graph.input("single_nucleus_sort", "single-nucleus sorting rule",
                            MaterialKind.SELECTION)
    sorted_tagged = graph.product("sorted_tagged", "single-nucleus r5-tagged DNA",
                                  MaterialKind.MOLECULAR)
    excluded = graph.product("sort_excluded", "excluded nuclei", MaterialKind.WASTE)
    graph.react("sort", "Sort single nuclei", StepKind.SAMPLING,
                ((first_tagged, InputRole.TEMPLATE), (sort_rule, InputRole.SELECTION)),
                (sorted_tagged, excluded), changes_molecule=False)
    graph.discard(excluded, "events outside the single-nucleus gate are not processed")
    polymerase = graph.input("polymerase", "gap-fill, dsDNA and PCR polymerase",
                             MaterialKind.CATALYST)
    filled = graph.product("gap_filled", "gap-filled T7/r5 insertion templates",
                           MaterialKind.MOLECULAR, state=ivt_template())
    graph.react("gap_fill", "Fill transposase-generated gaps", StepKind.EXTENSION,
                ((sorted_tagged, InputRole.TEMPLATE), (polymerase, InputRole.CATALYST)),
                (filled,))
    t7 = graph.input("t7_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    ivt_rna = graph.product("ivt_rna", "T7-amplified insertion-end RNA",
                            MaterialKind.MOLECULAR)
    graph.react("ivt", "Transcribe each insertion end from its T7 promoter",
                StepKind.TRANSCRIPTION,
                ((filled, InputRole.TEMPLATE), (t7, InputRole.CATALYST)), (ivt_rna,))
    random_primer = graph.input("random_rt_primer", "random reverse-transcription primer",
                                MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product("cdna", "random-primed insertion-end cDNA", MaterialKind.MOLECULAR)
    graph.react("reverse_transcription", "Random-prime the amplified RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((ivt_rna, InputRole.TEMPLATE), (random_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    connector_primer = graph.input("connector_primer", "IVT-connector dsDNA primer",
                                   MaterialKind.OLIGO)
    remade = graph.product("remade_dsdna", "connector-primed re-made dsDNA",
                           MaterialKind.MOLECULAR, state=ivt_template())
    graph.react("dsdna_conversion", "Make dsDNA from the IVT connector",
                StepKind.EXTENSION,
                ((cdna, InputRole.TEMPLATE), (connector_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (remade,))
    s7_transposome = graph.input("s7_transposome", "s7-only second transposome",
                                 MaterialKind.ADAPTER)
    second_tagged = graph.product("second_tagged", "r5-to-s7 doubly tagged fragments",
                                  MaterialKind.MOLECULAR)
    second_waste = graph.product("second_tagmentation_waste",
                                 "non-library second-round products", MaterialKind.WASTE)
    graph.react("second_tagmentation", "Tagment the re-made dsDNA with s7-only Tn5",
                StepKind.TAGMENTATION,
                ((remade, InputRole.TEMPLATE), (s7_transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (second_tagged, second_waste))
    graph.discard(second_waste, "fragments without an r5-to-s7 configuration do not amplify")
    p5 = graph.input("p5_primer", "scT7 S502 i5/P5 primer", MaterialKind.OLIGO)
    p7 = graph.input("p7_primer", "N701 i7/P7 primer", MaterialKind.OLIGO)
    final = final_library()
    library = graph.product("library", final.name, MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "Add i5/P5 at the connector end and i7/P7 at s7",
                StepKind.PCR,
                ((second_tagged, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The final construct is exact, while the page does not expose separate PCR landing objects."))
    graph.sequence(library, SequencingPlan(
        library=final, primers=SEQ_PRIMERS, acquired_roles=RUN_ROLES,
        read_lengths={"Read 1": 50, "Index 1 (i7)": 8, "Index 2 (i5)": 32},
        read_uses={"Read 1": "include", "Index 1 (i7)": "extract",
                   "Index 2 (i5)": "extract"}))
    graph.declare_reviewed_scope((
        "nuclei", "r5_transposome", "tn5", "first_tagmentation",
        "single_nucleus_sort", "sort", "polymerase", "gap_fill", "t7_polymerase",
        "ivt", "random_rt_primer", "reverse_transcriptase", "reverse_transcription",
        "connector_primer", "dsdna_conversion", "s7_transposome",
        "second_tagmentation", "p5_primer", "p7_primer", "library_pcr", "library"),
        source)
    return graph.require_publishable()
