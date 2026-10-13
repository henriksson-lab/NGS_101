"""Construct model for single-cell reduced-representation bisulfite sequencing.

The wet-lab order is from Guo et al. (Genome Research 2013). The paper identifies
standard premethylated indexed Illumina adapters; their single-index sequences are
published in Illumina's authoritative adapter-sequence document.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from endprep import dA_tailed_scene, repair_and_dA_tail
from adapter_ligation import reported_adapter_end
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

INDEX_NT = 6
MSPI_SITE = "CCGG"
MSPI_OVERHANG = "CG"
CELL_INDEX = feature("cell_i7", "cell_barcode", "whitelist",
                     whitelist="published indexed primer set")


def _seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


# The universal oligo is canonical P5 + Read 1 with their four-base ACAC overlap.
UNIVERSAL = il.TRUSEQ_P5_FULL
# The indexed oligo starts with the i7-primer site, then carries the cell index and P7'.
INDEXED_FIXED_5 = il.INDEX1_PRIMER


def universal_segments(*, inferred: bool = False) -> list[Segment]:
    """Universal premethylated adaptor, split without duplicating the ACAC overlap."""
    return [
        _seg("P5", il.P5, "p5", inferred=inferred),
        _seg("Read 1 arm", il.TRUSEQ_READ1[4:], "r1", inferred=inferred),
    ]


def indexed_segments(*, inferred: bool = False) -> list[Segment]:
    return [
        _seg("Index-primer / Read 2 arm", INDEXED_FIXED_5, "r2", inferred=inferred),
        _seg("i7 cell index", "N" * INDEX_NT, "cbc", placeholder=True,
             inferred=inferred, feature=CELL_INDEX),
        _seg("P7'", il.P7_RC, "p7", inferred=inferred),
    ]


def adapter_scene() -> Scene:
    """The TruSeq fork, placed by its 12-bp duplex stem.

    Scene.anneal validates the stem base by base, so changing either arm to a sequence
    which no longer forms the fork fails while the page is being built.
    """
    stem_n = len(il.STEM)
    universal = [
        _seg("universal fork", UNIVERSAL[:-(stem_n + 1)], "p5"),
        _seg("universal stem", UNIVERSAL[-(stem_n + 1):-1], "r1"),
        _seg("3' T", UNIVERSAL[-1], "r1"),
    ]
    indexed = [
        _seg("indexed stem", INDEXED_FIXED_5[:stem_n], "r2"),
        _seg("indexed fork", INDEXED_FIXED_5[stem_n:], "r2"),
        _seg("i7", "N" * INDEX_NT, "cbc", placeholder=True,
             feature=CELL_INDEX),
        _seg("P7'", il.P7_RC, "p7"),
    ]
    sc = Scene()
    sc.strand("universal", universal, label="universal")
    sc.anneal("indexed", indexed, to="universal",
              pair=("indexed stem", "universal stem"), label="indexed", mod5="p",
              unpaired=("indexed fork", "i7", "P7'"))
    sc.mark("universal", "universal stem", "12-bp TruSeq stem")
    return sc


def digested_fragment_scene() -> Scene:
    """A representative MspI fragment with a two-base 5' overhang at each end."""
    core = [_seg("paired fragment", "GXXXXXXXX...XXXXXXXXCC", placeholder=True)]
    sc = Scene()
    sc.strand("top", [_seg("left 5' CG", MSPI_OVERHANG), *core], label="fragment")
    sc.anneal("bottom", [_seg("right 5' CG", MSPI_OVERHANG),
                         *complement_segments(core)],
              to="top", pair=("paired fragment'", "paired fragment"), label="fragment",
              mod5="p")
    sc.mark("top", "left 5' CG", "MspI 5' overhang")
    sc.mark("bottom", "right 5' CG", "MspI 5' overhang")
    return sc


def a_tailed_fragment_scene() -> Scene:
    """Filled MspI fragment with one unpaired 3'-A at each end."""
    core = Construct([_seg("filled fragment", "CGGXXXXXXXX...XXXXXXXXCCG",
                           placeholder=True)])
    return dA_tailed_scene(repair_and_dA_tail(core), label="fragment")


def final_library() -> Construct:
    """PCR-completed, single-index TruSeq library, top strand 5' to 3'."""
    return Construct([
        *universal_segments(),
        _seg("bisulfite-converted insert", "YGGXXXXXXXX...XXXXXXXXYTG",
             placeholder=True),
        _seg("dA junction", "A"),
        *indexed_segments(),
    ], name="scRRBS final library")


SEQ_PRIMERS = [
    sp.TRUSEQ["R1"],
    sp.TRUSEQ["I1"],
    sp.TRUSEQ["R2"],
]


def _validate_model() -> None:
    if UNIVERSAL != il.P5 + il.TRUSEQ_READ1[4:]:
        raise ValueError("universal adaptor no longer assembles from canonical P5/Read 1")
    if revcomp(il.TRUSEQ_READ2) != "A" + INDEXED_FIXED_5:
        raise ValueError("A-tail plus indexed adaptor no longer forms the Read 2 site")
    # Constructing the Scene validates every paired base in the Y-adaptor stem.
    adapter_scene().rows()
    problems = sp.verify(final_library(), SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("invalid scRRBS sequencing layout: " + "; ".join(problems))


_validate_model()


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1101/gr.161679.113"
    graph = ProtocolGraph("scRRBS", authority=source)
    genome = graph.input("genome", "single-cell genomic DNA",
                         MaterialKind.BIOLOGICAL, primary=True)
    mspi = graph.input("mspi", "MspI restriction enzyme", MaterialKind.CATALYST)
    digested = graph.product("mspi_fragments", "MspI-digested genomic fragments",
                             MaterialKind.MOLECULAR, state=digested_fragment_scene())
    graph.react("mspi_digest", "Digest genomic DNA at C^CGG",
                StepKind.DIGESTION,
                ((genome, InputRole.TEMPLATE), (mspi, InputRole.CATALYST)), (digested,))
    klenow = graph.input("klenow", "Klenow exo- end-repair/dA-tailing enzyme",
                         MaterialKind.CATALYST)
    prepared = graph.product("a_tailed_fragments", "filled and dA-tailed MspI fragments",
                             MaterialKind.MOLECULAR, state=a_tailed_fragment_scene())
    graph.react("end_prepare", "Fill recessed ends and add one 3'-A",
                StepKind.END_REPAIR,
                ((digested, InputRole.TEMPLATE), (klenow, InputRole.CATALYST)),
                (prepared,))
    adapters = graph.input("methylated_adapters", "premethylated indexed TruSeq adapters",
                           MaterialKind.ADAPTER, state=adapter_scene())
    ligase = graph.input("ligase", "DNA ligase", MaterialKind.CATALYST)
    ligated = graph.product("ligated_fragments", "adapter-ligated genomic fragments",
                            MaterialKind.MOLECULAR)
    graph.react("adapter_ligation", "Ligate premethylated indexed Y-adapters",
                StepKind.LIGATION,
                ((prepared, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,),
                geometry=reported_adapter_end(
                    "premethylated TruSeq Y adapter", kind="cohesive",
                    polarity="3-prime", terminal_sequence="T",
                    substrate_scenes=(a_tailed_fragment_scene(),),
                    statement="The standard premethylated TruSeq adapter presents one 3'-dT to each dA-tailed insert end."))
    bisulfite = graph.input("bisulfite", "bisulfite-conversion conditions",
                            MaterialKind.CONDITION)
    converted = graph.product("converted_library", "bisulfite-converted ligation products",
                              MaterialKind.MOLECULAR)
    graph.react("bisulfite_conversion", "Bisulfite-convert adapter-ligated molecules",
                StepKind.CHEMICAL_CONVERSION,
                ((ligated, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)),
                (converted,))
    qp1 = graph.input("qp1", "QP1 / P5 PCR primer", MaterialKind.OLIGO)
    qp2 = graph.input("qp2", "QP2 / P7 PCR primer", MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "uracil-tolerant PCR polymerase",
                             MaterialKind.CATALYST)
    final = final_library()
    library = graph.product("library", final.name, MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "Amplify converted molecules and complete flow-cell arms",
                StepKind.PCR,
                ((converted, InputRole.TEMPLATE), (qp1, InputRole.STRUCTURAL),
                 (qp2, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The converted pre-PCR intermediate is not represented base-for-base; the final canonical arms are exact."))
    graph.sequence(library, SequencingPlan(
        library=final, primers=tuple(SEQ_PRIMERS),
        acquired_roles=("Read 1", "Index 1 (i7)", "Read 2"),
        read_lengths={"Index 1 (i7)": INDEX_NT},
        read_uses={"Read 1": "include", "Index 1 (i7)": "demultiplex",
                   "Read 2": "include"}))
    graph.declare_reviewed_scope((
        "genome", "mspi", "mspi_digest", "klenow", "end_prepare",
        "methylated_adapters", "ligase", "adapter_ligation", "bisulfite",
        "bisulfite_conversion", "qp1", "qp2", "pcr_polymerase", "library_pcr",
        "library"), source)
    return graph.require_publishable()
