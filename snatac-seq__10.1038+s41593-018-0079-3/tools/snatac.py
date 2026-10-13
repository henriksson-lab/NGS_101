"""Molecular model for snATAC-seq (Preissl et al. 2018)."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Segment, feature, revcomp
import illumina as il
import nextera as nx
import seqprimers as sp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)


# Constant portions of the two barcoded transposons, Supplementary Table 5.
P5_LINKER = "TCCACGC"
P7_LINKER = "CTGTCCCTGTCC"
READ1_SITE = "GCGATCGAGGACGGC"
READ2_SITE = "CACCGTCTCCGCCTC"
BARCODE_LEN = 8
P5_BARCODES = 8
P7_BARCODES = 12
PCR_INDEX_LEN = 8

PMENTS = nx.ME_RC
READ1_PRIMER = READ1_SITE + nx.ME
READ2_PRIMER = READ2_SITE + nx.ME
INDEX1_PRIMER = revcomp(READ2_PRIMER)

# Representative verbatim rows retained for the small transcription test.
P5_1 = "TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG"
P7_1 = "GTCTCGTGGGCTCGGCTGTCCCTGTCCCGAGTAATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG"
N701 = "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGG"
S502 = "AATGATACGGCGACCACCGAGATCTACACCTCTCTATTCGTCGGCAGCGTC"
P5_1_BARCODE = "TATAGCCT"
P7_1_BARCODE = "CGAGTAAT"
N701_ORDERED_INDEX = "TCGCCTTA"
S502_INDEX = "CTCTCTAT"


def _seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def barcode(name: str, value: str | None, letter: str) -> Segment:
    """A real or placeholder barcode whose required length is enforced here."""
    seq = letter * BARCODE_LEN if value is None else value
    if len(seq) != BARCODE_LEN:
        raise ValueError(f"{name} must be {BARCODE_LEN} nt")
    role = "cell_barcode" if name in ("p5 barcode", "p7 barcode") else "cell_barcode"
    part = name.replace(" barcode", "")
    return _seg(name, seq, "cbc", placeholder=value is None,
                feature=feature(f"cell_{part.replace(' ', '_')}", role, "combinatorial",
                                group="cell_id", part=part))


def p5_transposon(value: str | None = None) -> list[Segment]:
    return [_seg("s5", nx.S5, "s5"), _seg("p5 linker", P5_LINKER, "r2"),
            barcode("p5 barcode", value, "A"), _seg("Read 1 site", READ1_SITE, "r1"),
            _seg("ME", nx.ME, "me")]


def p7_transposon(value: str | None = None) -> list[Segment]:
    return [_seg("s7", nx.S7, "s7"), _seg("p7 linker", P7_LINKER, "r3"),
            barcode("p7 barcode", value, "B"), _seg("Read 2 site", READ2_SITE, "r1"),
            _seg("ME", nx.ME, "me")]


def tagged_fragment() -> Construct:
    """The amplifiable A...B tagmentation product after its 9-nt gaps are filled."""
    right = p7_transposon()
    return Construct([
        *p5_transposon(),
        _seg("genomic insert", "XXXXXXXX...XXXXXXXX", placeholder=True),
        *[Segment(s.name + "'", revcomp(s.top) if not s.placeholder else s.top.lower(),
                  s.tag, s.placeholder, feature=s.feature) for s in reversed(right)],
    ], name="gap-filled snATAC fragment")


SEQ_PRIMERS = (
    sp.custom("Read 1", "snATAC Read 1", READ1_PRIMER,
              "Preissl et al. Supplementary Table 5"),
    sp.custom("Index 1 (i7)", "snATAC Index 1", INDEX1_PRIMER,
              "Preissl et al. Supplementary Table 5"),
    sp.custom("Index 2 (i5)", "flow-cell P5 oligo", il.P5,
              "HiSeq 2500 forward-strand index workflow",
              "No separate Index 2 sequencing primer is listed."),
    sp.custom("Read 2", "snATAC Read 2", READ2_PRIMER,
              "Preissl et al. Supplementary Table 5"),
)


def final_library() -> Construct:
    """Dual-PCR-indexed, dual-transposon-indexed snATAC library.

    The B-side sequence is generated from the ordered p7 oligo rather than retyped.
    Primer placement is validated here so no caller can obtain an unsequenceable model.
    """
    tagged = tagged_fragment()
    con = Construct([
        _seg("P5", il.P5, "p5"),
        _seg("i5", "N" * PCR_INDEX_LEN, "cbc", placeholder=True,
             feature=feature("cell_i5", "cell_barcode", "combinatorial", group="cell_id", part="i5")),
        *tagged.segments,
        _seg("i7 read", "N" * PCR_INDEX_LEN, "cbc", placeholder=True,
             feature=feature("cell_i7", "cell_barcode", "combinatorial", group="cell_id", part="i7")),
        _seg("P7'", il.P7_RC, "p7"),
    ], name="snATAC-seq library")
    errors = sp.verify(con, SEQ_PRIMERS)
    if errors:
        raise ValueError("invalid snATAC library: " + "; ".join(errors))
    return con


def p5_index_primer(value: str | None = None) -> list[Segment]:
    idx = barcode("i5", value, "N")
    return [_seg("P5", il.P5, "p5"), idx,
            _seg("s5", nx.S5, "s5")]


def p7_index_primer(value: str | None = None) -> list[Segment]:
    idx = barcode("i7 as ordered", value, "N")
    return [_seg("P7", il.P7, "p7"),
            idx,
            _seg("s7", nx.S7, "s7")]


def protocol_graph() -> ProtocolGraph:
    source="Defining source: https://doi.org/10.1038/s41593-018-0079-3"
    graph=ProtocolGraph("snATAC-seq",authority=source)
    nuclei=graph.input("nuclei","isolated nuclei",MaterialKind.BIOLOGICAL,primary=True)
    well_a=graph.product("tagment_well_a","representative tagmentation well A nuclei",
                         MaterialKind.BIOLOGICAL)
    well_b=graph.product("tagment_well_b","representative tagmentation well B nuclei",
                         MaterialKind.BIOLOGICAL)
    graph.react("tagment_split","Split nuclei across indexed tagmentation wells",StepKind.SPLIT,
                ((nuclei,InputRole.TEMPLATE),),(well_a,well_b))
    transposomes=graph.input("indexed_transposomes","paired p5/p7 indexed transposome set",
                             MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged_a=graph.product("tagged_a","round-1-indexed nuclei from well A",MaterialKind.MOLECULAR)
    tagged_b=graph.product("tagged_b","round-1-indexed nuclei from well B",MaterialKind.MOLECULAR)
    waste_a=graph.product("tagment_waste_a","non-library products from well A",MaterialKind.WASTE)
    waste_b=graph.product("tagment_waste_b","non-library products from well B",MaterialKind.WASTE)
    graph.react("tagment_a","Tagment accessible chromatin with well-A p5/p7 indexes",
                StepKind.TAGMENTATION,
                ((well_a,InputRole.TEMPLATE),(transposomes,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged_a,waste_a))
    graph.react("tagment_b","Tagment accessible chromatin with well-B p5/p7 indexes",
                StepKind.TAGMENTATION,
                ((well_b,InputRole.TEMPLATE),(transposomes,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged_b,waste_b))
    graph.discard(waste_a,"same-ended and otherwise unusable transposition products do not amplify")
    graph.discard(waste_b,"same-ended and otherwise unusable transposition products do not amplify")
    pooled=graph.product("pooled_nuclei","pooled round-1-indexed nuclei",MaterialKind.MOLECULAR)
    graph.react("pool","Pool indexed tagmentation wells",StepKind.POOL,
                ((tagged_a,InputRole.TEMPLATE),(tagged_b,InputRole.TEMPLATE)),(pooled,))
    sorter=graph.input("pcr_well_sort","one-nucleus-per-PCR-well sorting rule",MaterialKind.SELECTION)
    wells=graph.product("pcr_wells","single round-1-indexed nuclei in PCR wells",MaterialKind.MOLECULAR)
    excluded=graph.product("excluded_events","empty and multiplet sorting events",MaterialKind.WASTE)
    graph.react("sort","Sort individual nuclei into second-index PCR wells",StepKind.SELECTION,
                ((pooled,InputRole.TEMPLATE),(sorter,InputRole.SELECTION)),(wells,excluded))
    graph.discard(excluded,"empty and multiplet sorting events are excluded")
    gap_pol=graph.input("gap_fill_polymerase","gap-fill and PCR polymerase",MaterialKind.CATALYST)
    filled=graph.product("gap_filled","gap-filled A-to-B tagmentation fragments",
                         MaterialKind.MOLECULAR,state=tagged_fragment())
    graph.react("gap_fill","Fill the 9-nt transposition gaps",StepKind.EXTENSION,
                ((wells,InputRole.TEMPLATE),(gap_pol,InputRole.CATALYST)),(filled,))
    i5=graph.input("i5_primer","second-round i5 PCR primer",MaterialKind.OLIGO,
                   state=tuple(p5_index_primer()))
    i7=graph.input("i7_primer","second-round i7 PCR primer",MaterialKind.OLIGO,
                   state=tuple(p7_index_primer()))
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final)
    graph.react("indexed_pcr","Add the second i5/i7 barcode pair by PCR",StepKind.PCR,
                ((filled,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(gap_pol,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.unavailable(final,
                    reason="The exact tagged intermediate is verified separately; the page shows both primer landings."))
    graph.sequence(library,SequencingPlan(
        library=final,primers=SEQ_PRIMERS,acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
        read_lengths={"Read 1":50,"Index 1 (i7)":43,"Index 2 (i5)":37,"Read 2":50},
        read_uses={"Read 1":"include","Index 1 (i7)":"extract",
                   "Index 2 (i5)":"extract","Read 2":"include"}))
    graph.declare_reviewed_scope((
        "nuclei","tagment_split","indexed_transposomes","tn5","tagment_a",
        "tagment_b","pool","pcr_well_sort","sort","gap_fill_polymerase",
        "gap_fill","i5_primer","i7_primer","indexed_pcr","library"),source)
    return graph.require_publishable()
