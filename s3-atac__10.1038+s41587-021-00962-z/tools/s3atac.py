"""Molecular construct model for s3-ATAC (Mulqueen et al., 2021).

The representative index sequences are transcribed from Supplementary Tables 2, 4 and
5.  Everything shared with Illumina/Nextera chemistry is assembled from ``lib/`` rather
than copied here.
"""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Disposition, Terminal, InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind)


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


# Representative wells from the source tables.  The protocol uses 96 Tn5, 32 i7 and
# 64 i5 sequences; one real member of each set keeps the molecular drawings legible.
TN5_INDEX = "GAACCGCG"       # SBS12_18_UME_sci_1
I7_OLIGO_INDEX = "TCGCCTTA"  # PCR_i7_P7.S701
I5_OLIGO_INDEX = "CCTTAAGA"  # PCR_A_i5_A
SBS12_PARTIAL = il.TRUSEQ_READ2[-18:]
LNA_ME_POSITIONS = (10, 12, 14, 16, 18)


def tn5_transfer() -> list[Segment]:
    """Transferred strand of the single indexed U-ME adaptor used in a well."""
    return [
        seg("partial TruSeq Read 2", SBS12_PARTIAL, "r2"),
        seg("Tn5 barcode", TN5_INDEX, "cbc"),
        seg("dU", "U", "w1", placeholder=True),
        seg("mosaic end", nx.ME, "me"),
    ]


def tn5_bottom() -> list[Segment]:
    """Phosphorylated, non-transferred strand paired to the mosaic end."""
    return [seg("mosaic end reverse complement", nx.ME_RC, "me")]


def switching_oligo() -> list[Segment]:
    """A14-ME oligo; the ME contains LNA at positions 10, 12, 14, 16 and 18."""
    return [seg("s5", nx.S5, "s5"), seg("LNA mosaic end", nx.ME, "me")]


def i7_pcr_primer() -> list[Segment]:
    return [seg("P7", il.P7, "p7"), seg("i7", I7_OLIGO_INDEX, "cbc"),
            seg("TruSeq Read 2", il.TRUSEQ_READ2, "r2")]


def i5_pcr_primer() -> list[Segment]:
    return [seg("P5", il.P5, "p5"), seg("i5", I5_OLIGO_INDEX, "cbc"),
            seg("s5", nx.S5, "s5")]


def transposome_scene() -> Scene:
    """Loaded adaptor.  ``Scene.anneal`` makes an incorrect ME bottom impossible."""
    sc = Scene()
    sc.strand("transfer", tn5_transfer(), label="transferred")
    sc.anneal("bottom", tn5_bottom(), to="transfer",
              pair=("mosaic end reverse complement", "mosaic end"),
              label="Tn5-bound", mod5="p")
    sc.mark("transfer", "dU", "polymerase stop")
    return sc


def tagmented_scene() -> Scene:
    """A symmetrically tagmented fragment before gap fill.

    Both transferred adaptors are attached to genomic 5' ends.  The 9-nt regions beside
    the joins are left single-stranded, as required by the staggered Tn5 cut.
    """
    gap = nx.TAGMENTATION_GAP
    core = seg("genomic core", "X" * 24, placeholder=True)
    top = [*tn5_transfer(), seg("left 9-nt gap", "X" * gap, placeholder=True), core]
    bottom = [*tn5_transfer(), seg("right 9-nt gap", "X" * gap, placeholder=True),
              *complement_segments([core])]
    sc = Scene()
    sc.strand("top", top, label="DNA")
    sc.anneal("bottom", bottom, to="top", pair=("genomic core'", "genomic core"),
              label="DNA")
    sc.anneal("left ME bottom", tn5_bottom(), to="top",
              pair=("mosaic end reverse complement", "mosaic end"), label="", mod5="p")
    sc.anneal("right ME bottom", tn5_bottom(), to="bottom",
              pair=("mosaic end reverse complement", "mosaic end"), label="", mod5="p",
              above=True)
    sc.mark("top", "left 9-nt gap", "unfilled on the opposite strand")
    sc.mark("bottom", "right 9-nt gap", "unfilled on the opposite strand")
    return sc


def gap_filled_strand(insert_nt: int = 28) -> Construct:
    """One strand after NPM closes its 9-nt gap and stops before the adaptor dU."""
    return Construct([
        *tn5_transfer(),
        seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
        seg("copied mosaic end", nx.ME_RC, "me"),
    ], name="gap-filled s3-ATAC strand")


def switching_scene(insert_nt: int = 28) -> Scene:
    """Blocked A14-ME annealed as the template for addition of copied s5."""
    target = gap_filled_strand(insert_nt)
    sc = Scene()
    sc.strand("target", list(target), label="target")
    sc.anneal("A14-ME", switching_oligo(), to="target",
              pair=("LNA mosaic end", "copied mosaic end"),
              label="switch oligo", mod3="/3InvdT/")
    sc.arrow("target", "NPM copies s5 from the blocked template")
    return sc


def switched_strand(insert_nt: int = 28) -> Construct:
    """Adapter-switched strand produced by the scene above."""
    return Construct([
        *gap_filled_strand(insert_nt),
        seg("copied s5", nx.S5_RC, "s5"),
    ], name="adapter-switched s3-ATAC strand")


SEQ_PRIMERS = tuple([
    sp.NEXTERA["R1"],
    sp.TRUSEQ["I1"],
    sp.NEXTERA["I2"],
    sp.TRUSEQ["R2"],
])


def final_library(insert_nt: int = 28) -> Construct:
    """Final duplex, represented by its P5-to-P7' strand.

    Primer verification is part of construction: a caller cannot obtain a final library
    whose declared standard sequencing primers do not land uniquely.
    """
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("i5", I5_OLIGO_INDEX, "cbc",
            feature=feature("cell_i5", "cell_barcode", "combinatorial", group="cell_id", part="i5")),
        seg("s5", nx.S5, "s5"),
        seg("mosaic end", nx.ME, "me"),
        seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
        seg("opposite mosaic end", nx.ME_RC, "me"),
        seg("base opposite dU", "A", "w1"),
        seg("Tn5 barcode, read orientation", revcomp(TN5_INDEX), "cbc",
            feature=feature("cell_tn5", "cell_barcode", "combinatorial", group="cell_id", part="Tn5")),
        seg("TruSeq Read 2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
        seg("i7, read orientation", revcomp(I7_OLIGO_INDEX), "cbc",
            feature=feature("cell_i7", "cell_barcode", "combinatorial", group="cell_id", part="i7")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="s3-ATAC sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid final s3-ATAC construct: " + "; ".join(problems))
    return lib


def primer_landings(lib: Construct | None = None) -> list[tuple[sp.SeqPrimer, sp.Landing]]:
    """All declared primer sites, refusing a missing or ambiguous landing."""
    lib = lib or final_library()
    out = []
    for primer in SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        if hit is None:
            raise ValueError(f"{primer.name} does not land on {lib.name}")
        out.append((primer, hit))
    return out


def read2_layout() -> list[tuple[str, str]]:
    """Cycle spans before Read 2 reaches genomic DNA, derived from the construct."""
    lib = final_library()
    ordered = [
        ("Tn5 barcode", lib.get("Tn5 barcode, read orientation")),
        ("T copied opposite dU", lib.get("base opposite dU")),
        ("mosaic end", lib.get("opposite mosaic end")),
    ]
    rows, cycle = [], 1
    for label, segment in ordered:
        end = cycle + len(segment) - 1
        rows.append((str(cycle) if cycle == end else f"{cycle}&ndash;{end}", label))
        cycle = end + 1
    rows.append((f"{cycle}&ndash;", "accessible genomic DNA"))
    return rows


def protocol_graph() -> ProtocolGraph:
    source="Defining source: https://doi.org/10.1038/s41587-021-00962-z"
    graph=ProtocolGraph("s3-ATAC",authority=source)
    nuclei=graph.input("nuclei","permeabilized nuclei",MaterialKind.BIOLOGICAL,primary=True)
    transposome=graph.input("indexed_transposome","indexed U-ME transposome",
                            MaterialKind.ADAPTER,state=transposome_scene())
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("tagmented_fragments","symmetrically U-ME-tagmented accessible chromatin",
                         MaterialKind.MOLECULAR,state=tagmented_scene())
    waste=graph.product("tagment_waste","candidate non-library tagmentation products",MaterialKind.MOLECULAR)
    graph.react("tagment","Tagment accessible chromatin with indexed U-ME transposomes",
                StepKind.TAGMENTATION,
                ((nuclei,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,waste))
    graph.mark_terminal(Terminal(waste.id,Disposition.UNRESOLVED,"Candidate insertion products and their physical recovery remain unresolved."))
    npm=graph.input("npm","uracil-intolerant Nextera polymerase mix",MaterialKind.CATALYST)
    dntps=graph.input("gap_fill_dntps","gap-fill dNTPs",MaterialKind.OLIGO)
    filled=graph.product("gap_filled",gap_filled_strand().name,MaterialKind.MOLECULAR,
                         state=gap_filled_strand())
    graph.react("gap_fill","Fill the Tn5 gap until polymerase stops at dU",StepKind.EXTENSION,
                ((tagged,InputRole.TEMPLATE),(dntps,InputRole.STRUCTURAL),
                 (npm,InputRole.CATALYST)),(filled,))
    switch=graph.input("switching_oligo","A14-ME LNA adapter-switching oligo",
                       MaterialKind.OLIGO,state=tuple(switching_oligo()))
    switched=graph.product("switched",switched_strand().name,MaterialKind.MOLECULAR,
                           state=switched_strand())
    graph.react("adapter_switch","Template addition of s5: 10 adapter-switching cycles",
                StepKind.EXTENSION,
                ((filled,InputRole.TEMPLATE),(switch,InputRole.STRUCTURAL),
                 (npm,InputRole.CATALYST)),(switched,))
    i5=graph.input("i5_primer","indexed i5 PCR primer",MaterialKind.OLIGO,
                   state=tuple(i5_pcr_primer()))
    i7=graph.input("i7_primer","indexed i7 PCR primer",MaterialKind.OLIGO,
                   state=tuple(i7_pcr_primer()))
    q5u=graph.input("q5u","uracil-tolerant Q5U polymerase",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final)
    graph.react("library_pcr","Add i5/P5 and i7/P7: 16–18 Q5U PCR cycles",
                StepKind.PCR,
                ((switched,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(q5u,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.exact(
                    switched_strand(),
                    PCRPrimer("i7 PCR primer","".join(s.top for s in i7_pcr_primer()),
                              anneal_nt=len(SBS12_PARTIAL)),
                    PCRPrimer("i5 PCR primer","".join(s.top for s in i5_pcr_primer()),
                              anneal_nt=len(nx.S5)),
                    statement="Supplementary Tables 4-5 print the indexed primers and their partial-Read-2/s5 landing regions."))
    lengths={"Read 1":85,"Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":85}
    graph.sequence(library,SequencingPlan(
        library=final,primers=SEQ_PRIMERS,acquired_roles=tuple(lengths),read_lengths=lengths,
        read_uses={"Read 1":"include","Index 1 (i7)":"extract",
                   "Index 2 (i5)":"extract","Read 2":"extract_then_include"}))
    graph.declare_reviewed_scope((
        "nuclei","indexed_transposome","tn5","tagment","npm","gap_fill_dntps",
        "gap_fill","switching_oligo","adapter_switch","i5_primer","i7_primer",
        "q5u","library_pcr","library"),source)
    return graph.require_publishable()
