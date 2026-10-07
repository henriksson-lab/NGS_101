"""Molecular construct model for s3-ATAC (Mulqueen et al., 2021).

The representative index sequences are transcribed from Supplementary Tables 2, 4 and
5.  Everything shared with Illumina/Nextera chemistry is assembled from ``lib/`` rather
than copied here.
"""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, revcomp


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
        seg("i5", I5_OLIGO_INDEX, "cbc"),
        seg("s5", nx.S5, "s5"),
        seg("mosaic end", nx.ME, "me"),
        seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
        seg("opposite mosaic end", nx.ME_RC, "me"),
        seg("base opposite dU", "A", "w1"),
        seg("Tn5 barcode, read orientation", revcomp(TN5_INDEX), "cbc"),
        seg("TruSeq Read 2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
        seg("i7, read orientation", revcomp(I7_OLIGO_INDEX), "cbc"),
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
