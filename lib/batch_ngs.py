"""Reusable constructs for the protocol-class gap batch.

The helpers here make final-library completeness a construction invariant: a short-read
library is returned only after its declared sequencing primers have been located.
"""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Row, Scene, Segment, complement_segments, feature


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name, top, tag, **kw)


def truseq_library(insert: list[Segment], name: str, *, dual_index: bool = True,
                    inferred_adapters: bool = False) -> tuple[Construct, tuple]:
    """Canonical TruSeq library with validated R1/I1/(I2)/R2 sites."""
    left = [seg("P5", il.P5, "p5", inferred=inferred_adapters)]
    if dual_index:
        left.append(seg("i5", "J" * 8, "cbc", placeholder=True,
                        feature=feature("sample_index_i5", "sample_index", "unknown"),
                        inferred=inferred_adapters))
    left.append(seg("Read 1 arm", il.TRUSEQ_READ1, "r1", inferred=inferred_adapters))
    lib = Construct([*left, *insert, seg("dA junction", "A", inferred=inferred_adapters),
        seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2", inferred=inferred_adapters),
        seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True,
            feature=feature("sample_index_i7", "sample_index", "unknown"),
            inferred=inferred_adapters),
        seg("P7 reverse complement", il.P7_RC, "p7", inferred=inferred_adapters)], name=name)
    primers = ((sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"], sp.TRUSEQ["R2"])
               if dual_index else (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"]))
    roles = tuple(p.role for p in primers)
    problems = sp.verify(lib, primers, required_roles=roles)
    if problems:
        raise ValueError(f"{name}: " + "; ".join(problems))
    return lib, primers


def nextera_library(insert: list[Segment], name: str) -> tuple[Construct, tuple]:
    """Heterologous S5/S7 Tn5 product after indexed PCR, with validated run sites."""
    lib = Construct([
        seg("P5", il.P5, "p5"), seg("i5", "J" * 8, "cbc", placeholder=True,
                                    feature=feature("sample_index_i5", "sample_index", "unknown")),
        seg("S5", nx.S5, "s5"), seg("left mosaic end", nx.ME, "me"), *insert,
        seg("right mosaic end reverse complement", nx.ME_RC, "me"),
        seg("S7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True,
            feature=feature("sample_index_i7", "sample_index", "unknown")),
        seg("P7 reverse complement", il.P7_RC, "p7")], name=name)
    primers = (sp.NEXTERA["R1"], sp.NEXTERA["I1"], sp.NEXTERA["I2"], sp.NEXTERA["R2"])
    problems = sp.verify(lib, primers)
    if problems:
        raise ValueError(f"{name}: " + "; ".join(problems))
    return lib, primers


def duplex_fragment(label: str = "genomic fragment", length: int = 34) -> Scene:
    return Scene.duplex([seg(label, "X" * length, placeholder=True)], label=label)


def joined_scene(parts: list[Segment], junctions: tuple[tuple[str, str, str], ...],
                 *, label: str, duplex: bool = True) -> Scene:
    """Draw covalently joined parts; ``Scene.junction`` validates every boundary."""
    if duplex:
        sc = Scene.duplex(parts, label=label)
        strand = "top"
    else:
        sc = Scene()
        strand = "molecule"
        sc.strand(strand, parts, label=label)
    for left, right, kind in junctions:
        sc.junction(strand, left, right, kind)
    sc.labels(strand)
    return sc


def spatial_rt_scene(*, barcode_parts: tuple[tuple[str, int], ...], umi: int,
                     surface: str) -> Scene:
    """Poly(A)-primed spatial RT; barcode provenance is encoded in named segments."""
    mrna = [seg("transcript", "X" * 28, placeholder=True), seg("poly(A)", "A" * 18)]
    primer = [seg(surface, "X" * 6, placeholder=True)]
    primer.extend(
        seg(name, "B" * n, "cbc", placeholder=True,
            feature=feature(f"spatial_barcode_{i}", "spatial_barcode", "combinatorial",
                            group="spatial_barcode", part=name))
        for i, (name, n) in enumerate(barcode_parts, 1)
    )
    if umi:
        primer.append(seg("UMI", "U" * umi, "umi", placeholder=True,
                          feature=feature("umi", "umi", "random")))
    primer.append(seg("poly(dT)", "T" * 18))
    sc = Scene(); sc.strand("mRNA", mrna, label="tissue mRNA")
    sc.anneal("capture primer", primer, to="mRNA", pair=("poly(dT)", "poly(A)"),
              label="spatial capture oligo")
    sc.arrow("capture primer", "reverse transcription")
    sc.labels("capture primer")
    return sc


def targeting_rows(enzyme: str, result: str) -> list[Row]:
    return [
        Row(chunks=[("chromatin ---- [target epitope] ---- chromatin", None, False)]),
        Row(chunks=[("                         ^ primary antibody", "cbc", False)]),
        Row(chunks=[(f"                         | Protein A–{enzyme}", "me", False)]),
        Row(chunks=[(f"                         └─ {result}", None, False)]),
    ]


def duplex_family_rows() -> list[Row]:
    return [
        Row(chunks=[("original duplex tags:       α — insert strand A — β", "umi", False)]),
        Row(chunks=[("complementary partner:      β — insert strand B — α", "umi", False)]),
        Row(chunks=[("PCR descendants → SSCS(A)       PCR descendants → SSCS(B)", None, False)]),
        Row(chunks=[("                    SSCS(A) + SSCS(B) → duplex consensus", None, False)]),
    ]
