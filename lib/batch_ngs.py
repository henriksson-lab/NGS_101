"""Reusable constructs for the protocol-class gap batch.

The helpers here make final-library completeness a construction invariant: a short-read
library is returned only after its declared sequencing primers have been located.
"""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Row, Scene, Segment, complement_segments


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name, top, tag, **kw)


def truseq_library(insert: list[Segment], name: str, *, dual_index: bool = True,
                    inferred_adapters: bool = False) -> tuple[Construct, tuple]:
    """Canonical TruSeq library with validated R1/I1/(I2)/R2 sites."""
    left = [seg("P5", il.P5, "p5", inferred=inferred_adapters)]
    if dual_index:
        left.append(seg("i5", "J" * 8, "cbc", placeholder=True,
                        inferred=inferred_adapters))
    left.append(seg("Read 1 arm", il.TRUSEQ_READ1, "r1", inferred=inferred_adapters))
    lib = Construct([*left, *insert, seg("dA junction", "A", inferred=inferred_adapters),
        seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2", inferred=inferred_adapters),
        seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True,
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
        seg("P5", il.P5, "p5"), seg("i5", "J" * 8, "cbc", placeholder=True),
        seg("S5", nx.S5, "s5"), seg("left mosaic end", nx.ME, "me"), *insert,
        seg("right mosaic end reverse complement", nx.ME_RC, "me"),
        seg("S7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True),
        seg("P7 reverse complement", il.P7_RC, "p7")], name=name)
    primers = (sp.NEXTERA["R1"], sp.NEXTERA["I1"], sp.NEXTERA["I2"], sp.NEXTERA["R2"])
    problems = sp.verify(lib, primers)
    if problems:
        raise ValueError(f"{name}: " + "; ".join(problems))
    return lib, primers


def duplex_fragment(label: str = "genomic fragment", length: int = 34) -> Scene:
    return Scene.duplex([seg(label, "X" * length, placeholder=True)], label=label)


def spatial_rt_scene(*, barcode_parts: tuple[tuple[str, int], ...], umi: int,
                     surface: str) -> Scene:
    """Poly(A)-primed spatial RT; barcode provenance is encoded in named segments."""
    mrna = [seg("transcript", "X" * 28, placeholder=True), seg("poly(A)", "A" * 18)]
    primer = [seg(surface, "X" * 6, placeholder=True),
              *(seg(name, "B" * n, "cbc", placeholder=True) for name, n in barcode_parts)]
    if umi:
        primer.append(seg("UMI", "U" * umi, "umi", placeholder=True))
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
