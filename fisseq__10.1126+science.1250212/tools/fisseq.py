"""Molecular model for the Church-lab transcriptome-wide RNA FISSEQ protocol."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import seqprimers as sp
from chemdraw import Construct, Scene, Segment
from circular import circularize_ssdna, rolling_circle


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


ADAPTER = "TCTCGGGAACGCTGAAGA"
RANDOM_HEXAMER = "N" * 6
CDNA = "N" * 36
RT_PRIMER_WRITTEN = "/5phos/TCTCGGGAACGCTGAAGANNNNNN"
RCA_PRIMER = "TCTTCAGCGTTCCCGAGA"
RCA_PRIMER_WRITTEN = "TCTTCAGCGTTCCCGA*G*A"
SEQUENCING_PRIMERS = (
    ("N", ADAPTER),
    ("N-1", ADAPTER[1:]),
    ("N-2", ADAPTER[2:]),
    ("N-3", ADAPTER[3:]),
    ("N-4", ADAPTER[4:]),
)
CIRCLIGASE_C = 60
CIRCLIGASE_MIN = 60


def rt_scene() -> Scene:
    target = [seg("RNA upstream", "N" * 18, placeholder=True),
              seg("random landing site", "N" * 6, placeholder=True),
              seg("RNA downstream", "N" * 18, placeholder=True)]
    primer = [seg("5-prime adapter", ADAPTER, "r1"),
              seg("random hexamer", RANDOM_HEXAMER, "umi", placeholder=True)]
    sc = Scene(); sc.strand("RNA", target, label="cellular RNA")
    sc.anneal("RT primer", primer, to="RNA",
              pair=("random hexamer", "random landing site"),
              unpaired=("5-prime adapter",), label="5'-p RT primer")
    sc.arrow("RT primer", "M-MuLV copies RNA; aminoallyl-dUTP marks cDNA")
    return sc


def linear_cdna() -> Construct:
    return Construct([
        seg("adapter", ADAPTER, "r1"),
        seg("random hexamer", RANDOM_HEXAMER, "umi", placeholder=True),
        seg("cDNA", CDNA, placeholder=True),
    ], name="cross-linked single-stranded cDNA")


CIRCLE = circularize_ssdna(linear_cdna(), five_prime_phosphate=True)
RCA = rolling_circle(CIRCLE, RCA_PRIMER, copies=3)


def rca_priming_scene() -> Scene:
    template = list(linear_cdna())
    primer = [seg("RCA primer", RCA.primer, "r2")]
    sc = Scene(); sc.strand("circle opened for drawing", template, label="cDNA circle")
    sc.anneal("RCA primer", primer, to="circle opened for drawing",
              pair=("RCA primer", "adapter"), label="RCA primer")
    sc.arrow("RCA primer", "Phi29 copies repeatedly around the circle")
    return sc


def rolony() -> Construct:
    p, c = len(RCA.primer), len(CDNA)
    unit = RCA.unit
    one = (unit[:p], unit[p:p + c], unit[p + c:])
    parts = []
    for i in range(1, RCA.copies + 1):
        parts.extend([
            seg(f"adapter complement {i}", one[0], "r1"),
            seg(f"transcript copy {i}", one[1], placeholder=True),
            seg(f"random hexamer copy {i}", one[2], "umi", placeholder=True),
        ])
    product = Construct(parts, name="RCA concatemer (three repeats drawn)")
    if product.top() != RCA.sequence:
        raise AssertionError("drawn rolony must come from the RCA product")
    return product


SEQ_N = sp.custom("SOLiD sequencing", "FISSEQ primer N", ADAPTER,
                  "Lee et al., Nature Protocols 2015")


def sequencing_scene() -> Scene:
    product = rolony()
    hits = sp.locate_all(product, SEQ_N)
    if len(hits) != RCA.copies:
        raise ValueError("primer N must bind once per RCA repeat")
    target = hits[len(hits) // 2].covers[0]
    sc = Scene(); sc.strand("rolony", list(product), label="RCA concatemer")
    for name, sequence in SEQUENCING_PRIMERS:
        strand = f"primer {name}"
        segment = f"sequencing primer {name}"
        sc.anneal(strand, [seg(segment, sequence, "r2")], to="rolony",
                  pair=(segment, target), label=strand, mod5="p")
    sc.arrow("primer N", "all five offsets ligate probes toward adjacent cDNA")
    return sc
