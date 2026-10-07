#!/usr/bin/env python3
"""Generate the diagram-centric Tang 2009 chemistry page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import tang2009 as T
from chemdraw import Scene, Segment, annotation_rows, complement_segments, oligo, panel
from page import head, info

OUT = HERE.parent / "tang-2009.html"


def seg(name, top, tag=None, **kw):
    return Segment(name, top, tag, **kw)


def preamble() -> str:
    return """<div class="wrap">
<h1>Tang 2009 &mdash; single-cell mRNA sequencing</h1>
<p><info>Defining protocol: <a href="https://doi.org/10.1038/nmeth.1315">Tang et al.
(2009)</a>. Oligo-dT reverse transcription and TdT tailing make amplified cDNA; an
unbarcoded SOLiD fragment library is then prepared from it.</info></p>
<div class="caveat"><b>Source boundary.</b> The 2009 supplementary oligo table and figures
define the unbarcoded construct. Step order is taken from the authors' detailed 2010
companion protocol; its later barcoded P2 extension is excluded below.</div>
"""


def oligos() -> str:
    rows = [
        oligo("UP1 (reverse transcription)",
              [seg("UP1", T.UP1_HANDLE, "r1"), seg("dT", "T" * T.DT_LEN, "tso")]),
        oligo("UP2 (second strand)",
              [seg("UP2", T.UP2_HANDLE, "r2"), seg("dT", "T" * T.DT_LEN, "tso")]),
        oligo("AUP1", [seg("", T.AUP1, "r1")], mods="NH2"),
        oligo("AUP2", [seg("", T.AUP2, "r2")], mods="NH2"),
        oligo("SOLiD P1, top", [seg("", T.P1_TOP, "p5")]),
        oligo("SOLiD P1, bottom", [seg("", T.P1_BOTTOM, "p5")]),
        oligo("SOLiD P2, top", [seg("", T.P2_TOP, "p7")]),
        oligo("SOLiD P2, bottom", [seg("", T.P2_BOTTOM, "p7")]),
    ]
    return "<h2>Key oligos</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>"


def first_strand() -> str:
    mrna = [seg("RNA", "XXXXXXXXXXXX", placeholder=True),
            seg("poly(A)", "A" * T.DT_LEN, "tso")]
    primer = [seg("UP1 handle", T.UP1_HANDLE, "r1"),
              seg("dT", "T" * T.DT_LEN, "tso")]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA")
    sc.anneal("UP1", primer, to="mRNA", pair=("dT", "poly(A)"), label="UP1")
    sc.arrow("UP1", "reverse transcriptase")
    return panel(sc.rows(), cls="small",
                 caption="(1) UP1 oligo-dT primes the poly(A) tail; reverse transcriptase copies the RNA.")


def tail_and_second_strand() -> str:
    first = [seg("UP1 handle", T.UP1_HANDLE, "r1"),
             seg("dT", "T" * T.DT_LEN, "tso"),
             seg("cDNA", "XXXXXXXXXXXX", placeholder=True),
             seg("TdT poly(A)", "A" * 12, "tso")]
    up2 = [seg("UP2 handle", T.UP2_HANDLE, "r2"), seg("dT", "T" * 12, "tso")]
    sc = Scene()
    sc.strand("first", first, label="first strand")
    sc.anneal("UP2", up2, to="first", pair=("dT", "TdT poly(A)"), label="UP2")
    sc.arrow("UP2", "second-strand synthesis")
    return panel(sc.rows(), cls="small",
                 caption="(2–4) Exonuclease I removes free UP1. TdT adds a second poly(A) tail; UP2 primes it and is extended.")


def amplified() -> str:
    lib = T.amplified_cdna()
    sc = Scene()
    sc.strand("top", list(lib), label="", mod5="NH2")
    sc.anneal("bottom", complement_segments(list(lib)), to="top",
              pair=("UP2 handle'", "UP2 handle"), label="", mod5="NH2")
    sc.labels("top")
    return panel(sc.rows(), cls="long",
                 caption="(5–7) PCR with UP1/UP2, then 5'-amine-blocked AUP1/AUP2. Products of 0.5–3 kb are selected. The amines prevent intact cDNA ends from accepting SOLiD adapters.")


def solid_adapters() -> str:
    p1 = Scene()
    p1.strand("top", [seg("paired", T.P1_TOP, "p5")], label="P1 top")
    p1.anneal("bottom", [seg("paired'", T.P1_BOTTOM[:-2], "p5"),
                          seg("TT", T.P1_BOTTOM[-2:], "p5")],
              to="top", pair=("paired'", "paired"), label="P1 bottom")
    p1.mark("bottom", "TT", "distal 3'-TT overhang")

    p2 = Scene()
    p2.strand("top", [seg("paired", T.P2_TOP[:-2], "p7"),
                       seg("TT", T.P2_TOP[-2:], "p7")], label="P2 top")
    p2.anneal("bottom", [seg("paired'", T.P2_BOTTOM, "p7")],
              to="top", pair=("paired'", "paired"), label="P2 bottom")
    p2.mark("top", "TT", "distal 3'-TT overhang")
    return (panel(p1.rows(), cls="small", caption="The unbarcoded 2009 SOLiD P1 adapter.")
            + panel(p2.rows(), cls="small", caption="The unbarcoded 2009 SOLiD P2 adapter."))


def solid_library() -> str:
    lib = T.final_library()
    duplex = Scene.duplex(lib.segments)
    duplex.strands["top"].label = duplex.strands["bottom"].label = ""
    rows = [*duplex.rows(), *annotation_rows(lib, prefix_width=4)]
    return ("<h2>Fragment library</h2>\n"
            + info("The amplified cDNA is sheared to about 100–110 bp, end-repaired, ligated to unbarcoded SOLiD P1/P2 adapters and nick-translated.")
            + solid_adapters()
            + panel(rows, cls="long",
                    caption="2009 final library: one cell per unbarcoded library. Adapter ligation is not directional, so either cDNA strand may face P1.")
            + "<h3>Sequencing</h3>\n"
            + info("SOLiD single-end reads are 35 or 50 bases and enter the cDNA from P1. The sequencing-primer sequence was not published in the sources read, so none is drawn or declared.")
            + '<details class="sources"><summary>Related source</summary><ul>'
              '<li><a href="https://doi.org/10.1038/nprot.2009.236">Tang et al. 2010 '
              'protocol</a>: adds sixteen 6-nt barcoded P2 adapters for multiplexing; '
              'this is not part of the 2009 schematic.</li>'
              '</ul></details></div>')


def main() -> None:
    parts = [head("Tang 2009 library chemistry"), preamble(), oligos(),
             "<h2>cDNA amplification</h2>", first_strand(), tail_and_second_strand(),
             amplified(), solid_library()]
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
