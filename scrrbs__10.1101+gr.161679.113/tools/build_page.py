#!/usr/bin/env python3
"""Build the diagram-centric scRRBS chemistry page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import illumina as il
import scrrbs as S
import seqprimers as sp
from chemdraw import (Construct, Scene, Segment, annotation_rows, junction_row, oligo, panel,
                      strand_row)
from page import head, info

OUT = HERE.parent / "scrrbs.html"


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def oligos() -> str:
    rows = [
        oligo("premethylated universal TruSeq adaptor", S.universal_segments(),
              mods="all C = 5mC"),
        oligo("premethylated indexed TruSeq adaptor", S.indexed_segments(),
              mods="/5Phos/; all C = 5mC"),
        oligo("QP1 / P5 PCR primer", [seg("P5", il.P5, "p5")], tm_segments="P5"),
        oligo("QP2 / P7 PCR primer", [seg("P7", il.P7, "p7")], tm_segments="P7"),
    ]
    return f"<h2>Key oligos</h2><seq>{''.join(rows)}</seq>"


def digest_and_repair() -> str:
    substrate = Scene.duplex([
        seg("left MspI site", "CCGG"),
        seg("retained genomic fragment", "XXXXXXXX...XXXXXXXX", placeholder=True),
        seg("right MspI site", "CCGG"),
    ])
    substrate.mark("top", "left MspI site", "C^CGG")
    substrate.mark("top", "right MspI site", "C^CGG")
    return f"""<h2>Library construction</h2>
<h3>(1) Lyse one cell and digest with MspI</h3>
{panel(substrate.rows(), cls="small",
       caption="MspI selects fragments bounded by C^CGG sites.")}
{panel(S.digested_fragment_scene().rows(), cls="small",
       caption="The retained fragment has a two-base 5'-CG overhang at each end.")}
<h3>(2) Fill the recessed ends and add one 3'-A</h3>
{panel(S.a_tailed_fragment_scene().rows(), cls="small",
       caption="Klenow exo- fills the MspI ends and dA-tails both 3' ends.")}"""


def ligate_and_convert() -> str:
    before = [
        seg("protected left adaptor", "GATCT", "r1"),
        seg("MspI end", "CGG"),
        seg("genomic DNA", "CXCXXXXXXXX...XXXXXXXXC", placeholder=True),
        seg("filled end", "CGA"),
        seg("protected right adaptor", "GATC", "r2"),
    ]
    after = [
        seg("protected left adaptor", "GATCT", "r1"),
        seg("converted MspI end", "UGG", placeholder=True),
        seg("converted genomic DNA", "UXUXXXXXXXX...XXXXXXXXU", placeholder=True),
        seg("converted filled end", "UGA", placeholder=True),
        seg("protected right adaptor", "GATC", "r2"),
    ]
    before_con = Construct(before)
    after_con = Construct(after)
    return f"""<h3>(3) Ligate the premethylated indexed Y-adaptor</h3>
{panel(S.adapter_scene().rows(), cls="small",
       caption="The standard premethylated universal and indexed oligos meet only in "
               "the 12-bp TruSeq stem; "
               "3'-T joins each dA-tailed insert end.")}
<h3>(4) Bisulfite-convert the ligated molecules</h3>
{panel([strand_row(before_con, prefix="before  5'- ", suffix=" -3'"),
        junction_row(before_con, "protected left adaptor", "MspI end", prefix_width=12),
        junction_row(before_con, "filled end", "protected right adaptor", prefix_width=12),
        strand_row(after_con, prefix="after   5'- ", suffix=" -3'")],
       cls="small", caption="Representative unmethylated C becomes U; methylated genomic C "
                            "stays C. Premethylated adaptor C is protected. PCR later copies "
                            "U as T.")}"""


def amplify_and_sequence() -> str:
    lib = S.final_library()
    sc = Scene.duplex(list(lib))
    return f"""<h3>(5) Uracil-tolerant PCR, then indexed PCR</h3>
{panel([*sc.rows(), *annotation_rows(lib)], cls="long",
       caption="PCR-completed single-index TruSeq library. "
               "The six-base i7 identifies the cell; there is no UMI or i5 index.")}
{sp.section(lib, S.SEQ_PRIMERS, intro="Paired-end TruSeq sequencing. Read 1 starts at the "
            "MspI end; the i7 read identifies the cell. Y denotes a bisulfite call (C/T).",
            required_roles=("Read 1", "Index 1 (i7)", "Read 2"))}"""


def render() -> str:
    return "\n".join([
        head("scRRBS library chemistry"), '<div class="wrap">',
        '<h1>scRRBS &mdash; single-cell reduced-representation bisulfite sequencing</h1>',
        info('Defining source: <a href="https://doi.org/10.1101/gr.161679.113">Guo et al., '
             '<i>Genome Research</i> (2013)</a>. Detailed protocol: '
             '<a href="https://doi.org/10.1038/nprot.2015.039">Guo et al., '
             '<i>Nature Protocols</i> (2015)</a>.'),
        info('The paper specifies standard premethylated indexed Illumina adapters. '
             'The displayed six-base single-index sequences are the corresponding '
             'vendor-published TruSeq adapters; there is no i5 index.'),
        oligos(), digest_and_repair(), ligate_and_convert(), amplify_and_sequence(),
        '</div>',
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
