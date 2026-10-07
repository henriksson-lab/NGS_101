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
from chemdraw import Construct, Scene, Segment, annotation_rows, oligo, panel, strand_row
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
        oligo("QP1 / P5 PCR primer", [seg("P5", il.P5, "p5")]),
        oligo("QP2 / P7 PCR primer", [seg("P7", il.P7, "p7")]),
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
        seg("protected left adaptor", "GATCT", "r1", inferred=True),
        seg("MspI end", "CGG"),
        seg("genomic DNA", "CXCXXXXXXXX...XXXXXXXXC", placeholder=True),
        seg("filled end", "CGA"),
        seg("protected right adaptor", "GATC", "r2", inferred=True),
    ]
    after = [
        seg("protected left adaptor", "GATCT", "r1", inferred=True),
        seg("converted MspI end", "UGG", placeholder=True),
        seg("converted genomic DNA", "UXUXXXXXXXX...XXXXXXXXU", placeholder=True),
        seg("converted filled end", "UGA", placeholder=True),
        seg("protected right adaptor", "GATC", "r2", inferred=True),
    ]
    return f"""<h3>(3) Ligate the premethylated indexed Y-adaptor</h3>
{panel(S.adapter_scene().rows(), cls="small",
       caption="INFERRED — exact adaptor bases come from the secondary reconstruction. "
               "The universal and indexed oligos meet only in the 12-bp TruSeq stem; "
               "3'-T joins each dA-tailed insert end.")}
<h3>(4) Bisulfite-convert the ligated molecules</h3>
{panel([strand_row(Construct(before), prefix="before  5'- ", suffix=" -3'"),
        strand_row(Construct(after), prefix="after   5'- ", suffix=" -3'")],
       cls="small", caption="Representative unmethylated C becomes U; methylated genomic C "
                            "stays C. Premethylated adaptor C is protected. PCR later copies "
                            "U as T.")}"""


def amplify_and_sequence() -> str:
    lib = S.final_library()
    sc = Scene.duplex(list(lib))
    return f"""<h3>(5) Uracil-tolerant PCR, then indexed PCR</h3>
{panel([*sc.rows(), *annotation_rows(lib)], cls="long",
       caption="INFERRED — PCR-completed library assembled from the reconstructed adaptors. "
               "The six-base i7 identifies the cell; there is no UMI or i5 index.")}
{sp.section(lib, S.SEQ_PRIMERS, intro="Paired-end TruSeq sequencing. Read 1 starts at the "
            "MspI end; the i7 read identifies the cell. Y denotes a bisulfite call (C/T).")}"""


def render() -> str:
    return "\n".join([
        head("scRRBS library chemistry"), '<div class="wrap">',
        '<h1>scRRBS &mdash; single-cell reduced-representation bisulfite sequencing</h1>',
        info('Defining source: <a href="https://doi.org/10.1101/gr.161679.113">Guo et al., '
             '<i>Genome Research</i> (2013)</a>. Detailed protocol: '
             '<a href="https://doi.org/10.1038/nprot.2015.039">Guo et al., '
             '<i>Nature Protocols</i> (2015)</a>.'),
        '<div class="caveat"><b>Adaptor sequence.</b> The defining paper specifies '
        'premethylated indexed Illumina adaptors but does not print their bases. Dotted '
        'adaptor regions use the TruSeq sequences recorded by the upstream '
        'scg_lib_structs reconstruction.</div>',
        oligos(), digest_and_repair(), ligate_and_convert(), amplify_and_sequence(),
        '</div>',
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
