#!/usr/bin/env python3
"""Generate scbs-seq.html -- the scBS-seq chemistry page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import illumina as il
import scbs as S
import seqprimers as sp
from chemdraw import (Construct, Scene, Segment, annotation_rows, complement_segments,
                      oligo, panel, strand_row)
from page import head

OUT = HERE.parent / "scbs-seq.html"


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


def preamble() -> str:
    return """<div class="wrap">
<h1>scBS-seq &mdash; post-bisulfite random priming</h1>
<p><info>Protocol: <a href="https://doi.org/10.1038/nmeth.3035">Smallwood et al.,
<i>Nature Methods</i> (2014)</a>. Bisulfite conversion comes first; two random-primer
handles are then added without ligation.</info></p>
<p><info>The indexed arm is the eight-base iPCRTag design cited by the paper and printed
in Quail et al. 2012 Supplementary Table&nbsp;1. This is the obsolete paired-end layout,
not a modern TruSeq dual-index library.</info></p>
"""


def oligos() -> str:
    pe1 = [seg("P5", il.P5, "p5"),
           seg("Read 1 remainder", il.TRUSEQ_READ1[4:], "r1")]
    rows = [
        oligo("oligo1", S.oligo1_segments(), mods="[Btn]"),
        oligo("oligo2", S.oligo2_segments()),
        oligo("PE1.0 forward PCR primer", pe1),
        oligo("indexed iPCRTag reverse primer", S.ipcrtag_primer_segments()),
    ]
    return f"""<h2>Oligos</h2><seq>{''.join(rows)}</seq>"""


def step_bisulfite() -> str:
    con = Construct([seg("converted genomic DNA", "XXXXXXXXXX...XXXXXXXXXX",
                         placeholder=True)], name="converted DNA")
    return f"""<h3>(1) Bisulfite-convert the single-cell genome</h3>
<p><info>Conversion changes unmethylated C to U and fragments the genome. Adapter
addition happens afterwards, avoiding loss of pre-ligated molecules.</info></p>
{panel([strand_row(con)], caption="Single-stranded, fragmented bisulfite-converted DNA.")}"""


def step_oligo1() -> str:
    product = Construct([*S.oligo1_segments(),
                         seg("copied insert", "XXXXXXXXXX...XXXXXXXXXX", placeholder=True)],
                        name="oligo1 product")
    return f"""<h3>(2) Five rounds of oligo1 priming</h3>
<p><info>The random 3' N9 anneals across converted DNA; Klenow exo&minus; extends it.
Denaturation and fresh primer repeat this four more times.</info></p>
{panel([strand_row(product), *annotation_rows(product)],
       caption="One displaced, biotin-tagged oligo1 product, 5' to 3'.")}"""


def step_oligo2() -> str:
    top = [seg("oligo1 handle", S.OLIGO1_HANDLE, "r1"),
           seg("oligo1 N9", "N" * S.RANDOM_NT, placeholder=True),
           seg("insert", "XXXXXXXXXX...X", placeholder=True),
           seg("oligo2 site", "X" * S.RANDOM_NT, placeholder=True)]
    bottom = [seg("oligo2 handle", S.OLIGO2_HANDLE, "r2"),
              seg("oligo2 N9", "N" * S.RANDOM_NT, placeholder=True),
              *complement_segments(top[:-1])]
    sc = Scene()
    sc.strand("captured product", top, label="biotin product", mod5="[Btn]")
    sc.anneal("oligo2 product", bottom, to="captured product",
              pair=("oligo1 handle'", "oligo1 handle"), label="new strand")
    sc.mark("oligo2 product", "oligo2 handle", "second PCR handle")
    drawing = panel(sc.rows(), caption="Duplex after oligo2 extension. N9 sequences are "
                    "random-primer sequence, not molecular identifiers.")
    return f"""<h3>(3) Capture, wash, then one round of oligo2 priming</h3>
<p><info>Streptavidin beads retain the oligo1 products. Oligo2 primes on-bead and
extends through the captured strand, putting a different handle on the other end.</info></p>
{drawing}"""


def final_library() -> str:
    lib = S.final_library()
    sc = Scene.duplex(list(lib))
    drawing = panel([*sc.rows(), *annotation_rows(lib)],
                    caption="Final single-index paired-end library using the historical "
                    "eight-base iPCRTag arm.")
    return f"""<h3>(4) Indexed PCR on the beads</h3>
<p><info>PE1.0 adds the P5 / Read&nbsp;1 end. The indexed iPCRTag primer adds the
single-cell sample index and P7 end.</info></p>
{drawing}"""


def sequencing() -> str:
    primer_table = sp.section(
        S.final_library(), S.SEQ_PRIMERS,
        intro="The dedicated Quail iPCRTag index primer reads one eight-base i7 index; "
              "there is no i5 index.",
        required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    return f"""<h2>Sequencing</h2>
<p><info>100-bp paired-end sequencing. Each read begins with a random N9, which is
clipped before non-directional bisulfite alignment; the index read identifies the cell.</info></p>
{primer_table}
</div>"""


def main() -> None:
    html = "\n".join([head("scBS-seq library chemistry"), preamble(), oligos(),
                      "<h2>Library construction</h2>", step_bisulfite(), step_oligo1(),
                      step_oligo2(), final_library(), sequencing()])
    OUT.write_text(html, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
