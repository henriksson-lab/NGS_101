#!/usr/bin/env python3
"""Build the sparse itChIP-seq custom-Nextera chemistry page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]

import itchip as I
import seqprimers as sp
from chemdraw import Scene, annotation_rows, oligo, panel
from page import head, info

OUT = HERE.parent / "itchip-seq.html"


def key_oligos() -> str:
    return ("<h2>Key oligos</h2><seq>" +
            oligo("representative T5-1 transposon", I.t5_segments()) +
            oligo("representative T7-1 transposon", I.t7_segments()) +
            oligo("common annealing primer", [I.seg("ME'", I.COMMON_ANNEALING, "me")],
                  mods="/5Phos/") +
            oligo("i5-N501 PCR primer", I.indexed_p5_segments()) +
            oligo("i7-N701 PCR primer", I.indexed_p7_segments()) + "</seq>")


def construction() -> str:
    gap = I.gap_filled_fragment()
    final = I.final_library()
    return f"""<h2>Library construction</h2>
<h3>(1) Anneal each barcoded transposon to the common ME strand and load Tn5</h3>
{panel(I.loaded_end_scene("T5").rows(), cls="long", caption="Representative T5 end.")}
{panel(I.loaded_end_scene("T7").rows(), cls="long", caption="Representative T7 end.")}
<h3>(2) Tagment opened chromatin with one T5/T7 barcode pair per well</h3>
<h3>(3) Pool, release chromatin, immunoprecipitate, reverse-crosslink and purify DNA</h3>
<h3>(4) Gap-fill and amplify with one i5 and one i7 primer</h3>
{panel([*Scene.duplex(list(gap)).rows(), *annotation_rows(gap)], cls="long",
       caption="Representative heterotypic T5--T7 product after gap fill. Only a fragment "
               "with one end of each type accepts both PCR primers.")}
{panel([*Scene.duplex(list(final)).rows(), *annotation_rows(final)], cls="long",
       caption="PCR-completed custom-Nextera library. The T5/T7 barcode pair identifies "
               "the tagmentation well; i5/i7 index the amplified pool.")}
{sp.section(final, I.SEQ_PRIMERS, heading="Final library and read layout",
            intro="Paired-end 150-bp sequencing uses four custom primers. Read 1 and Read "
                  "2 enter genomic DNA; Index 1 crosses the T7 barcode and Index 2 crosses "
                  "the T5 barcode. There is no UMI.")}"""


def render() -> str:
    return "\n".join([
        head("itChIP-seq library chemistry"), '<div class="wrap">',
        '<h1>itChIP-seq &mdash; barcoded tagmentation before pooled ChIP</h1>',
        info('Defining source: <a href="https://doi.org/10.1038/s41556-019-0383-5">'
             'Ai et al., <i>Nature Cell Biology</i> (2019)</a>.'),
        '<div class="caveat"><b>Source boundary.</b> This page follows the custom-Nextera '
        'library documented by Supplementary Table 1 and Supplementary Fig. 4. The exact '
        'oligos are published; unavailable main-text Methods leave some reaction ordering '
        'and conditions outside this schematic.</div>',
        key_oligos(), construction(), '</div>',
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
