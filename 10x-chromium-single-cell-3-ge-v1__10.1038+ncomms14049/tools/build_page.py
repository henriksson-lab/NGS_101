#!/usr/bin/env python3
"""Build the 10x Chromium Single Cell 3' Gene Expression v1 page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import chromium_v1 as C
import seqprimers as sp
from chemdraw import annotation_rows, panel, strand_row
from page import head, info, table

OUT = HERE.parent / "10x-chromium-single-cell-3-ge-v1.html"


def construction() -> str:
    bead = C.bead_oligo()
    lib = C.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Release the barcoded oligo-dT primer in each GEM</h3>
{panel([strand_row(bead), *annotation_rows(bead)], cls="small",
       caption="INFERRED — the v1 guide establishes P7–cell-barcode–Read-2–UMI–poly(T) order but does not print the oligo sequence.")}
<h3>(2) Reverse-transcribe and template-switch in the GEM</h3>
{panel(C.rt_scene().rows(), cls="small",
       caption="Cell barcode and UMI are installed next to the transcript's poly(A) end; the opposite end receives a template-switch handle.")}
<h3>(3) Amplify cDNA, shear, end-repair, A-tail and ligate the Read-1 adaptor</h3>
{panel([strand_row(C.final_library()), *annotation_rows(C.final_library())], cls="small",
       caption="Only the 3'-end fragment retains the bead-derived P7 side and becomes a P5–P7 library after sample-index PCR.")}
<h3>(4) INFERRED &mdash; final sequence-level reconstruction</h3>
{panel([*C.final_scene().rows(), *annotation_rows(lib)], cls="small",
       caption="INFERRED — dotted standard-adapter bases reconstruct the primary-source library architecture; no v1 oligo sequence is presented as published.")}'''


def sequencing() -> str:
    lib = C.final_library()
    rows = []
    for primer in C.SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Read 1: 98 cycles of cDNA; Index 1: 14-nt cell barcode; Index 2: 8-nt sample index; Read 2: 10-nt UMI. The v1 guide requires standard Illumina primers.')}
{sp.diagram(lib, C.SEQ_PRIMERS)}
{table(("Read", "Primer", "Primer site", "First bases"), rows)}'''


def main() -> None:
    page = "\n".join([head("10x Chromium Single Cell 3' GE v1 chemistry"), '<div class="wrap">',
        '<h1>10x Chromium Single Cell 3-prime Gene Expression v1</h1>',
        info('Defining source: <a href="https://doi.org/10.1038/ncomms14049">Zheng et al., <i>Nature Communications</i> (2017)</a>; protocol authority: 10x Genomics CG00026 Rev B.'),
        '<div class="caveat"><b>Sequence boundary.</b> The paper and v1 guide publish the library architecture and read layout, but no oligo sequences; reconstructed adapter bases are dotted.</div>',
        construction(), sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
