#!/usr/bin/env python3
"""Build the diagram-centric paper-described PIP-seq page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import nextera as nx
import pip_seq as P
import seqprimers as sp
from chemdraw import annotation_rows, oligo, panel, strand_row
from page import head, info, table

OUT = HERE.parent / "pip-seq.html"


def oligos() -> str:
    rows = [
        oligo("Paper-described bead oligo", list(P.bead_oligo()), mods="/5Acryd/"),
        oligo("PIPS template-switching oligo", P.tso_oligo()),
        oligo("PIPS WTA primer", [P.seg("SMART handle", P.rt.SMART_HANDLE, "tso")]),
        oligo("PIPs P5 library primer", [P.seg("P5", P.il.P5, "p5"),
                                          P.seg("spacer", P.P5_SPACER),
                                          P.seg("SMART handle", P.rt.SMART_HANDLE, "tso"),
                                          P.seg("AC", "AC")], three="-*A*C-3'"),
        oligo("N7xx library-PCR primer", [P.seg("P7", P.il.P7, "p7"),
                                           P.seg("i7", "N" * 8, "cbc", placeholder=True),
                                           P.seg("s7", nx.S7, "s7")]),
    ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def construction() -> str:
    bead = P.bead_oligo()
    lib = P.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Build split-pool-barcoded hydrogel beads</h3>
{panel([strand_row(bead), *annotation_rows(bead)], cls="small",
       caption="INFERRED — the paper gives the bead's fixed ends, UMI and poly(T), but does not publish the intervening barcode blocks or splints.")}
<h3>(2) Vortex cells and beads into particle-templated droplets, then lyse</h3>
{panel(P.rt_scene().rows(), cls="small",
       caption="mRNA hybridises in droplets; after emulsion break, reverse transcription and template switching occur in bulk on the beads.")}
<h3>(3) Amplify handle-flanked cDNA with one WTA primer</h3>
{panel([strand_row(bead), *annotation_rows(bead)], cls="small",
       caption="The SMART handle at both ends supports single-primer whole-transcriptome amplification.")}
<h3>(4) Tagment and select the bead-end/s7 fragment</h3>
{panel([*P.final_scene().rows(), *annotation_rows(lib)], cls="small",
       caption="INFERRED — barcode/linker bases remain dotted. Published P5-SMART and named N70x primers determine the outer library structure.")}'''


def sequencing() -> str:
    lib = P.final_library()
    rows = []
    for primer in P.SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Read 1 reports the bead barcode cassette and 12-nt UMI; Index 1 reports the sample index; Read 2 reports cDNA. Exact barcode-cycle parsing is intentionally left to the research note because the paper does not publish the blocks.')}
{table(("Read", "Primer", "Primer site", "First bases"), rows)}'''


def main() -> None:
    page = "\n".join([head("PIP-seq library chemistry"), '<div class="wrap">',
        '<h1>PIP-seq &mdash; particle-templated single-cell RNA sequencing</h1>',
        info('Defining source: <a href="https://doi.org/10.1038/s41587-023-01685-z">Clark et al., <i>Nature Biotechnology</i> (2023)</a>.'),
        '<div class="caveat"><b>Barcode boundary.</b> The defining paper does not print the split-pool barcode blocks or splints; those regions remain dotted instead of adopting upstream reconstructions.</div>',
        oligos(), construction(), sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
