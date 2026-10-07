#!/usr/bin/env python3
"""Build the diagram-centric scifi-ATAC-seq chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import nextera as nx
import scifi_atac as S
import seqprimers as sp
from chemdraw import annotation_rows, oligo, panel, strand_row
from page import head, info, table

OUT = HERE.parent / "scifi-atac-seq.html"


def oligos() -> str:
    rows = [
        oligo("Tn5 A adapter", S.adapter_a()),
        oligo("Tn5 B adapter", S.adapter_b()),
        oligo("Shared ME bottom", [S.seg("ME reverse complement", nx.ME_RC, "me")], mods="/5Phos/"),
        oligo("10x scATAC v1.1 gel-bead oligo", S.bead_oligo(), five="bead-|-5'-"),
    ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def construction() -> str:
    frag = S.indexed_fragment()
    lib = S.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Tagment nuclei in 96 wells with one A/B barcode pair</h3>
{panel([strand_row(frag), *annotation_rows(frag)], cls="small",
       caption="Only fragments with one A and one B end can later receive P5 from the bead and P7 from sample-index PCR.")}
<h3>(2) Pool and overload one unmodified 10x scATAC v1.1 channel</h3>
{panel([strand_row(S.bead_oligo()), *annotation_rows(S.bead_oligo())], cls="small",
       caption="Gap fill precedes linear extension. The bead's s5 end primes A adapters and adds the 16-nt GEM barcode plus P5.")}
<h3>(3) Sample-index PCR selects the B end</h3>
{panel([*S.final_scene().rows(), *annotation_rows(lib)], cls="small",
       caption="Final library combines two inline 5-nt Tn5 barcodes, the 16-nt GEM barcode and the i7 sample index.")}'''


def sequencing() -> str:
    lib = S.final_library()
    rows = []
    for primer in S.SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Four custom primers are used. Read 1 and Read 2 each begin with a 5-nt Tn5 barcode then the 19-nt ME; Index 1 reports sample index and Index 2 the 16-nt GEM barcode.')}
{table(("Read", "Primer", "Primer site", "First bases"), rows)}'''


def main() -> None:
    page = "\n".join([head("scifi-ATAC-seq library chemistry"), '<div class="wrap">',
        '<h1>scifi-ATAC-seq &mdash; plate pre-indexing plus overloaded 10x scATAC</h1>',
        info('Defining source: <a href="https://doi.org/10.1101/2023.09.17.558155">Wang et al., bioRxiv (2023)</a>; droplet chemistry follows the cited 10x scATAC v1.1 guide.'),
        oligos(), construction(), sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
