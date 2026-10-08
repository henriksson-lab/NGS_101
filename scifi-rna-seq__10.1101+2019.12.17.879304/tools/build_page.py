#!/usr/bin/env python3
"""Build the diagram-centric scifi-RNA-seq chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import nextera as nx
import scifi_rna as S
import seqprimers as sp
from chemdraw import annotation_rows, oligo, panel, strand_row
from page import head, info, table

OUT = HERE.parent / "scifi-rna-seq.html"


def oligos() -> str:
    rows = [
        oligo("Round-1 barcoded RT primer", S.rt_primer(), mods="/5Phos/"),
        oligo("10x scATAC gel-bead oligo", S.bead_oligo(), five="bead-|-5'-"),
        oligo("Bridge oligo", S.bridge_oligo(), three="-[ddC]-3'"),
        oligo("Template-switching oligo", [S.seg("SMART handle + GAAT", S.TSO[:-3], "tso"),
                                             S.seg("rGrGrG", "rGrGrG", placeholder=True)]),
        oligo("Partial-P5", [S.seg("partial P5", S.PARTIAL_P5, "p5")]),
        oligo("i7-only Tn5 adaptor", [S.seg("s7", nx.S7, "s7"),
                                       S.seg("mosaic end", nx.ME, "me")]),
    ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def construction() -> str:
    first = S.ligated_first_strand()
    lib = S.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Reverse-transcribe in 384 wells</h3>
{panel([strand_row(S.rt_primer()), *annotation_rows(S.rt_primer())], cls="small",
       caption="The phosphorylated RT primer writes an 8-nt UMI and 13-nt round-1 barcode next to the poly(A) junction.")}
<h3>(2) Ligate the overloaded droplet's bead barcode across a blocked bridge</h3>
{panel(S.ligation_scene().rows(), cls="small",
       caption="The bridge places bead s5 directly against the RT primer's Read-1 handle; Ampligase seals that nick.")}
{panel([strand_row(first), *annotation_rows(first)], cls="small",
       caption="One two-level-barcoded first strand after droplet thermoligation.")}
<h3>(3) Template-switch and enrich full cDNA in bulk</h3>
{panel(S.tso_scene().rows(), cls="small",
       caption="The TSO is added after emulsion break; its SMART handle supplies the opposite cDNA-PCR end.")}
<h3>(4) Tagment with i7-only Tn5 and select bead-end fragments</h3>
{panel([*S.final_scene().rows(), *annotation_rows(lib)], cls="small",
       caption="Partial-P5 selects the bead end and indexed P7 selects an s7-Tn5 end. Tn5/Tn5 and TSO/Tn5 fragments lack both PCR sites.")}'''


def sequencing() -> str:
    lib = S.final_library()
    rows = []
    for primer in S.SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('NovaSeq layout: Read 1 21 cycles; Index 1 8; Index 2 16; Read 2 78. Read 1 carries the plate barcode and UMI; Index 2 carries the bead barcode.')}
{sp.diagram(lib, S.SEQ_PRIMERS)}
{table(("Read", "Primer", "Primer site", "First bases"), rows)}
<h3>Read 1 layout</h3>
{table(("Cycles", "Content"), S.read1_layout(), scroll=False)}'''


def main() -> None:
    page = "\n".join([head("scifi-RNA-seq library chemistry"), '<div class="wrap">',
        '<h1>scifi-RNA-seq &mdash; combinatorial fluidic indexing</h1>',
        info('Defining source: <a href="https://doi.org/10.1101/2019.12.17.879304">Datlinger et al., bioRxiv (2019)</a>; expanded protocol: <a href="https://doi.org/10.1038/s41592-021-01153-z">Nature Methods (2021)</a>.'),
        oligos(), construction(), sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
