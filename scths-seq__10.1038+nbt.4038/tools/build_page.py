#!/usr/bin/env python3
"""Build the diagram-centric scTHS-seq chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import nextera as nx
import scths as S
import seqprimers as sp
from chemdraw import annotation_rows, oligo, panel, strand_row
from page import head, info, table

OUT = HERE.parent / "scths-seq.html"


def oligos() -> str:
    rows = [
        oligo("Barcoded r5/T7 transposon", S.r5_transposon()),
        oligo("Shared phosphorylated ME bottom", [S.seg("ME reverse complement", nx.ME_RC, "me")], mods="/5Phos/"),
        oligo("s7-only second transposon", S.s7_transposon()),
        oligo("scT7 S502 PCR primer", [S.seg("P5", S.il.P5, "p5"),
                                       S.seg("i5", S.I5_OLIGO_INDEX, "cbc"),
                                       S.seg("connector", S.CONNECTOR, "t7")]),
        oligo("N701 PCR primer", [S.seg("P7", S.il.P7, "p7"),
                                   S.seg("i7", S.I7_OLIGO_INDEX, "cbc"),
                                   S.seg("s7", nx.S7, "s7")]),
    ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def construction() -> str:
    ivt = S.ivt_template()
    lib = S.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Tagment nuclei in 384 wells with a barcoded T7/r5/s5 homodimer</h3>
{panel(S.first_tagmented_scene().rows(), cls="small",
       caption="Each accessible-site insertion carries a T7 promoter, 6-nt round-1 barcode and s5/ME sequence.")}
<h3>(2) Sort single nuclei, fill gaps and transcribe each insertion end</h3>
{panel([strand_row(ivt), *annotation_rows(ivt)], cls="small",
       caption="Gap fill makes the promoter double-stranded. T7 IVT starts in the connector and amplifies each insertion independently.")}
<h3>(3) Random-prime RNA, make dsDNA from the connector, then tagment with s7-only Tn5</h3>
{panel([strand_row(ivt), *annotation_rows(ivt)], cls="small",
       caption="The connector primer restores dsDNA whose original accessible-site end remains identifiable by the r5 barcode.")}
<h3>(4) Add i5/P5 at the connector end and i7/P7 at the new s7 end</h3>
{panel([*S.final_scene().rows(), *annotation_rows(lib)], cls="small",
       caption="Final single-end library. Each molecule contains one original r5-tagged insertion and one downstream s7 cut.")}'''


def sequencing() -> str:
    lib = S.final_library()
    rows = []
    for primer in S.SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('HiSeq layout: Read 1 50 cycles; Index 1 8; Index 2 32; no Read 2. The long Index 2 read covers i5, connector, r5 cell barcode and three s5 bases.')}
{sp.diagram(lib, S.SEQ_PRIMERS)}
{table(("Read", "Primer", "Primer site", "First bases"), rows)}'''


def main() -> None:
    page = "\n".join([head("scTHS-seq library chemistry"), '<div class="wrap">',
        '<h1>scTHS-seq &mdash; T7-amplified single-cell accessibility</h1>',
        info('Defining source: <a href="https://doi.org/10.1038/nbt.4038">Lake et al., <i>Nature Biotechnology</i> (2018)</a>.'),
        oligos(), construction(), sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
