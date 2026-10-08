#!/usr/bin/env python3
"""Build the diagram-centric Quartz-Seq2 chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import quartz_seq2 as Q
import seqprimers as sp
from chemdraw import Scene, Segment, annotation_rows, oligo, panel, strand_row
from page import head, info, table

OUT = HERE.parent / "quartz-seq2.html"


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


def oligos() -> str:
    rows = [
        oligo("1536-well RT primer", Q.rt_primer()),
        oligo("Tagging primer", Q.tagging_primer()),
        oligo("gM suppression primer", [seg("G", "G"), seg("M handle", Q.M, "tso")]),
        oligo("rYshapeP5 short strand", [seg("stem", Q.RYSHAPE_P5[:12], "r1"),
                                          seg("unpaired arm", Q.RYSHAPE_P5[12:], "r1")]),
        oligo("rYshapeP7LT06 long strand",
              [seg("P7", Q.RYSHAPE_P7_LT06[:24], "p7"),
               seg("i7", Q.LT06_OLIGO_INDEX, "cbc"),
               seg("TruSeq Read 2", Q.RYSHAPE_P7_LT06[-34:], "r2")]),
        oligo("P5-gMac selective PCR primer",
              [seg("truncated P5", Q.P5_GMAC[:28], "p5"), seg("TT", "TT"),
               seg("gM + AC", Q.P5_GMAC[30:], "tso")],
              mods="two 3'-terminal phosphorothioates"),
        oligo("Read1DropQuartz sequencing primer", [seg("custom R1", Q.READ1_DROPQUARTZ, "r1")]),
    ]
    return "<h2>Key oligos</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>"


def construction() -> str:
    tailed = Q.tailed_first_strand()
    lib = Q.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Reverse-transcribe with a well-specific barcode and UMI</h3>
{panel(Q.mrna_rt_scene().rows(), cls="small",
       caption="The RT primer writes the cell barcode and UMI next to the transcript's poly(A) end. Wells can be pooled after this step.")}

<h3>(2) Add a poly(A) tail and synthesize the second strand</h3>
{panel([strand_row(tailed), *annotation_rows(tailed)], cls="small",
       caption="Terminal transferase creates the tagging-primer site at the cDNA 3' end.")}
{panel(Q.tagging_scene().rows(), cls="small",
       caption="The M-tagging primer starts the second strand; G+M then suppression-amplifies the pooled cDNA.")}

<h3>(3) Shear, end-repair, A-tail, and ligate the truncated Y adapter</h3>
{panel(Q.truncated_adapter_scene().rows(), cls="small",
       caption="The P7 strand supplies the 3'-T overhang. The short strand makes a 12-bp stem and leaves an 8-nt arm unpaired.")}

<h3>(4) Select the RT-primer end with P5-gMac and complete the library with TPC2</h3>
{panel([*Scene.duplex(list(lib), label="library").rows(), *annotation_rows(lib)], cls="small",
       caption="Final Quartz-Seq2 v3.2 library, P5 to P7'. P5-gMac ends in M+AC, so only the cell-barcode/UMI end is exponentially amplified.")}
'''


def sequencing() -> str:
    lib = Q.final_library()
    rows = []
    for primer, hit in Q.primer_landings():
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Quartz-Seq2 v3.2 uses Read 1: 23 cycles, Index 1: 6, and Read 2: 63. Read 1 reports the 15-base cell barcode followed by the 8-base UMI; Read 2 reports cDNA.')}
{sp.diagram(lib, Q.SEQ_PRIMERS)}
{table(("Read", "Primer", "Primer site", "First bases"), rows)}
<h3>Read 1 barcode layout</h3>
{table(("Cycles", "Content"), Q.read1_layout(), scroll=False)}
'''


def main() -> None:
    page = "\n".join([head("Quartz-Seq2 library chemistry"), '<div class="wrap">',
        '<h1>Quartz-Seq2 &mdash; barcoded 3-prime single-cell RNA sequencing</h1>',
        info('Defining source: <a href="https://doi.org/10.1186/s13059-018-1407-3">Sasagawa et al., <i>Genome Biology</i> (2018)</a>.'),
        oligos(), construction(), sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
