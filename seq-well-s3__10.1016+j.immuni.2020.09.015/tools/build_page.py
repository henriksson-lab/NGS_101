#!/usr/bin/env python3
"""Build the diagram-centric Seq-Well S3 chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import nextera as nx
import seqprimers as sp
import seqwell_s3 as S
from chemdraw import annotation_rows, oligo, panel, strand_row
from page import head, info, table

OUT = HERE.parent / "seq-well-s3.html"


def oligos() -> str:
    rows = [
        oligo("Seq-Well bead oligo", list(S.bead_oligo()), five="bead-|-5'-"),
        oligo("Template-switching oligo", S.tso_oligo()),
        oligo("S3 randomer", S.randomer()),
        oligo("SMART PCR primer", [S.seg("SMART handle", S.rt.SMART_HANDLE, "tso")]),
        oligo("P5-SMART hybrid", [S.seg("P5", S.il.P5, "p5"),
                                   S.seg("spacer", S.P5_SPACER),
                                   S.seg("SMART handle", S.rt.SMART_HANDLE, "tso"),
                                   S.seg("AC", S.BEAD_SUFFIX)],
              three="-*A*C-3'"),
        oligo("N7xx library-PCR primer", [S.seg("P7", S.il.P7, "p7"),
                                           S.seg("i7", "N" * 8, "cbc", placeholder=True),
                                           S.seg("s7", nx.S7, "s7")]),
    ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def construction() -> str:
    bead = S.bead_oligo()
    lib = S.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Capture mRNA in nanowells and reverse-transcribe on recovered beads</h3>
{panel(S.rt_scene().rows(), cls="small",
       caption="The bead writes the 12-nt cell barcode and 8-nt UMI next to the transcript's poly(A) end; the TSO is present but need not succeed.")}
<h3>(2) Strip RNA and make a handle-tagged second strand</h3>
{panel(S.random_second_strand_scene().rows(), cls="small",
       caption="After NaOH removes RNA, the S3 randomer lands on first-strand cDNA and Klenow exo-minus copies toward the bead.")}
<h3>(3) Amplify with the single SMART primer</h3>
{panel([strand_row(bead), *annotation_rows(bead)], cls="small",
       caption="Both template switching and random second-strand synthesis lead to SMART-handle-flanked WTA products.")}
<h3>(4) Tagment and select the bead-end/s7 fragment</h3>
{panel([*S.final_scene().rows(), *annotation_rows(lib)], cls="small",
       caption="The P5-SMART hybrid selects the bead end; N7xx selects an s7 Tn5 end. The S3 randomer changes capture yield, not the final read layout.")}'''


def sequencing() -> str:
    lib = S.final_library()
    rows = []
    for primer in S.SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Read 1 reports the 12-nt cell barcode and 8-nt UMI; Index 1 reports the N700 sample index; Read 2 reports cDNA from the Tn5 end.')}
{table(("Read", "Primer", "Primer site", "First bases"), rows)}'''


def main() -> None:
    page = "\n".join([head("Seq-Well S3 library chemistry"), '<div class="wrap">',
        '<h1>Seq-Well S3 &mdash; randomly primed second-strand capture</h1>',
        info('Defining source: <a href="https://doi.org/10.1016/j.immuni.2020.09.015">Hughes et al., <i>Immunity</i> (2020)</a>.'),
        oligos(), construction(), sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
