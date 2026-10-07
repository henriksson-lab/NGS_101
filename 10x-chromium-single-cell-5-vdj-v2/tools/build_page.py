#!/usr/bin/env python3
"""Build the sparse 10x Next GEM 5' V(D)J v2 schematic."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
AUDIT_TOOLS = HERE.parents[1] / "10x-chromium-single-cell-5-vdj" / "tools"
sys.path[:0] = [str(HERE), str(AUDIT_TOOLS), str(HERE.parents[1] / "lib")]

import seqprimers as sp
from chemdraw import Scene, Segment, annotation_rows, oligo, panel
from page import head, info, table
from vdj_common import VDJ_FORWARD
from vdj_v2 import CHEM as C

OUT = HERE.parent / "10x-5prime-vdj-v2.html"


def sequencing(lib) -> str:
    rows = []
    for primer in C.seq_primers:
        hit = sp.locate(lib, primer)
        if hit is None:
            raise ValueError(f"{primer.name} has no site")
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return ('<h2>Final library and read layout</h2>' +
            info('Paired-end, dual-index sequencing: Read 1 26 cycles, i7 10, i5 10, '
                 'Read 2 90. Read 1 reports the 16-base cell barcode and 10-base UMI.') +
            table(("Read", "Primer", "Primer site", "First bases"), rows) +
            table(("Read 1 cycles", "Content"),
                  [("1–16", "cell barcode"), ("17–26", "UMI")], scroll=False))


def render() -> str:
    amplicon = C.vdj_amplicon()
    lib = C.final_library()
    return "\n".join([
        head("10x Single Cell 5' V(D)J v2 chemistry"), '<div class="wrap">',
        '<h1>10x Chromium Next GEM Single Cell 5&#39; V(D)J v2</h1>',
        info('Primary source: <a href="https://www.10xgenomics.com/support">10x Genomics '
             'User Guide CG000331 Rev F</a>.'),
        '<h2>Key oligos</h2><seq>',
        oligo("gel-bead barcoded TSO", C.gel_bead_segments()[:-1], three="rGrGrG-3'"),
        oligo("poly(dT) RT primer", C.polydt_segments()),
        oligo("V(D)J universal forward primer", [Segment("", VDJ_FORWARD)]),
        oligo("TT i5 primer", C.index_p5_segments()),
        oligo("TT i7 primer", C.index_p7_segments()), '</seq>',
        '<h2>Library construction</h2>',
        '<h3>(1) Reverse-transcribe poly(A) RNA and template-switch onto the barcoded TSO</h3>',
        panel(C.capture_scene().rows(), cls="small"),
        panel(C.template_switch_scene().rows(), cls="long",
              caption="The gel-bead TSO installs the cell barcode, UMI and Read 1 handle."),
        '<h3>(2) Amplify cDNA, then enrich V(D)J with outer and inner constant-region pools</h3>',
        panel([*Scene.duplex(list(amplicon)).rows(), *annotation_rows(amplicon)], cls="long",
              caption="Both nested reactions use the universal forward primer; reverse "
                      "primer pools land in receptor constant regions."),
        '<h3>(3) Fragment, end-repair, A-tail and ligate the Read 2 adaptor</h3>',
        panel(C.adapter_scene().rows(), cls="small",
              caption="The half adaptor has a 12-base stem and 3'-T overhang."),
        '<h3>(4) Select the barcode-bearing fragment by dual-index PCR</h3>',
        panel([*Scene.duplex(list(lib)).rows(), *annotation_rows(lib)], cls="long",
              caption="Only fragments retaining the Read 1/barcode end acquire P5; Read 2 "
                      "starts at variable fragmentation points across V(D)J."),
        sequencing(lib), '</div>',
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
