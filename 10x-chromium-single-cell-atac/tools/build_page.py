#!/usr/bin/env python3
"""Build the 10x Chromium Single Cell ATAC v1 chemistry page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import chromium_atac as C
import nextera as nx
import seqprimers as sp
from chemdraw import annotation_rows, oligo, panel, strand_row
from page import head, info

OUT = HERE.parent / "10x-chromium-single-cell-atac.html"


def oligos() -> str:
    rows = [
        oligo("Read 1N transposon top", [C.seg("s5", nx.S5, "s5"), C.seg("ME", nx.ME, "me")]),
        oligo("Read 2N transposon top", [C.seg("s7", nx.S7, "s7"), C.seg("ME", nx.ME, "me")]),
        oligo("Gel-bead oligo PN-2000132", C.bead_oligo(), five="bead-|-5'-"),
        oligo("SI-PCR Primer B", [C.seg("partial P5", C.il.P5[:22], "p5")]),
        oligo("i7 sample-index primer", [C.seg("P7", C.il.P7, "p7"),
                                          C.seg("i7", "N" * 8, "cbc", placeholder=True),
                                          C.seg("s7", nx.S7, "s7")]),
    ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def construction() -> str:
    template = C.gap_filled_template()
    lib = C.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Tagment accessible chromatin in intact nuclei</h3>
{panel(C.transposed_scene().rows(), cls="small",
       caption="Bulk Nextera transposition creates s5/s7 fragments with the characteristic two 9-nt gaps.")}
<h3>(2) Partition nuclei with dissolving barcoded gel beads</h3>
{panel([strand_row(template), *annotation_rows(template)], cls="small",
       caption="The first 72-degree step fills the Tn5 gaps before denaturation.")}
{panel(C.capture_scene().rows(), cls="small",
       caption="The bead oligo's 3' s5 end primes linearly from s5' and attaches P5 plus the 16-nt cell barcode.")}
<h3>(3) Sample-index PCR completes the s7 end</h3>
{panel([*C.final_scene().rows(), *annotation_rows(lib)], cls="small",
       caption="Final dual-index Nextera library: cell barcode in the i5 position and sample index in i7.")}'''


def sequencing() -> str:
    lib = C.final_library()
    return sp.section(
        lib, C.SEQ_PRIMERS, heading="Sequencing",
        intro=("Read 1 and Read 2 cover accessible genomic DNA; Index 1 carries the "
               "sample index and Index 2 carries the 10x cell barcode."),
        read_lengths=C.READ_LENGTHS, show_read_content=True)


def main() -> None:
    page = "\n".join([head("10x Chromium Single Cell ATAC chemistry"), '<div class="wrap">',
        '<h1>10x Chromium Single Cell ATAC v1</h1>',
        info('Authoritative source: 10x Genomics <a href="https://www.10xgenomics.com/support/single-cell-atac/documentation/steps/library-prep/chromium-single-cell-atac-reagent-kits-user-guide-v-1-chemistry">CG000168 Rev A user guide</a>.'),
        oligos(), construction(), sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
