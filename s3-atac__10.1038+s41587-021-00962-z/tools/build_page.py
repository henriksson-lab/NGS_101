#!/usr/bin/env python3
"""Build the diagram-centric s3-ATAC chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import s3atac as S
import seqprimers as sp
from chemdraw import Scene, annotation_rows, oligo, panel, strand_row
from page import head, info, table

OUT = HERE.parent / "s3-atac.html"


def key_oligos() -> str:
    rows = [
        oligo("Indexed U-ME transposon strand", S.tn5_transfer()),
        oligo("Mosaic-end bottom strand", S.tn5_bottom(), mods="/5Phos/"),
        oligo("A14-ME adapter-switching oligo", S.switching_oligo(),
              mods="ME positions " + ", ".join(map(str, S.LNA_ME_POSITIONS)) + " are LNA",
              three="/3InvdT/-3'"),
        oligo("i7 PCR primer", S.i7_pcr_primer()),
        oligo("i5 PCR primer", S.i5_pcr_primer()),
    ]
    return "<h2>Key oligos</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>"


def single_strand(con, caption: str) -> str:
    return panel([strand_row(con), *annotation_rows(con)], cls="small", caption=caption)


def construction() -> str:
    gap = S.gap_filled_strand()
    switched = S.switched_strand()
    lib = S.final_library()
    duplex = Scene.duplex(list(lib), label="library")
    return f'''<h2>Library construction</h2>
<h3>(1) Tagment accessible chromatin with one indexed U-ME adaptor species</h3>
{panel(S.transposome_scene().rows(), cls="small",
       caption="Each first-round well uses one of 96 Tn5 barcodes. The same adaptor is loaded on both Tn5 monomers.")}
{panel(S.tagmented_scene().rows(), cls="small",
       caption="Both strands of every accessible fragment receive the indexed adaptor. Tn5 leaves a 9-nt gap beside each transferred end.")}

<h3>(2) Gap-fill with uracil-intolerant NPM</h3>
{single_strand(gap, "NPM copies the opposite mosaic end and then stops before dU, leaving the Tn5 barcode and partial Read-2 handle single-copy.")}

<h3>(3) Switch the opposite end to s5</h3>
{panel(S.switching_scene().rows(), cls="small",
       caption="The LNA mosaic end anneals to the new ME copy. Its inverted-dT block prevents extension of the oligo itself.")}
{single_strand(switched, "NPM instead extends the target over the oligo's s5 template. Ten switching cycles convert both fragment strands.")}

<h3>(4) Add i5/P5 and i7/P7 by uracil-tolerant Q5U PCR</h3>
{panel([*duplex.rows(), *annotation_rows(lib)], cls="small",
       caption="Final s3-ATAC library, P5 to P7'. Q5U reads through dU; the strand shown carries A opposite that original uracil.")}
'''


def sequencing() -> str:
    lib = S.final_library()
    landing_rows = []
    for primer, hit in S.primer_landings(lib):
        landing_rows.append((html.escape(primer.role), html.escape(primer.name),
                             html.escape(", ".join(hit.covers)),
                             f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    layout_rows = [(cycles, content) for cycles, content in S.read2_layout()]
    return f'''<h2>Sequencing</h2>
{info('Paired-end 85 + 85 with 10-cycle index reads. Cell identity is the i5 + i7 + Tn5-barcode combination.')}
{sp.diagram(lib, S.SEQ_PRIMERS)}
{table(("Read", "Standard primer", "Primer site", "First bases"), landing_rows)}
<h3>Read 2 before genomic DNA</h3>
{table(("Cycles", "Content"), layout_rows, scroll=False)}
'''


def render() -> str:
    return "\n".join([
        head("s3-ATAC library chemistry"), '<div class="wrap">',
        '<h1>s3-ATAC &mdash; symmetrical-strand single-cell chromatin accessibility</h1>',
        info('Defining source: <a href="https://doi.org/10.1038/s41587-021-00962-z">Mulqueen et al., <i>Nature Biotechnology</i> (2021)</a>.'),
        key_oligos(), construction(), sequencing(), '</div>'
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
