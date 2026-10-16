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
from chemdraw import Scene, annotation_rows, oligo, panel, protocol_graph_panels, strand_row
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
    return '<h2>Library construction</h2>\n'+protocol_graph_panels(S.protocol_graph(),cls="long")


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
{sp.read_content_table(lib, S.SEQ_PRIMERS,
    {"Read 1":85,"Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":85})}
'''


def render() -> str:
    return "\n".join([
        head("s3-ATAC library chemistry"), '<div class="wrap" data-protocol-complete="true">',
        '<h1>s3-ATAC &mdash; symmetrical-strand single-cell chromatin accessibility</h1>',
        info('Defining source: <a href="https://doi.org/10.1038/s41587-021-00962-z">Mulqueen et al., <i>Nature Biotechnology</i> (2021)</a>.'),
        info('Source discrepancy: published Supplementary Table 5 is headed “Oligonucleotides used for PCR, i7 Molecule side,” the same wording as Table 4. Table 5 rows are named PCR_A..H_i5_A..H and carry P5/i5/s5 arms, while Table 4 rows carry P7/i7/Read 2 arms. The model therefore treats Table 5 as the i5 primer set. The authors have not confirmed whether its heading was a copy error.'),
        key_oligos(), construction(), sequencing(), '</div>'
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
