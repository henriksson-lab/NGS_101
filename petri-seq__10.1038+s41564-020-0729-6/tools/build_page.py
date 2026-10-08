#!/usr/bin/env python3
"""Build the diagram-centric PETRI-seq chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import petri_seq as P
import seqprimers as sp
from chemdraw import annotation_rows, junction_row, oligo, panel, strand_row
from page import head, info, table

OUT = HERE.parent / "petri-seq.html"


def oligos() -> str:
    rows = [
        oligo("Round-1 random-hexamer RT primer", P.round1_primer(), mods="/5Phos/"),
        oligo("Round-2 ligation oligo", P.round2_oligo(), mods="/5Phos/"),
        oligo("Round-2 linker", [P.seg("splint", P.LINKER_R2, placeholder=True)]),
        oligo("Round-3 ligation oligo", P.round3_oligo()),
        oligo("Round-3 linker", [P.seg("splint", P.LINKER_R3, placeholder=True)]),
    ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def construction() -> str:
    r1, r2, r3 = P.after_round1(), P.after_round2(), P.after_round3()
    lib = P.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Reverse-transcribe bacterial RNA in 96 wells</h3>
{panel([strand_row(r1), *annotation_rows(r1)], cls="small",
       caption="A phosphorylated random-hexamer primer installs the first 7-nt barcode in fixed, permeabilised cells.")}
<h3>(2) Ligate the second barcode across the round-2 linker</h3>
{panel([strand_row(r2), junction_row(r2, "round-2 junction", "round-1 constant"),
        *annotation_rows(r2)], cls="small",
       caption="The 16-nt splint spans the round-2/round-1 junction; blocking oligos sequester unused linker and arms.")}
<h3>(3) Ligate the UMI and third barcode across the round-3 linker</h3>
{panel([strand_row(r3), junction_row(r3, "round-3 junction", "round-2 constant"),
        *annotation_rows(r3)], cls="small",
       caption="The third oligo supplies the 7-nt UMI, barcode 3 and the future Read-1-side library handle.")}
<h3>(4) Lyse, make the second strand, tagment and selectively PCR</h3>
{panel([*P.final_scene().rows(), *annotation_rows(lib)], cls="small",
       caption="NEB i50x selects the round-3 handle and Nextera N7xx selects an s7 Tn5 end; s5-ended fragments do not amplify.")}'''


def sequencing() -> str:
    lib = P.final_library()
    rows = []
    for primer in P.SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('NextSeq layout: Read 1 58 cycles; Read 2 17; Index 1 and Index 2 8 cycles each. Read 1 carries the UMI and three cell barcodes.')}
{sp.diagram(lib, P.SEQ_PRIMERS)}
{table(("Read", "Primer", "Primer site", "First bases"), rows)}
<h3>Read 1 layout</h3>
{table(("Cycles", "Content"), P.read1_layout(), scroll=False)}'''


def main() -> None:
    page = "\n".join([head("PETRI-seq library chemistry"), '<div class="wrap">',
        '<h1>PETRI-seq &mdash; bacterial split-pool RNA indexing</h1>',
        info('Defining source: <a href="https://doi.org/10.1038/s41564-020-0729-6">Blattman et al., <i>Nature Microbiology</i> (2020)</a>.'),
        oligos(), construction(), sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
