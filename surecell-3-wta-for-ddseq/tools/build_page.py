#!/usr/bin/env python3
"""Build the commercial SureCell WTA 3' schematic."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import surecell as S
from chemdraw import Construct, Scene, Segment, annotation_rows, complement_segments, panel, strand_row
from page import head, info, table

OUT = HERE.parent / "surecell-wta3.html"


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def preamble() -> str:
    return f'''<div class="wrap">
<h1>SureCell WTA 3' for ddSEQ &mdash; commercial droplet 3' RNA sequencing</h1>
{info('Protocol source: Illumina <i>SureCell WTA 3\' Library Prep Reference Guide</i>, document #1000000021452 v01 (June 2017).')}
<div class="caveat"><b>Proprietary oligos.</b> The authoritative Illumina guide names the bead mix, TPP1, N7xx adapters and custom Sequencing Primer SP but prints no sequence or barcode anatomy. Every dotted bracketed region is therefore a non-sequence placeholder, not an upstream reconstruction; spacing is not to scale.</div>
'''


def bead_panel() -> str:
    bead = S.bead_primer()
    return panel([strand_row(bead), *annotation_rows(bead)], cls="small",
                 caption="INFERRED — minimal capture-oligo roles consistent with the assay workflow. The proprietary barcode region and oligo-dT length are not published.")


def rt_scene() -> Scene:
    mrna = [seg("mRNA", "XXXXXXXX...XXXXXXXX", placeholder=True), seg("poly(A)", "A" * 12)]
    primer = [S.inferred("3' Barcode Mix oligo", "[proprietary barcode region]", "cbc"),
              S.inferred_poly_t("dT")]
    sc = Scene(); sc.strand("mRNA", mrna); sc.anneal("bead primer", primer, to="mRNA", pair=("dT", "poly(A)")); sc.arrow("bead primer", "reverse transcriptase")
    return sc


def second_strand() -> str:
    first = S.first_strand()
    sc = S.schematic_duplex(first)
    return panel([*sc.rows(), *annotation_rows(first)],
                 caption="INFERRED — replacement second strand after the guide's RNA-removal / second-strand reaction. Proprietary sequence regions remain dotted.")


def tagmented() -> str:
    con = Construct([
        S.inferred("bead end", "[bead-derived end]", "r1"),
        seg("cDNA fragment", "XXXXXXXX...XXXXXXXX", placeholder=True),
        S.inferred("transposome end", "[SureCell transposome end]", "me"),
    ])
    return panel([strand_row(con), *annotation_rows(con)],
                 caption="INFERRED — bead-end fragment after Nextera SureCell tagmentation. The guide does not disclose the transposome adapter load or sequence.")


def final_panel() -> str:
    lib = S.selected_fragment()
    sc = S.schematic_duplex(lib)
    return panel([*sc.rows(), *annotation_rows(lib)],
                 caption="INFERRED — qualitative final library after TPP1 + N7xx PCR. Bracketed regions encode reagent roles only, not bases or lengths.")


def steps() -> str:
    return f'''<h2>Library generation</h2>
<h3>(1) Co-encapsulate cells, RT mix and barcode beads</h3>
{bead_panel()}
<h3>(2) Lyse and reverse-transcribe inside droplets</h3>
{panel(rt_scene().rows(), caption="INFERRED — oligo-dT priming installs a 3' Barcode Mix-derived region on first-strand cDNA; the guide does not disclose its molecular layout.")}
<h3>(3) Break droplets and remove unbound barcode reagent</h3>
{info('Purification beads retain the pooled first-strand products while free barcode reagent is removed.')}
<h3>(4) Remove the RNA template and synthesize the replacement strand</h3>
{second_strand()}
<h3>(5) Tagment double-stranded cDNA with the Nextera SureCell transposome</h3>
{tagmented()}
<h3>(6) Select bead-end fragments by PCR with TPP1 and one N7xx adapter</h3>
{final_panel()}
'''


def sequencing() -> str:
    boundaries = S.read_boundaries()
    rows = [(read, boundary) for read, boundary in boundaries]
    return f'''<h2>Supported read boundaries</h2>
{table(("Read", "Starts on"), rows)}
{info('The guide requires the supplied custom Sequencing Primer SP for Read 1 and names the N7xx sample index. It does not publish read lengths or oligo sequences, so no base-resolved layout is asserted here.')}
</div>'''


def main() -> None:
    OUT.write_text("\n".join([head("SureCell WTA 3' library chemistry"), preamble(),
                              steps(), sequencing()]), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
