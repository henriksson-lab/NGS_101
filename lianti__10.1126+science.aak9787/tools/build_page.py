#!/usr/bin/env python3
"""Build the diagram-centric LIANTI chemistry page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import illumina as il
import lianti as L
import seqprimers as sp
from chemdraw import Scene, Segment, annotation_rows, oligo, panel, strand_row
from page import head, info

OUT = HERE.parent / "lianti.html"


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def oligos() -> str:
    rows = [
        oligo("LIANTI hairpin transposon", L.transposon_segments(), mods="/5Phos/"),
        oligo("Second-strand primer", L.second_strand_primer_segments()),
        oligo("NEBNext hairpin adaptor",
              [seg("Read 2 arm", il.NEBNEXT_ARM_READ2, "r2"),
               seg("dU", "U", placeholder=True),
               seg("Read 1 arm + T", il.NEBNEXT_ARM_READ1, "r1")],
              mods="/5Phos/"),
    ]
    return f"<h2>Key oligos</h2><seq>{''.join(rows)}</seq>"


def hairpin() -> str:
    stem_5 = [seg("5' ME reverse complement", L.TRANSPOSON[:19], "me", inferred=True)]
    stem_3 = [seg("3' mosaic end", L.TRANSPOSON[-19:], "me", inferred=True)]
    sc = Scene()
    sc.strand("5' arm", stem_5, label="same oligo", mod5="p")
    sc.anneal("3' arm", stem_3, to="5' arm",
              pair=("3' mosaic end", "5' ME reverse complement"), label="same oligo")
    sc.mark("5' arm", "5' ME reverse complement", "19-bp Tn5 mosaic-end stem")
    return panel(sc.rows(), cls="small",
                 caption="INFERRED — fold implied by the main paper and reconstructed with the "
                         "secondary-source sequence. The 30-nt loop between these two arms "
                         "contains the T7 promoter; both arms belong to one covalent strand.")


def tagment_and_fill() -> str:
    lib = L.gap_filled_fragment()
    sc = Scene.duplex(list(lib))
    sc.mark("top", "T7 promoter", "T7 promoter →")
    sc.mark("top", "T7 promoter, opposite", "← opposite T7 promoter")
    return f"""<h3>(1) Load Tn5 with the hairpin transposon</h3>
{hairpin()}
<h3>(2) Tagment single-cell genomic DNA, then fill the 9-nt gaps</h3>
{panel([*sc.rows(), *annotation_rows(lib)], cls="small",
       caption="INFERRED — representative filled fragment assembled from the secondary-source "
               "oligo. Gap fill makes a double-stranded T7 promoter at each end; the defining "
               "paper establishes the two inward-facing promoters.")}"""


def ivt() -> str:
    tx = L.transcript()
    return f"""<h3>(3) Amplify linearly by T7 in-vitro transcription</h3>
{panel([strand_row(tx), *annotation_rows(tx)], cls="small",
       caption="INFERRED — one RNA product, shown with T in place of U so the sequence "
               "relationships remain legible. Each genomic fragment is copied directly into "
               "many RNA molecules; copied molecules are not templates for another IVT round.")}"""


def self_primed_rt() -> str:
    body = [seg("RNA body", "GGGXXXXXXXX...XXXXXXXX", placeholder=True),
            seg("self-primer site", L.TRANSPOSON[:19], "me", inferred=True)]
    end = [seg("loop", "L" * 30, placeholder=True, inferred=True),
           seg("3' self-primer", L.TRANSPOSON[-19:], "me", inferred=True)]
    sc = Scene()
    sc.strand("RNA body", body, label="RNA 5' body")
    sc.anneal("same RNA 3' end", end, to="RNA body",
              pair=("3' self-primer", "self-primer site"), label="RNA 3' end",
              unpaired=("loop",))
    sc.arrow("same RNA 3' end", "reverse transcriptase copies back across the genomic RNA")
    sc.mark("same RNA 3' end", "3' self-primer", "intramolecular 19-bp stem")
    amp = L.umi_amplicon()
    return f"""<h3>(4) Self-prime reverse transcription from the RNA 3' end</h3>
{panel(sc.rows(), cls="small",
       caption="INFERRED — structural reconstruction of the self-priming described by the "
               "paper. The transcript's terminal mosaic end folds onto its complementary "
               "internal copy; the same RNA is shown on two aligned rows.")}
<h3>(5) Remove RNA and synthesize the second strand with the UMI primer</h3>
{panel(Scene.duplex(list(amp)).rows(), cls="small",
       caption="INFERRED — LIANTI amplicon orientation reconstructed from the "
               "secondary-source primer. The paper independently establishes that the "
               "second-strand primer carries a unique molecular barcode.")}"""


def final_library() -> str:
    lib = L.final_library()
    drawing = panel([*Scene.duplex(list(lib)).rows(), *annotation_rows(lib)], cls="small",
                    caption="INFERRED — one of two ligation orientations. Dotted LIANTI "
                            "regions and the 6-nt i7 assignment come from the upstream page; "
                            "solid Illumina/NEBNext regions use canonical vendor sequences.")
    return f"""<h3>(6) End repair, dA-tail, ligate the NEBNext hairpin adaptor, open it and PCR</h3>
{drawing}
<p><info>The opposite ligation orientation places the UMI-bearing LIANTI end next to P7.
Read&nbsp;1 or Read&nbsp;2 therefore encounters either genomic sequence immediately or
8-nt UMI + GGG + 19-nt mosaic end before the genome.</info></p>"""


def render() -> str:
    return "\n".join([
        head("LIANTI library chemistry"), '<div class="wrap">',
        '<h1>LIANTI &mdash; linear whole-genome amplification by T7 transcription</h1>',
        info('Defining source: <a href="https://doi.org/10.1126/science.aak9787">Chen et al., '
             '<i>Science</i> (2017)</a>.'),
        '<div class="caveat"><b>Sequence source.</b> The accessible paper establishes the '
        'hairpin-transposon, gap-fill, IVT, self-primed RT and UMI architecture, but does not '
        'print the oligos. Dotted LIANTI-specific bases and steps labelled <b>INFERRED</b> use '
        'the upstream scg_lib_structs reconstruction because the supplement is unavailable.</div>',
        oligos(), '<h2>Library construction</h2>', tagment_and_fill(), ivt(),
        self_primed_rt(), final_library(),
        sp.section(L.final_library(), L.SEQ_PRIMERS,
                   intro="Standard TruSeq primers address the NEBNext-built library. The 6-nt "
                         "i7 placement follows the secondary upstream reconstruction."),
        '</div>'
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
