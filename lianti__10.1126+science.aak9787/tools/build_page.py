#!/usr/bin/env python3
"""Build the diagram-centric LIANTI chemistry page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import lianti as L
import seqprimers as sp
from chemdraw import (Scene, annotation_rows, oligo, panel,
                      protocol_graph_panels, strand_row)
from page import head, info

OUT = HERE.parent / "lianti.html"


def oligos() -> str:
    rows = [
        oligo("LIANTI hairpin transposon", L.transposon_segments()),
        oligo("Second-strand primer", L.second_strand_primer_segments()),
    ]
    return ("<h2>Key oligos</h2>"
            "<p><info>INFERRED throughout. The paper prints no oligo: the mosaic ends are "
            "the canonical Tn5 ME, the promoter the canonical T7 class III promoter, the "
            "spacer between them an unpublished placeholder of unknown length, and the UMI "
            "width a drawing choice. No library-prep adapter is shown because the paper "
            "names no kit.</info></p>"
            f"<seq>{''.join(rows)}</seq>")


def hairpin() -> str:
    return L.TRANSPOSON_TOPOLOGY.folded_panel(
        "INFERRED — the paper establishes a 19-bp double-stranded transposase binding "
        "site closed by a single-stranded T7-promoter loop (Fig. 1B), and the drawing "
        "fills that with the canonical mosaic end and T7 promoter. The spacer inside the "
        "loop is a placeholder: neither its bases nor its length is published, so the "
        "loop's drawn extent is not a measurement.")


def tagment_and_fill() -> str:
    lib = L.gap_filled_fragment()
    sc = Scene.duplex(list(lib))
    sc.mark("top", "T7 promoter", "T7 promoter →")
    sc.mark("top", "T7 promoter, opposite", "← opposite T7 promoter")
    return f"""<h3>(1) Load Tn5 with the hairpin transposon</h3>
{hairpin()}
<h3>(2) Tagment single-cell genomic DNA, then fill the 9-nt gaps</h3>
{panel([*sc.rows(), *annotation_rows(lib)], cls="small",
       caption="INFERRED — representative filled fragment. Gap fill makes a double-stranded "
               "T7 promoter at each end; the defining paper establishes the two inward-facing "
               "promoters, and every base shown is inferred or a placeholder.")}"""


def ivt() -> str:
    tx = L.transcript()
    return f"""<h3>(3) Amplify linearly by T7 in-vitro transcription</h3>
{panel([strand_row(tx), *annotation_rows(tx)], cls="small",
       caption="INFERRED — one RNA product, shown with T in place of U so the sequence "
               "relationships remain legible. Each genomic fragment is copied directly into "
               "many RNA molecules; copied molecules are not templates for another IVT round.")}"""


def self_primed_rt() -> str:
    amp = L.umi_amplicon()
    return f"""<h3>(4) Self-prime reverse transcription from the RNA 3' end</h3>
{L.transcript_topology().folded_panel(
    "INFERRED — nucleotide-resolved reconstruction of the self-priming described by "
    "the paper. The RNA 3′ mosaic end folds onto its complementary internal copy; "
    "N denotes the unpublished genomic sequence. The 5′ RNA body remains a straight "
    "tail, and the unpaired promoter/spacer interval is the flexible turn.")}
<h3>(5) Remove RNA and synthesize the second strand with the UMI primer</h3>
{panel(Scene.duplex(list(amp)).rows(), cls="small",
       caption="INFERRED — LIANTI amplicon orientation. The paper establishes that the "
               "second-strand primer carries a unique molecular barcode; its width and "
               "bases are not published.")}"""


def final_library() -> str:
    lib = L.final_library()
    drawing = panel([*Scene.duplex(list(lib)).rows(), *annotation_rows(lib)], cls="small",
                    caption="INFERRED throughout — one of two ligation orientations, drawn "
                            "with canonical single-index Illumina arms because the paper "
                            "names no library-prep kit, adapter or primer. The i7 width is a "
                            "placeholder.")
    return f"""<h3>(6) End repair, dA-tail, ligate an indexed Illumina adapter and PCR</h3>
{drawing}
<p><info>The opposite ligation orientation places the UMI-bearing LIANTI end next to P7.
Read&nbsp;1 or Read&nbsp;2 therefore encounters either genomic sequence immediately or
8-nt UMI + GGG + 19-nt mosaic end before the genome.</info></p>"""


def render() -> str:
    return "\n".join([
        head("LIANTI library chemistry"), '<div class="wrap" data-protocol-complete="true">',
        '<h1>LIANTI &mdash; linear whole-genome amplification by T7 transcription</h1>',
        info('Defining source: <a href="https://doi.org/10.1126/science.aak9787">Chen et al., '
             '<i>Science</i> (2017)</a>.'),
        '<div class="caveat"><b>Sequence boundary.</b> The accessible paper establishes the '
        'hairpin-transposon, gap-fill, IVT, self-primed RT and UMI architecture but prints no '
        'oligo, and its Supplementary Materials and Methods is behind the PMC download gate. '
        'Every dotted base is therefore inferred: the canonical Tn5 mosaic end and T7 promoter '
        'stand in for parts the paper names by role, the transposon spacer, UMI and index are '
        'placeholders of unpublished length, and no library-prep kit, adapter geometry, primer '
        'or read length is established.</div>',
        oligos(), '<h2>Library construction</h2>',
        protocol_graph_panels(L.protocol_graph(),cls="long"),
        sp.section(L.final_library(), L.SEQ_PRIMERS,
                   intro="INFERRED — the paper reports Illumina sequencing but neither a "
                         "library-prep kit nor a read length, so standard single-index TruSeq "
                         "primer sites are assumed and the 6-nt i7 width is a placeholder.",
                   required_roles=("Read 1", "Index 1 (i7)", "Read 2"),
                   read_lengths={"Index 1 (i7)": L.I7_NT},
                   show_read_content=True),
        '</div>'
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
