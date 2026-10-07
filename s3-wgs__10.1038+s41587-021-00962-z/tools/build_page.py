#!/usr/bin/env python3
"""Build the diagram-centric s3-WGS chemistry page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import illumina as il
import nextera as nx
import s3wgs as S
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, annotation_rows, oligo, panel, strand_row
from page import head, info

OUT = HERE.parent / "s3-wgs.html"


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def key_oligos() -> str:
    rows = [
        oligo("Indexed U-ME transposon strand", S.tn5_transfer_segments()),
        oligo("Mosaic-end bottom strand",
              [seg("mosaic end reverse complement", nx.ME_RC, "me")], mods="/5Phos/"),
        oligo("A14-LNA-ME adapter-switching oligo", S.switching_oligo_segments(),
              mods="ME bases 9, 11, 13, 15, 17 and 19 are LNA; /3InvdT/"),
        oligo("i7 PCR primer",
              [seg("P7", il.P7, "p7"), seg("i7", S.I7_PRIMER_INDEX, "cbc"),
               seg("TruSeq Read 2", il.TRUSEQ_READ2, "r2")]),
        oligo("i5 PCR primer",
              [seg("P5", il.P5, "p5"), seg("i5", S.I5_INDEX, "cbc"),
               seg("s5", nx.S5, "s5")]),
    ]
    return f"<h2>Key oligos</h2><seq>{''.join(rows)}</seq>"


def transposome() -> str:
    sc = Scene()
    sc.strand("transfer", S.tn5_transfer_segments(), label="transferred strand")
    sc.anneal("bottom", [seg("ME reverse complement", nx.ME_RC, "me")],
              to="transfer", pair=("ME reverse complement", "mosaic end"),
              label="Tn5-bound strand", mod5="p")
    sc.mark("transfer", "dU", "polymerase stop")
    return panel(sc.rows(), cls="small",
                 caption="Each of the 96 first-round barcodes is loaded on the same single "
                         "adaptor species. The internal dU sits immediately before the mosaic end.")


def tagmentation() -> str:
    left = S.tn5_transfer_segments()
    right = [seg("opposite mosaic end", nx.ME_RC, "me"),
             seg("opposite dU", "A", "w1"),
             seg("opposite Tn5 barcode", "b" * 8, "cbc", placeholder=True),
             seg("opposite partial Read 2", "r" * len(S.SBS12_PARTIAL), "r2", placeholder=True)]
    fragment = Construct([
        *left, seg("genomic DNA", "X" * 28, placeholder=True), *right
    ], name="symmetrically tagmented fragment")
    return panel([strand_row(fragment), *annotation_rows(fragment)], cls="small",
                 caption="Representative strand after symmetric tagmentation and gap closure. "
                         "Both fragment ends derive from the same indexed U-ME adaptor.")


def gap_fill_and_switch() -> str:
    pre = S.pre_switch_strand()
    sc = Scene()
    sc.strand("target", list(pre), label="gap-filled strand")
    sc.anneal("switch", S.switching_oligo_segments(), to="target",
              pair=("LNA mosaic end", "copied mosaic end"), label="A14-LNA-ME",
              mod3="/3InvdT/")
    sc.mark("target", "dU", "gap-fill polymerase stops here")
    sc.arrow("switch", "polymerase copies s5 onto the target strand")
    after = S.switched_strand()
    return f"""<h3>(2) Gap fill stops at dU</h3>
{panel([strand_row(pre), *annotation_rows(pre)], cls="small",
       caption="Nextera polymerase closes the tagmentation gap and copies the mosaic end, "
               "but does not copy the dU, barcode or partial Read-2 handle.")}
<h3>(3) Switch one end to s5 with A14-LNA-ME</h3>
{panel(sc.rows(), cls="small",
       caption="The blocked LNA oligo anneals to the newly copied mosaic end and acts only as "
               "a template. Ten denature/anneal/extend cycles maximize switching.")}
{panel([strand_row(after), *annotation_rows(after)], cls="small",
       caption="One converted strand: the original indexed end remains at the left and a copied "
               "s5 sequence now supplies the other PCR entry point.")}"""


def final_library() -> str:
    lib = S.final_library()
    return f"""<h3>(4) Uracil-tolerant indexed PCR completes the flow-cell adapters</h3>
{panel([*Scene.duplex(list(lib)).rows(), *annotation_rows(lib)], cls="small",
       caption="Final s3-WGS library, P5 to P7'. Q5U reads through the original dU, which "
               "therefore appears as A on the strand shown and T on its complement.")}
<p><info>Read&nbsp;1 enters genomic DNA directly. Read&nbsp;2 first reports the 8-nt Tn5 barcode,
then T + the 19-nt mosaic end, and reaches genomic DNA at cycle&nbsp;29.</info></p>"""


def render() -> str:
    return "\n".join([
        head("s3-WGS library chemistry"), '<div class="wrap">',
        '<h1>s3-WGS &mdash; symmetrical-strand combinatorial-indexed whole-genome sequencing</h1>',
        info('Defining source: <a href="https://doi.org/10.1038/s41587-021-00962-z">Mulqueen '
             'et al., <i>Nature Biotechnology</i> (2021)</a>.'),
        key_oligos(), '<h2>Library construction</h2>',
        '<h3>(1) Tagment with one indexed U-ME transposome per first-round well</h3>',
        transposome(), tagmentation(), gap_fill_and_switch(), final_library(),
        sp.section(S.final_library(), S.SEQ_PRIMERS,
                   intro="The protocol uses stock Nextera Read 1 / Index 2 and stock TruSeq "
                         "Read 2 / Index 1 primers."),
        '</div>'
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
