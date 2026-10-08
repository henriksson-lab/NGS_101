#!/usr/bin/env python3
"""Build the molecule-focused Drop-ChIP schematic."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import dropchip as D
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, annotation_rows, oligo, panel, strand_row
from page import head, info, table

OUT = HERE.parent / "drop-chip.html"


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def oligos() -> str:
    return f"""<h2>Key oligos</h2><seq>
{oligo("representative barcode adaptor, forward strand", list(D.barcode_adaptor()))}
{oligo("SC-PCR1", [seg("", D.SC_PCR1)])}
{oligo("SC-PCR2", [seg("", D.SC_PCR2)])}
</seq>
{panel(D.adaptor_scene().rows(), cls="long",
       caption="The published 60-nt strand and its derived complement. The same 8-nt "
               "identity is encoded at both ends, with the second copy reversed.")}"""


def ligation() -> str:
    molecule = D.blunt_ligated_fragment()
    sc = Scene.duplex(list(molecule))
    sc.junction("top", "left right blunt end", "nucleosomal DNA")
    sc.junction("top", "nucleosomal DNA", "right left blunt end")
    return f"""<h2>Library construction</h2>
<h3>(1) End-repair nucleosomal DNA and blunt-ligate barcode adaptors</h3>
{panel(sc.rows(), cls="long",
       caption="Representative ligation orientation. A double-stranded barcode adaptor is "
               "joined to both blunt fragment ends; either adaptor may ligate in either "
               "orientation.")}"""


def paci_and_pcr() -> str:
    paci = D.paci_product()
    return f"""<h3>(2) PacI removes adaptor concatemers and the outer adaptor halves</h3>
{panel([*Scene.duplex(list(paci)).rows(), *annotation_rows(paci)], cls="long",
       caption="PacI cuts the central TTAATTAA site. Each retained half still contains one "
               "barcode and an SC-PCR priming site.")}
<h3>(3) SC-PCR amplifies the barcode-labelled fragment</h3>
{panel(D.sc_pcr_scene().rows(), cls="long",
       caption="The two published 17-nt primers land on opposite strands of the PacI product.")}"""


def bcivi_and_boundary() -> str:
    scene, _ = D.bcivi_product()
    final = D.inferred_final_library()
    return f"""<h3>(4) BciVI removes the SC-PCR ends and leaves 3'-A overhangs</h3>
{panel(scene.rows(), cls="long",
       caption="Type-IIS cutting retains, at each insert end, an 11-nt constant, the 8-nt "
               "cell barcode and the TTAA PacI half-site.")}
<h3>(5) Ligate Illumina adaptors and perform library PCR</h3>
{panel([strand_row(final), *annotation_rows(final)], cls="long",
       caption="INFERRED — final architecture only. The paper does not publish the Illumina "
               "adaptor or library-primer sequences, so the dotted outer arms stop at that "
               "authoritative boundary.")}"""


def reads() -> str:
    final = D.inferred_final_library()
    rows = [
        ("Read 1", "8-nt cell barcode at bases 1–8, then TTAA and genomic DNA"),
        ("Read 2", "11-nt constant, 8-nt cell barcode at bases 12–19, then TTAA and genomic DNA"),
        ("Index", "sample index distinguishes the pooled library; sequence not published"),
    ]
    r1 = Construct([seg("cell barcode", "N" * 8, "cbc", placeholder=True),
                    seg("PacI half", "TTAA"),
                    seg("genomic DNA", "XXXXXXXX...", placeholder=True)], name="Read 1")
    r2 = Construct([seg("11 constant cycles", "XXXXXXXXXXX", placeholder=True),
                    seg("cell barcode", "N" * 8, "cbc", placeholder=True),
                    seg("PacI half", "TTAA"),
                    seg("genomic DNA", "XXXXXXXX...", placeholder=True)], name="Read 2")
    return f"""<h2>Read layout</h2>
{sp.unavailable_diagram(final, "the paper does not publish the Illumina library-primer sequences")}
{panel([strand_row(r1, prefix="Read 1  ", suffix=""),
        strand_row(r2, prefix="Read 2  ", suffix="")],
       cls="small", caption="Published barcode coordinates. The sequencing run uses 11 "
                            "initial dark cycles before collecting base calls.")}
{table(("Read", "Published layout"), rows)}"""


def render() -> str:
    return "\n".join([
        head("Drop-ChIP library chemistry"), '<div class="wrap">',
        '<h1>Drop-ChIP &mdash; barcoded nucleosomal DNA for single-cell ChIP-seq</h1>',
        info('Defining source: <a href="https://doi.org/10.1038/nbt.3383">Rotem et al., '
             '<i>Nature Biotechnology</i> (2015)</a>.'),
        '<div class="caveat"><b>Downstream adaptor boundary.</b> The paper publishes the '
        'barcode adaptor, SC-PCR primers, restriction steps and read coordinates, but not '
        'the Illumina adaptor or library-primer sequences. Those outer final-library regions '
        'are dotted placeholders rather than a guessed TruSeq construct.</div>',
        oligos(), ligation(), paci_and_pcr(), bcivi_and_boundary(), reads(), '</div>',
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
