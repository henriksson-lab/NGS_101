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
from chemdraw import (Construct, Scene, Segment, annotation_rows, oligo, panel,
                      protocol_graph_panels, strand_row)
from endprep import ligation_product_scene
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
    ligated = ligation_product_scene(
        final, "left 11-nt constant", "right 3' A", label="adapter-ligated library",
        adapter_group="same unresolved Illumina adapter")
    return f"""<h3>(4) BciVI removes the SC-PCR ends and leaves 3'-A overhangs</h3>
{panel(scene.rows(), cls="long",
       caption="Type-IIS cutting retains, at each insert end, an 11-nt constant, the 8-nt "
               "cell barcode and the TTAA PacI half-site.")}
<h3>(5) Ligate Illumina adaptors and perform library PCR</h3>
{panel([*ligated.rows(), *annotation_rows(final)], cls="long",
       caption="INFERRED — the matching background at both ends denotes two copies of the "
               "same Illumina adapter ligated to the displayed 3'-A-ended product. The paper "
               "and its supplements leave the selected Illumina adapter and library-primer bases unresolved, so the dotted outer "
               "arms show that source limit.")}"""


def reads() -> str:
    final = D.inferred_final_library()
    rows = [
        ("Read 1", "8-nt cell barcode at bases 1–8, then TTAA and genomic DNA"),
        ("Read 2", "11-nt constant, 8-nt cell barcode at bases 12–19, then TTAA and genomic DNA"),
        ("Index", "sample index distinguishes the pooled library; selected sequence unresolved from checked article and supplements"),
    ]
    r1 = Construct([seg("cell barcode", "N" * 8, "cbc", placeholder=True),
                    seg("PacI half", "TTAA"),
                    seg("genomic DNA", "XXXXXXXX...", placeholder=True)], name="Read 1")
    r2 = Construct([seg("11 constant cycles", "XXXXXXXXXXX", placeholder=True),
                    seg("cell barcode", "N" * 8, "cbc", placeholder=True),
                    seg("PacI half", "TTAA"),
                    seg("genomic DNA", "XXXXXXXX...", placeholder=True)], name="Read 2")
    return f"""<h2>Read layout</h2>
{sp.unavailable_diagram(final, "the selected Illumina library-primer bases remain unresolved from the checked article and supplements")}
{panel([strand_row(r1, prefix="Read 1  ", suffix=""),
        strand_row(r2, prefix="Read 2  ", suffix="")],
       cls="small", caption="Published barcode coordinates and Supplementary Fig. 2B read boundaries: "
                            "collected Read 1 begins at the barcode after the 11-nt constant; "
                            "Read 2 includes that constant before its barcode. The Methods states "
                            "11 initial dark cycles but does not specify their per-read run setup.")}
{table(("Read", "Published layout"), rows)}
{sp.read_spec_table([
    sp.read_spec("Read 1", None, "8-nt cell barcode at cycles 1–8, TTAA, then genomic DNA", "extract_then_include"),
    sp.read_spec("Read 2", None, "11-nt constant, 8-nt cell barcode at cycles 12–19, TTAA, then genomic DNA", "extract_then_include"),
    sp.read_spec("Index 1 (i7)", None, "sample index", "demultiplex"),
], "Barcode coordinates are published; total read and index lengths and exact primer sequences are not.")}"""


def render() -> str:
    graph = D.protocol_graph()
    return "\n".join([
        head("Drop-ChIP library chemistry"), '<div class="wrap" data-protocol-complete="true">',
        '<h1>Drop-ChIP &mdash; barcoded nucleosomal DNA for single-cell ChIP-seq</h1>',
        info('Defining source: <a href="https://doi.org/10.1038/nbt.3383">Rotem et al., '
             '<i>Nature Biotechnology</i> (2015)</a>.'),
        '<div class="caveat"><b>Downstream adaptor boundary.</b> The paper publishes the '
        'barcode adaptor, SC-PCR primers, restriction steps and read coordinates, but not '
        'the Illumina adaptor or library-primer sequences. Those outer final-library regions '
        'are dotted placeholders rather than a guessed TruSeq construct.</div>',
        '<div class="caveat"><b>Published source discrepancies.</b> Fig. 3A calls the '
        'barcode adaptors 64 bp, while Supplementary Table 2 lists 1,152 forward strands '
        'of 60 nt; the model uses the listed sequences and cannot establish a four-base '
        'extension. Online Methods calls the BciVI product a 5′-end A overhang, while '
        'the cut geometry of the printed BciVI motif yields a 3′-A overhang in the model. '
        'The intended published end orientation remains uncertain.</div>',
        '<div class="caveat"><b>Published Supplementary Fig. 2A pairing discrepancy.</b> '
        'In the original publisher PDF, page 2, the top green priming blocks are '
        '<code>GGGGGACC</code> and <code>GGTGGGGG</code>, while both lower rows print '
        '<code>GGGGGTGG</code> and <code>CCAGGGGG</code>. Those lower strings do not pair '
        'with the top; their complements would be <code>CCCCCTGG</code> and '
        '<code>CCACCCCC</code>. Supplementary Table 2 supplies the 1,152 forward '
        'oligos with the top-strand constant, so this model uses a listed top strand '
        'and computes its complement. The intended lower figure motifs remain uncertain. '
        'Fig. 2B labels the outer ends only “Illumina adapter”, leaving those bases unresolved.</div>',
        oligos(), '<h2>Reaction workflow</h2>', protocol_graph_panels(graph, cls="long"),
        reads(), '</div>',
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
