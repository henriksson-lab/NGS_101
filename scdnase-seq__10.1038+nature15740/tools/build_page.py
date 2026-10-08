#!/usr/bin/env python3
"""Build the diagram-centric scDNase-seq chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]

import scdnase_seq as S
import seqprimers as sp
from chemdraw import Scene, oligo, panel, annotation_rows
from page import head, info, table

OUT = HERE.parent / "scdnase-seq.html"


def sequencing() -> str:
    lib = S.final_library()
    rows = []
    for primer, hit in S.primer_landings():
        if hit is None:
            raise ValueError(f"{primer.name} has no site")
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return "".join([
        '<h2>Sequencing</h2>',
        info("TruSeq Read 1 and Read 2 interrogate the genomic fragment from opposite ends; the single i7 read identifies the library."),
        sp.diagram(lib, S.SEQ_PRIMERS),
        table(("Read", "Primer", "Primer site", "First bases"), rows),
    ])


def render() -> str:
    lib = S.final_library()
    return "\n".join([
        head("scDNase-seq library chemistry"), '<div class="wrap">',
        '<h1>scDNase-seq (Pico-Seq) &mdash; DNase-accessible DNA from one cell</h1>',
        info('Defining source: <a href="https://doi.org/10.1038/nature15740">Jin et al., <i>Nature</i> (2015)</a>.'),
        '<div class="caveat"><b>Adapter sequence.</b> The defining paper names Illumina kits and indexed primers but does not print their bases. Dotted regions and steps labelled <b>INFERRED</b> use the TruSeq reconstruction recorded by upstream scg_lib_structs.</div>',
        '<h2>Inferred key oligos</h2><seq>',
        oligo("Y-adapter, 3'-T strand", S.adapter_bottom()),
        oligo("Y-adapter, phosphorylated strand", S.adapter_top(), mods="/5Phos/"),
        oligo("PCR Primer 1.0", S.p5_primer()),
        oligo("indexed multiplexing primer", S.indexed_p7_primer()),
        '</seq>',
        '<h2>Library construction</h2>',
        '<h3>(1) Digest one lysed cell with DNase I</h3>',
        panel(S.digested_fragment_scene().rows(), cls="small",
              caption="DNase I releases short fragments from accessible genomic DNA; circular carrier DNA is added before purification."),
        '<h3>(2) INFERRED &mdash; end-repair and add one 3-prime A at each end</h3>',
        panel(S.a_tailed_fragment_scene().rows(), cls="small",
              caption="INFERRED — the paper names end repair; the dA tail follows from the reconstructed 3'-T adapter."),
        '<h3>(3) INFERRED &mdash; ligate the forked adapter</h3>',
        panel(S.adapter_scene().rows(), cls="small",
              caption="INFERRED — the two adapter oligos form a 12-bp stem and leave a 3'-T ligation overhang."),
        '<h3>(4) Index-amplify, size-select, then amplify with P5 and P7</h3>',
        panel([*Scene.duplex(list(lib), label="library").rows(), *annotation_rows(lib)],
              cls="small", caption="INFERRED — PCR-completed scDNase-seq library, P5 to P7'. The i7 sample index identifies the cell."),
        sequencing(),
        '</div>',
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
