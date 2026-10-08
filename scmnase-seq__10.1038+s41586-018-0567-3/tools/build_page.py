#!/usr/bin/env python3
"""Build the diagram-centric scMNase-seq chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]

import scmnase_seq as S
import seqprimers as sp
from chemdraw import Scene, annotation_rows, oligo, panel
from page import head, info, table

OUT = HERE.parent / "scmnase-seq.html"


def sequencing(lib) -> str:
    rows = []
    for primer in S.SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        if hit is None:
            raise ValueError(f"{primer.name} has no site")
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return "".join([
        '<h2>Sequencing</h2>',
        info("Paired-end TruSeq reads interrogate the protected fragment from opposite ends; the single i7 read identifies the library."),
        sp.diagram(lib, S.SEQ_PRIMERS),
        table(("Read", "Primer", "Primer site", "First bases"), rows),
    ])


def render() -> str:
    lib = S.final_library()
    return "\n".join([
        head("scMNase-seq library chemistry"), '<div class="wrap">',
        '<h1>scMNase-seq &mdash; MNase-protected DNA from one cell</h1>',
        info('Defining source: <a href="https://doi.org/10.1038/s41586-018-0567-3">Lai et al., <i>Nature</i> (2018)</a>.'),
        '<div class="caveat"><b>Adapter sequence.</b> The defining paper names universal adapters and indexed primers but does not print their bases. Dotted regions and steps labelled <b>INFERRED</b> use the TruSeq reconstruction recorded by upstream scg_lib_structs.</div>',
        '<h2>Inferred key oligos</h2><seq>',
        oligo("Y-adapter, 3'-T strand", S.adapter_bottom()),
        oligo("Y-adapter, phosphorylated strand", S.adapter_top(), mods="/5Phos/"),
        oligo("PCR Primer 1.0", S.p5_primer()),
        oligo("indexed multiplexing primer", S.indexed_p7_primer()),
        '</seq>',
        '<h2>Library construction</h2>',
        '<h3>(1) Digest one lysed cell with micrococcal nuclease</h3>',
        panel(S.digested_fragment_scene().rows(), cls="small",
              caption="MNase removes exposed linker DNA and leaves protected genomic fragments."),
        '<h3>(2) INFERRED &mdash; end-repair and add one 3-prime A at each end</h3>',
        panel(S.a_tailed_fragment_scene().rows(), cls="small",
              caption="INFERRED — the dA-tail intermediate follows from the reconstructed 3'-T adapter."),
        '<h3>(3) INFERRED &mdash; ligate the forked adapter</h3>',
        panel(S.adapter_scene().rows(), cls="small",
              caption="INFERRED — the two adapter oligos form a 12-bp stem and leave a 3'-T ligation overhang."),
        '<h3>(4) Amplify with indexed primers and size-select</h3>',
        panel([*Scene.duplex(list(lib), label="library").rows(), *annotation_rows(lib)],
              cls="small", caption="INFERRED — PCR-completed scMNase-seq library, P5 to P7'. The i7 sample index identifies the cell."),
        sequencing(lib),
        '</div>',
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
