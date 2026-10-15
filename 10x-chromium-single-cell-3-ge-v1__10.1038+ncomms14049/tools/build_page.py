#!/usr/bin/env python3
"""Build the 10x Chromium Single Cell 3' Gene Expression v1 page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import chromium_v1 as C
from protocol_views import PrimerBindingSection
from chemdraw import annotation_rows, panel, protocol_graph_panels, strand_row
from page import head, info

OUT = HERE.parent / "10x-chromium-single-cell-3-ge-v1.html"


def construction() -> str:
    bead = C.bead_oligo()
    lib = C.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Release the barcoded oligo-dT primer in each GEM</h3>
{panel([strand_row(bead), *annotation_rows(bead)], cls="small",
       caption="INFERRED — the v1 guide's final-library figure establishes P7–cell-barcode–Read-2–UMI–poly(T) order; the complete synthesized bead oligo is unresolved.")}
<h3>(2) Reverse-transcribe and template-switch in the GEM</h3>
{panel(C.rt_scene().rows(), cls="small",
       caption="Cell barcode and UMI are installed next to the transcript's poly(A) end; the opposite end receives a template-switch handle.")}
<h3>(3) Amplify cDNA, shear, end-repair, A-tail and ligate the Read-1 adaptor</h3>
{panel([strand_row(C.final_library()), *annotation_rows(C.final_library())], cls="small",
       caption="Only the 3'-end fragment retains the bead-derived P7 side and becomes a P5–P7 library after sample-index PCR.")}
<h3>(4) INFERRED &mdash; final sequence-level reconstruction</h3>
{panel([*C.final_scene().rows(), *annotation_rows(lib)], cls="small",
       caption="INFERRED — dotted standard-adapter bases reconstruct the primary-source library architecture; no v1 oligo sequence is presented as published.")}'''


def main() -> None:
    graph = C.protocol_graph()
    page = "\n".join([head("10x Chromium Single Cell 3' GE v1 chemistry"), '<div class="wrap" data-protocol-complete="true">',
        '<h1>10x Chromium Single Cell 3-prime Gene Expression v1</h1>',
        '<p class="research-notes"><a href="01_10x-3p-v1.html">Research notes and source scope</a></p>',
        info('Defining source: <a href="https://doi.org/10.1038/ncomms14049">Zheng et al., <i>Nature Communications</i> (2017)</a>; protocol authority: 10x Genomics CG00026 Rev B.'),
        '<div class="caveat"><b>Sequence boundary.</b> The paper prints anchored (T)30VN and the 10x barcode whitelist supplies individual barcode members. The checked source set does not give a complete synthesized v1 bead oligo. Reconstructed adapter bases are dotted.</div>',
        '<p>Source review scope: two independent readers checked the Zheng article, its bioRxiv v1 preprint, both preprint supplements, all five published supplements, the original 10x CG00026 Rev B guide, and the 10x-origin 737K barcode whitelist. They followed the article and guide links to those materials. The complete synthesized v1 bead oligo was not found in that set; this conclusion does not identify which oligo was used on any particular bead. <a href="01_10x-3p-v1.html">Source checklist and locations</a>.</p>',
        '<h2>Reaction workflow</h2>', protocol_graph_panels(graph, cls="long"),
        construction(), PrimerBindingSection(
            "library", "Library sequencing", "Sequencing",
            'Read 1: 98 cycles of cDNA; Index 1: 14-nt cell barcode; Index 2: 8-nt sample index; Read 2: 10-nt UMI in the displayed main variant. The paper reports 14 earlier libraries with a 5-nt UMI and 5-cycle Read 2; the bioRxiv v1 Supplementary Table ST2 prints both 5- and 10-nt UMI variants. The v1 guide requires standard Illumina primers.').render(graph),
        '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
