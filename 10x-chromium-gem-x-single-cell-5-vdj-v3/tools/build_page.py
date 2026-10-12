#!/usr/bin/env python3
"""Build the sparse 10x GEM-X 5' V(D)J v3 schematic."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
AUDIT_TOOLS = HERE.parents[1] / "10x-chromium-single-cell-5-vdj" / "tools"
sys.path[:0] = [str(HERE), str(AUDIT_TOOLS), str(HERE.parents[1] / "lib")]

import seqprimers as sp
from chemdraw import (MolecularState, Segment, Workflow, duplex_rows, oligo, panel,
                      workflow_panel)
from page import head, info, table
from vdj_common import VDJ_FORWARD
from vdj_v3 import CHEM as C

OUT = HERE.parent / "10x-5prime-vdj-v3.html"


def sequencing(lib) -> str:
    rows = []
    for primer in C.seq_primers:
        hit = sp.locate(lib, primer)
        if hit is None:
            raise ValueError(f"{primer.name} has no site")
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return ('<h2>Final library and read layout</h2>' +
            info('Paired-end, dual-index sequencing: Read 1 28 cycles, i7 10, i5 10, '
                 'Read 2 90. Read 1 reports the 16-base cell barcode and 12-base UMI.') +
            sp.diagram(lib, C.seq_primers) +
            table(("Read", "Primer", "Primer site", "First bases"), rows) +
            table(("Read 1 cycles", "Content"),
                  [("1–16", "cell barcode"), ("17–28", "UMI")], scroll=False) +
            sp.read_content_table(lib, C.seq_primers,
                {"Read 1":28,"Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":90}))


def render() -> str:
    amplicon = C.vdj_amplicon()
    lib = C.final_library()
    cdna = Workflow(MolecularState("Amplification-ready full-length cDNA",
                                  tuple(duplex_rows(C.cdna_template()))))
    cdna.pcr("Amplify full-length cDNA", C.cdna_pcr(),
             duplex_rows(C.cdna_template()), name="amplified full-length cDNA")
    nested = Workflow(MolecularState("Full-length cDNA selected for V(D)J enrichment",
                                     tuple(duplex_rows(C.vdj_template()))))
    nested.pcr("Outer V(D)J PCR (Mix 1)", C.vdj_pcr(inner=False),
               duplex_rows(C.vdj_template(inner=True)), name="outer V(D)J amplicon")
    nested.pcr("Inner V(D)J PCR (Mix 2)", C.vdj_pcr(inner=True),
               duplex_rows(amplicon), name="nested V(D)J amplicon")
    indexed = Workflow(MolecularState("Adapter-ligated barcode-bearing fragment",
                                      tuple(duplex_rows(C.pre_index_library()))))
    indexed.pcr("Dual-index PCR", C.index_pcr(), duplex_rows(lib),
                name="cluster-ready V(D)J library")
    return "\n".join([
        head("10x GEM-X Single Cell 5' V(D)J v3 chemistry"), '<div class="wrap">',
        '<h1>10x Chromium GEM-X Single Cell 5&#39; V(D)J v3</h1>',
        info('Primary source: <a href="https://www.10xgenomics.com/support">10x Genomics '
             'User Guide CG000734 Rev A</a>.'),
        '<h2>Key oligos</h2><seq>',
        oligo("gel-bead barcoded TSO", C.gel_bead_segments()[:-1], three="rGrGrG-3'"),
        oligo("poly(dT) RT primer", C.polydt_segments()),
        oligo("V(D)J universal forward primer", [Segment("", VDJ_FORWARD)]),
        oligo("TT i5 primer", C.index_p5_segments()),
        oligo("TT i7 primer", C.index_p7_segments()), '</seq>',
        '<h2>Library construction</h2>',
        '<h3>(1) Reverse-transcribe poly(A) RNA and template-switch onto the barcoded TSO</h3>',
        panel(C.capture_scene().rows(), cls="small"),
        panel(C.template_switch_scene().rows(), cls="long",
              caption="The GEM-X TSO installs the cell barcode, 12-base UMI and Read 1 handle."),
        '<h3>(2) Amplify full-length cDNA</h3>',
        workflow_panel(cdna, caption="cDNA amplification primer binding and product"),
        '<h3>(3) Enrich V(D)J with outer and inner constant-region pools</h3>',
        workflow_panel(nested, caption="Nested V(D)J PCR primer binding and products"),
        info('The guide publishes every outer and inner reverse-primer sequence and places '
             'the pools in constant regions, but does not identify the exact transcript '
             'coordinate for each oligo. The shared forward primer lands exactly on the '
             'partial Read 1 handle. Individual pool sequences are retained in the research note.'),
        '<h3>(4) Fragment, end-repair, A-tail and ligate the Read 2 adaptor</h3>',
        panel(C.adapter_scene().rows(), cls="small",
              caption="The half adaptor has a 12-base stem and 3'-T overhang."),
        '<h3>(5) Select the barcode-bearing fragment by dual-index PCR</h3>',
        workflow_panel(indexed, caption="Dual-index primer binding and selected product"),
        info('Only fragments retaining the Read 1/barcode end acquire P5; Read 2 starts at '
             'variable fragmentation points across V(D)J.'),
        '<p class="quiet"><a href="01_10x-5p-vdj-v3.html">Research notes</a></p>',
        sequencing(lib), '</div>',
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
