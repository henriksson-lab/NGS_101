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
                      protocol_graph_panels, workflow_panel)
from page import head, info, table
from vdj_common import VDJ_FORWARD
from vdj_v3 import ADAPTER_CAVEAT, CHEM as C, RT_PRIMER_CAVEAT, protocol_graph

OUT = HERE.parent / "10x-5prime-vdj-v3.html"


def sequencing(graph) -> str:
    plan = graph.sequencing["index_pcr_product"]
    lib = plan.library
    rows = []
    for primer in plan.primers:
        hit = sp.locate(lib, primer)
        if hit is None:
            raise ValueError(f"{primer.name} has no site")
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return ('<h2>Final library and read layout</h2>' +
            info('Paired-end, dual-index sequencing: Read 1 28 cycles, i7 10, i5 10, '
                 'Read 2 90. Read 1 reports the 16-base cell barcode and 12-base UMI.') +
            sp.diagram(lib, plan.primers) +
            table(("Read", "Primer", "Primer site", "First bases"), rows) +
            table(("Read 1 cycles", "Content"),
                  [("1–16", "cell barcode"), ("17–28", "UMI")], scroll=False) +
            sp.read_content_table(lib, plan.primers, plan.read_lengths,
                                  acquired_roles=plan.acquired_roles))


def render() -> str:
    amplicon = C.vdj_amplicon()
    lib = C.final_library()
    graph = protocol_graph()
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
        head("10x GEM-X Single Cell 5' V(D)J v3 chemistry"), '<div class="wrap" data-protocol-complete="true">',
        '<h1>10x Chromium GEM-X Single Cell 5&#39; V(D)J v3</h1>',
        info('Primary source: <a href="https://www.10xgenomics.com/support">10x Genomics '
             'User Guide CG000734 Rev A</a>.'),
        info('Inferred region: ' + html.escape(RT_PRIMER_CAVEAT) +
             ' Inferred bases are styled as inferred in the oligo list and step 1.'),
        info('Inferred region: ' + html.escape(ADAPTER_CAVEAT) +
             ' Its bases are styled as inferred in step 4 and the workflow.'),
        '<h2>Reaction workflow</h2>', protocol_graph_panels(graph, cls="long"),
        '<h2>Key oligos</h2><seq>',
        oligo("gel-bead barcoded TSO", C.gel_bead_segments()[:-1], three="rGrGrG-3'"),
        oligo("Poly-dT RT Primer B (inferred from v2 PN-2000007)", C.polydt_segments()),
        oligo("V(D)J universal forward primer", [Segment("", VDJ_FORWARD)]),
        oligo("TT i5 primer", C.index_p5_segments()),
        oligo("TT i7 primer", C.index_p7_segments()), '</seq>',
        '<h2>Library construction</h2>',
        '<h3>(1) Reverse-transcribe poly(A) RNA and template-switch onto the barcoded TSO</h3>',
        panel(C.capture_scene().rows(), cls="small",
              caption="INFERRED — RT primer sequence: Poly-dT RT Primer B is drawn with the "
                      "v2 PN-2000007 sequence because CG000734 does not print its own."),
        panel(C.template_switch_scene().rows(), cls="long",
              caption="The GEM-X TSO installs the cell barcode, 12-base UMI and Read 1 handle."),
        '<h3>(2) Amplify full-length cDNA</h3>',
        workflow_panel(cdna, caption="cDNA amplification primer binding and product"),
        info('CG000734 Rev A names cDNA Primers PN-2000089, the same part as v2; the '
             'drawn primer sequences are printed under that PN in CG000331 Rev E Appendix p.72.'),
        '<h3>(3) Enrich V(D)J with outer and inner constant-region pools</h3>',
        workflow_panel(nested, caption="Nested V(D)J PCR primer binding and products"),
        info('The v3 guide names the v2-labelled enrichment mixes by part number and its '
             'image-only Additional Sequences page (CG000734 Rev A Appendix p.113) prints '
             'the universal forward primer and every outer and inner reverse primer; they '
             'are identical to CG000331 Rev E pp.73&ndash;74. The guides place the pools in '
             'constant regions but do not identify the exact transcript '
             'coordinate for each oligo. The shared forward primer lands exactly on the '
             'partial Read 1 handle. Individual pool sequences are retained in the research note.'),
        '<h3>(4) Fragment, end-repair, A-tail and ligate the Read 2 adaptor</h3>',
        panel(C.adapter_scene().rows(), cls="small",
              caption="INFERRED — Read 2 adaptor: CG000734 premixes it in Ligation Mix "
                      "PN-2001109; the v2 PN-2000094 half adaptor (12-base stem, 3'-T "
                      "overhang) is drawn."),
        '<h3>(5) Select the barcode-bearing fragment by dual-index PCR</h3>',
        workflow_panel(indexed, caption="Dual-index primer binding and selected product"),
        info('Only fragments retaining the Read 1/barcode end acquire P5; Read 2 starts at '
             'variable fragmentation points across V(D)J.'),
        '<p class="quiet"><a href="01_10x-5p-vdj-v3.html">Research notes</a></p>',
        sequencing(graph), '</div>',
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
