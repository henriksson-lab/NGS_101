#!/usr/bin/env python3
from __future__ import annotations
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_5p_v3 as V
from protocol_views import PrimerBindingSection
from chemdraw import Scene,annotation_rows,oligo,panel,protocol_graph_panels
from page import head,info
OUT=HERE.parent/"10x-5prime-v3.html"
def render():
    lib=V.final_library(); cdna=V.amplified_cdna(); graph=V.protocol_graph()
    return "\n".join([head("10x GEM-X 5' Gene Expression v3 chemistry"),'<div class="wrap" data-protocol-complete="true">','<h1>10x Chromium GEM-X Single Cell 5&#39; Gene Expression v3</h1>',info('Primary chemistry source: <a href="https://www.10xgenomics.com/support/universal-five-prime-gene-expression/documentation/steps/library-prep/chromium-gem-x-single-cell-5-v3-gene-expression-user-guide">10x Genomics User Guide CG000733 Rev A</a>.'),
      '<h2>Reaction workflow</h2>'+protocol_graph_panels(graph,cls="long"),
      '<div class="caveat"><b>Intermediate oligos.</b> The gel-bead TSO and complete final library follow CG000733. Exact current-kit RT Primer B and ligation-mix structures are not established in this declaration; source review remains open. Its kit contents list cDNA Primers PN-2000089 (p. 9); the same part number has forward and reverse sequences in the <a href="https://cdn.10xgenomics.com/image/upload/v1660261285/support-documents/CG000331_ChromiumNextGEMSingleCell5-v2_UserGuide_RevE.pdf">earlier CG000331 Rev E guide</a> (p. 73). Assigning those sequences to GEM-X v3 depends on the part number denoting the same formulation across guides; CG000733 does not print their bases. Dotted regions and steps labelled <b>INFERRED</b> reconstruct intermediates that remain unresolved.</div>',
      '<h2>Key oligo</h2><seq>'+oligo("barcoded GEM-X gel-bead TSO",V.bead_tso(),mods="bead; 3' rGrGrG")+'</seq>',
      '<h2>Catalog-linked cDNA primers</h2><p>Same PN-2000089 in the current kit and predecessor guide; this is a catalog-based assignment.</p><seq>'+''.join(oligo(p.name,list(p)) for p in V.cdna_primers())+'</seq>',
      '<h2>Library construction</h2><h3>(1) INFERRED &mdash; prime reverse transcription from the mRNA poly(A) tail</h3>'+panel(V.capture_scene().rows(),cls="small",caption="INFERRED — the primer illustrates the predecessor architecture; its sequence and length are not an assignment to current Poly-dT RT Primer B."),
      '<h3>(2) INFERRED &mdash; template-switch onto the barcoded gel-bead oligo</h3>'+panel(V.switch_scene().rows(),cls="small",caption="INFERRED — the source-listed bead TSO writes the cell barcode and UMI at the transcript 5' end; dotted bases belong to the conditionally reconstructed RT-primer side.")+panel([*Scene.duplex(list(cdna),label="amplified cDNA").rows(),*annotation_rows(cdna)],cls="small",caption="INFERRED — the bead-derived end is source-listed; the far cDNA end is a conditional predecessor-based reconstruction."),
      '<h3>(3) Fragment, end-repair, A-tail, and ligate adapters</h3><p>Fragmentation (32 &deg;C) and end repair/A-tailing (65 &deg;C) run in one Fragmentation Mix on 25% of the purified cDNA; the workflow above also shows the source bead cleanups and both double-sided SPRIselect 0.6X/0.8X size selections with their discarded fractions. The workflow above records unresolved physical geometry for Ligation Mix PN-2001109. The final library sequence does not identify the supplied adapter strands or their ligating ends.</p>',
      '<h3>(4) Add dual sample indexes by PCR</h3>'+panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="Final GEM-X v3 library as printed by CG000733, with 10-base i5 and i7 sample indexes."),
      PrimerBindingSection("library", "Library sequencing", "Sequencing",
          'Paired-end, dual-index sequencing: Read 1 28 cycles, i7 10 cycles, i5 10 cycles, Read 2 90 cycles. Read 1 contains the 16-base cell barcode and 12-base UMI.').render(graph),'</div>'])
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
