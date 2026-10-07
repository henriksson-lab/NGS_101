#!/usr/bin/env python3
from __future__ import annotations
import html,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_v4 as V
import seqprimers as sp
from chemdraw import Scene,annotation_rows,oligo,panel
from page import head,info,table
OUT=HERE.parent/"10x-3prime-v4.html"
def sequencing(lib):
    rows=[]
    for p in V.SEQ_PRIMERS:
        h=sp.locate(lib,p)
        if h is None: raise ValueError(f"{p.name} has no site")
        rows.append((html.escape(p.role),html.escape(p.name),html.escape(", ".join(h.covers)),f"<code>{html.escape(h.reads)}</code>&hellip;"))
    return '<h2>Sequencing</h2>'+info('Paired-end, dual-index sequencing: Read 1 28 cycles, i7 10 cycles, i5 10 cycles, Read 2 90 cycles. Read 1 contains the 16-base cell barcode and 12-base UMI.')+table(("Read","Primer","Primer site","First bases"),rows)+table(("Read 1 cycles","Content"),V.read1_layout(),scroll=False)
def render():
    lib=V.final_library(); cdna=V.amplified_cdna()
    return "\n".join([head("10x GEM-X 3' Gene Expression v4 chemistry"),'<div class="wrap">','<h1>10x Chromium GEM-X Single Cell 3&#39; Gene Expression v4</h1>',info('Primary chemistry source: <a href="https://www.10xgenomics.com/support/universal-three-prime-gene-expression/documentation/steps/library-prep/chromium-gem-x-single-cell-3-v4-gene-expression-user-guide">10x Genomics User Guide CG000731</a>.'),
      '<div class="caveat"><b>Intermediate oligos.</b> CG000731 prints the bead oligo and complete final library, but not the TSO or ligation-adapter oligos. Dotted intermediate regions and steps labelled <b>INFERRED</b> reconstruct those roles from the final product.</div>',
      '<h2>Key oligo</h2><seq>'+oligo("GEM-X gel-bead oligo",V.bead_oligo())+'</seq>',
      '<h2>Library construction</h2><h3>(1) Capture and barcode mRNA during GEM reverse transcription</h3>'+panel(V.capture_scene().rows(),cls="small",caption="Each molecule receives a 16-base cell barcode and 12-base UMI."),
      '<h3>(2) INFERRED &mdash; template-switch, break GEMs, and amplify cDNA</h3>'+panel(V.switch_scene().rows(),cls="small",caption="INFERRED — CG000731 does not print the v4 TSO sequence; dotted bases show the v3-compatible reconstruction.")+panel([*Scene.duplex(list(cdna),label="amplified cDNA").rows(),*annotation_rows(cdna)],cls="small",caption="INFERRED — amplified cDNA retains the source-listed bead end and an unpublished TSO-derived end."),
      '<h3>(3) INFERRED &mdash; fragment, end-repair, A-tail, and ligate the Read 2 adapter</h3>'+panel(V.adapter_scene().rows(),cls="small",caption="INFERRED — adapter geometry reconstructed from the complete final product printed in CG000731."),
      '<h3>(4) Add dual sample indexes by PCR</h3>'+panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="Final v4 library as printed by CG000731: a 10-base i5 before Read 1 and a 10-base i7 before P7'."),sequencing(lib),'</div>'])
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
