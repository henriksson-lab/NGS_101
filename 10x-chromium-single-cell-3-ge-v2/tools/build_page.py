#!/usr/bin/env python3
from __future__ import annotations
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_v2 as V
from protocol_views import PrimerBindingSection
from chemdraw import Scene,annotation_rows,oligo,panel,protocol_graph_panels
from page import head,info
OUT=HERE.parent/"10x-3prime-v2.html"

def render():
    lib=V.final_library(); cdna=V.amplified_cdna(); graph=V.protocol_graph()
    return "\n".join([head("10x Chromium 3' Gene Expression v2 chemistry"),'<div class="wrap" data-protocol-complete="true">',
      '<h1>10x Chromium Single Cell 3&#39; Gene Expression v2</h1>',
      info('Molecular structure: <a href="https://www.10xgenomics.com/support">10x Genomics Technical Note CG000108 Rev A</a>. The later <a href="https://assets.ctfassets.net/an68im79xiti/UhAMGmlaEMmYMaA4A4Uwa/274a813b81e42cba81345d49380432d7/CG00052_SingleCell3_ReagentKitv2UserGuide_RevD.pdf">CG00052 Rev D user guide</a> gives v2 reaction programs (including 67 °C cDNA-PCR annealing) and a 26/8/98 run; differences from the 2017 guide revision remain unverified.'),
      info('The displayed read layout follows CG00052 Rev D: R1 26 / i7 8 / R2 98. '
           '<a href="https://cdn.10xgenomics.com/image/upload/v1660261286/support-documents/CG000120_10x_Technical_Note_SequencingMetrics_NovaSeq_SC3v2_RevA.pdf">CG000120 Rev A, Table 2</a> uses R2 91 with a NovaSeq 100-cycle kit. '
           '<a href="https://cdn.10xgenomics.com/image/upload/v1660261285/support-documents/CG000080_10x_Technical_Note_Base_Composition_SC3_v2_RevB.pdf">CG000080 Rev B, Table 1</a> also describes 150 / 8 / 150 runs; Cell Ranger uses only the first 26 R1 cycles.'),
      '<h2>Reaction workflow</h2>'+protocol_graph_panels(graph,cls="long"),
      '<h2>Key oligos</h2><seq>'+oligo("gel-bead oligo",V.bead_oligo())+oligo("template-switch oligo",V.tso(),mods="3' rGrGrG")+'</seq>',
      '<h2>Library construction</h2><h3>(1) Capture and barcode mRNA during GEM reverse transcription</h3>'+panel(V.capture_scene().rows(),cls="small",caption="Each dissolved gel bead supplies one 16-base cell barcode; each captured molecule receives a 10-base UMI."),
      '<h3>(2) Template-switch, break GEMs, and amplify cDNA</h3>'+panel(V.switch_scene().rows(),cls="small",caption="The TSO supplies the second cDNA-amplification handle.")+panel([*Scene.duplex(list(cdna),label="amplified cDNA").rows(),*annotation_rows(cdna)],cls="small",caption="Amplified full-length cDNA retains barcode and UMI at the transcript's 3' end."),
      '<h3>(3) Fragment, end-repair, A-tail, and ligate the Read 2 adapter</h3>'+panel(V.adapter_scene().rows(),cls="small",caption="The short strand forms a 12-bp stem and leaves a 3'-T overhang for the dA-tailed cDNA fragment."),
      '<h3>(4) Select the barcode-bearing end by sample-index PCR</h3>'+panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="Final v2 library. P5-side priming selects the fragment that retains Partial Read 1, cell barcode and UMI."),
      PrimerBindingSection("library", "Library sequencing", "Sequencing",
          'Read 1 is 26 cycles: 16-base cell barcode then 10-base UMI. Index 1 identifies the sample; Read 2 sequences cDNA.').render(graph),'</div>'])

def main():
    OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
