#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import starmap as M
import seqprimers as sp
from chemdraw import panel, protocol_graph_panels
from page import head, info
OUT=HERE.parent/"starmap.html"

def render():
    return "\n".join([head("STARmap chemistry"), '<div class="wrap" data-protocol-complete="true">',
      '<h1>STARmap</h1>', '<p class="research-notes"><a href="01_starmap.html">Research notes</a></p>',
      info('Defining source: <a href="https://doi.org/10.1126/science.aat5691">Wang et al., <i>Science</i> (2018)</a>. Its supplementary methods cite sequence Tables S1 and S2, but the available public PDF stops before those tables. Symbolic bases mark original-panel sequences that remain unverified until the tables are recovered. The workflow follows the supplement&#8217;s thin-section/cell-culture procedure (hybridize, ligate, amplify, acrylate, embed, clear); its thick-tissue procedure embeds and digests before ligation and RCA, uses 5&#8242;-amino primer probes and 1-base decoders, and is not drawn separately.'),
      info('Published discrepancy: the supplement&#8217;s probe-design methods split each 40&#8211;46-nt target into halves of 20&#8211;25 nt, whereas its Fig. S1A legend gives 19&#8211;25 nt per primer or padlock target arm. The drawn target arms are symbolic placeholders of schematic width, not a chosen length; the exact arms await the original Table S1.'),
      '<h2>Paired SNAIL recognition</h2>', panel(M.snail_scene().rows(), cls="long"),
      '<h3>Nucleotide-resolved SNAIL fold</h3>', M.SNAIL_TOPOLOGY.folded_panel(
        'The same complex with every modelled nucleotide shown. N denotes panel-specific '
        'bases awaiting the original Tables S1 and S2. The primer probe remains bound to '
        'RNA and supplies the 12-nt splint/RCA primer; both 6-nt padlock ends dock on it, '
        'and the marked bond is the sole T4-ligase closure.'),
      '<h2>Reaction workflow</h2>', protocol_graph_panels(M.protocol_graph(), cls="long"),
      '<h2>SEDAL interrogation</h2>', panel(M.sedal_scene().rows(), cls="long"),
      sp.read_spec_table([sp.read_spec("Sequencing output", None,
        "six microscopy cycles interrogate the 5-base target code in each anchored rolony", "include")],
        "STARmap has no Illumina R1, I1, I2 or R2 reads."), '</div>'])

if __name__=="__main__": OUT.write_text(render(), encoding="utf-8"); print(f"wrote {OUT}")
