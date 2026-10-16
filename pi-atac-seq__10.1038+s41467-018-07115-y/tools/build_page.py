#!/usr/bin/env python3
from __future__ import annotations
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import pi_atac as P
import seqprimers as sp
from chemdraw import Scene,annotation_rows,panel,protocol_graph_panels
from page import caveat,head,info
OUT=HERE.parent/"pi-atac-seq.html"
def main():
 lib=P.final_library()
 graph=P.protocol_graph()
 page="\n".join([head("Pi-ATAC-seq library chemistry"),'<div class="wrap" data-protocol-complete="true">','<h1>Pi-ATAC-seq &mdash; joint protein indexing and chromatin accessibility</h1>',info('Defining source: <a href="https://doi.org/10.1038/s41467-018-07115-y">Chen et al., <i>Nature Communications</i> (2018)</a>.'),
 caveat('<b>Inferred adapter structure.</b> The paper reports 96 × 90 barcoding-primer combinations; the exact Pi-ATAC well-to-primer assignment remains unassigned in this model. The <a href="https://doi.org/10.17504/protocols.io.q6gdzbw">authors’ detailed protocol</a> names Ad1/Ad2 stocks and links the <a href="https://doi.org/10.1038/nature14590">Buenrostro scATAC predecessor</a>, whose Supplementary Table 1 prints a Nextera Ad1/Ad2 primer family. Which members Pi-ATAC used for its 96 × 90 wells remains unverified; dotted regions show this source boundary.'),
 caveat('<b>Published descriptions differ.</b> The <a href="https://doi.org/10.1038/s41467-018-07115-y">article Methods, “Index sorting and library preparation for Pi-ATAC”</a> states 20 µl reverse-crosslink buffer, 65°C overnight and a 25 mM barcoding-primer stock. The <a href="https://doi.org/10.17504/protocols.io.q6gdzbw">author protocol, 23 October 2018, steps 26, 32 and 34</a> states 25 µl buffer, 60°C and 25 µM Ad1/Ad2 stocks. Protocol step 25 also prints proteinase K “200 µg/µl” alongside a 1:100 dilution of 20 mg/ml stock, which calculates to 200 µg/ml; the article’s buffer description instead prints 5 ng/ml. These are distinct source values, including an apparent internal unit error. This molecular schematic assigns none of these disputed concentrations, volumes or temperatures; it does not establish which conditions were performed.'),
 '<h2>Reaction workflow</h2>'+protocol_graph_panels(graph,cls="long"),'<h2>Molecular details</h2><h3>(1) Fix, antibody-stain, and tagment cells in bulk</h3>'+panel(P.tagmented_scene().rows(),cls="small",caption="INFERRED — standard ATAC s5/s7 geometry. Bound Tn5 keeps the tagged DNA inside each fixed cell for sorting."),
 '<h3>(2) Index-sort one cell per well and reverse crosslinks</h3>'+info('Protein abundance is recorded by FACS; it is not encoded in the DNA library.'),
 '<h3>(3) Gap-fill and amplify with the well’s barcoding-primer pair</h3>'+panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="INFERRED — Nextera-compatible final library, P5 to P7'. The Pi-ATAC well-to-primer assignments remain unassembled here; a cited predecessor publishes related Ad1/Ad2 primer families."),
 sp.section(lib,P.SEQ_PRIMERS,heading="Sequencing",intro='Paired-end 2 × 75 with two per-well index reads. Primer landing sites below are computed for the inferred Nextera-compatible structure.',read_lengths={"Read 1":75,"Index 1 (i7)":8,"Index 2 (i5)":8,"Read 2":75},show_read_content=True),'</div>'])
 OUT.write_text(page,encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
