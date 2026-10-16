#!/usr/bin/env python3
"""Build the diagram-centric CH-ATAC-seq page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import ch_atac as C
import seqprimers as sp
from chemdraw import Scene, annotation_rows, oligo, panel, protocol_graph_panels
from page import head, info

OUT = HERE.parent / "ch-atac-seq.html"


def render() -> str:
    lib = C.final_library()
    graph = C.protocol_graph()
    duplex = Scene.duplex(list(lib), label="library")
    oligos = "\n".join([
        oligo("Indexed Tn5 transfer strand", C.indexed_tn5_transfer()),
        oligo("Tn5 Primer C", C.primer_c()),
        oligo("HY head", C.hy_head()),
        oligo("Barcoded HY oligo", C.barcoded_hy()),
    ])
    return "\n".join([
        head("CH-ATAC-seq library chemistry"), '<div class="wrap" data-protocol-complete="true">',
        '<h1>CH-ATAC-seq</h1>',
        info('Defining source: <a href="https://doi.org/10.1016/j.devcel.2024.01.015">Zhang et al., <i>Developmental Cell</i> (2024)</a>.'),
        info('Every base drawn below is from the defining paper&rsquo;s Supplementary Table 1 (its ordered oligo table) or from <code>lib/nextera.py</code>; the genomic insert is a placeholder. The STAR Methods are behind the publisher paywall, so reaction volumes, cycle numbers and well-handling details are not asserted here.'),
        '<h2>Three indexing rounds</h2>',
        info('Tn5 barcode (384 wells) &rarr; hybridization barcode (768 wells) &rarr; indexed PCR (96 wells). Pools are split again between rounds.'),
        '<h2>Key oligos</h2>',
        info('Supplementary Table 1 writes the indexed Tn5 oligo as <code>A*AGCAGTGG&hellip;</code> with <code>*:Phosphorthioate</code>: the <code>*</code> marks a phosphorothioate <i>bond</i> between the first two bases, so the oligo is 52 bases beginning with the 23-nt SMART handle. The table also marks the PCR P5 primer and the Tn5 ME oligo 5′-phosphorylated. Bonds and 5′ groups are not drawn as bases.'),
        '<seq>', oligos, '</seq>',
        '<h2>Reaction workflow</h2>', protocol_graph_panels(graph, cls="long"),
        '<h2>Molecular details</h2>',
        '<h3>(1) Tagment with indexed-handle and Primer-C transposomes</h3>',
        info('Productive fragments receive the SMART handle and Tn5 barcode at one end, and s7 from Primer C at the other.'),
        '<h3>(2) Hybridize the second barcode</h3>',
        info('A pre-annealed HY duplex uses its SMART-handle overhang to capture the indexed Tn5 end. The HY head supplies s5; a blocking oligo occupies free capture overhangs.'),
        '<h3>(3) Fill, amplify and add MGI ends</h3>',
        panel([*duplex.rows(), *annotation_rows(lib)], cls="small", caption="Final DNBSEQ library assembled from Supplementary Table 1. The HY/Tn5 junction is the unique 23-nt overlap between the reverse complement of Barcoded_HY_oligo and the indexed Tn5 transfer strand, computed by ch_atac.hy_tn5_junction(); only the genomic insert is a placeholder."),
        sp.section(lib, C.SEQ_PRIMERS, heading="Sequencing",
                   intro='All three primers are the standard Nextera ones and are located here by computation, not transcription. The 10-cycle index read follows from the 10-nt i7. <b>The CH-ATAC cycle counts remain unverified</b> because the defining STAR Methods and linked preprint have not been read; paired 100 is assumed, on the grounds that Read 1 must pass 62 constant cycles before reaching genomic DNA and that the sibling CH-RNA-seq paper reports 51 cycles of Read 1 with a dark reaction from 11 to 33 plus 100 cycles of Read 2. The dark windows below (11&ndash;33 over the shared handle, 44&ndash;62 over the Tn5 barcode and mosaic end) are computed from the segment lengths.',
                   required_roles=("Read 1", "Index 1 (i7)", "Read 2"),
                   read_lengths={"Read 1": 100, "Index 1 (i7)": 10, "Read 2": 100},
                   show_read_content=True),
        '</div>'
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
