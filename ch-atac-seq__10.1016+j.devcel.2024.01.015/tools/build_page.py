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
        info('Every base drawn below is from the defining paper&rsquo;s Supplementary Table 1 (its ordered oligo table) or from <code>lib/nextera.py</code>; the genomic insert is a placeholder. Reaction order, enzymes, the 12-cycle PCR and the read recipe follow the paper&rsquo;s STAR Methods (&ldquo;Library construction of CH-ATAC-seq&rdquo; and &ldquo;MGI library preparation and sequencing&rdquo;).'),
        '<h2>Three indexing rounds</h2>',
        info('Tn5 barcode (384 wells: four 96-well plates) &rarr; hybridization barcode (768 wells: eight 96-well plates) &rarr; indexed PCR (one or more 96-well plates). Pools are split again between rounds.'),
        info('Published discrepancy &mdash; number of PCR indices: the Results (p. 795) say the third round uses &ldquo;384 indexed P7 adapters&rdquo;, but Supplementary Table 1 prints a sheet named &ldquo;96 MGI P7 indexed primers&rdquo; with 96 sequences, and the STAR Methods split nuclei into &ldquo;one or more 96-well plates&rdquo;. The model draws the 96 printed primers. How 384 PCR indices would be made from them is not stated.'),
        info('Published discrepancy &mdash; nuclei per PCR well: the Results (p. 795) collect &ldquo;3,000&ndash;6,000 nuclei for each well&rdquo; after the PNK reaction, while the STAR Methods (p. e4) load &ldquo;3,000-10,000 nuclei&rdquo; per well. Cell loading is not part of the molecular model, so neither value changes the drawn chemistry; which range applies to which experiment is not stated.'),
        info('Reduced-index variant (not drawn): for the fresh/fixed human&ndash;mouse cell-line mixtures the STAR Methods state that only 96 Tn5, 384 hybridization and 32 PCR indices were used. The workflow below shows the full 384 &times; 768 &times; 96-primer design.'),
        '<h2>Key oligos</h2>',
        info('Supplementary Table 1 writes the indexed Tn5 oligo as <code>A*AGCAGTGG&hellip;</code> with <code>*:Phosphorthioate</code>: the <code>*</code> marks a phosphorothioate <i>bond</i> between the first two bases, so the oligo is 52 bases beginning with the 23-nt SMART handle. The table also marks the PCR P5 primer and the Tn5 ME oligo 5′-phosphorylated. Bonds and 5′ groups are not drawn as bases.'),
        '<seq>', oligos, '</seq>',
        '<h2>Reaction workflow</h2>', protocol_graph_panels(graph, cls="long"),
        '<h2>Molecular details</h2>',
        '<h3>(1) Tagment with indexed-handle and Primer-C transposomes</h3>',
        info('Productive fragments receive the SMART handle and Tn5 barcode at one end, and s7 from Primer C at the other.'),
        '<h3>(2) Hybridize the second barcode</h3>',
        info('A pre-annealed HY duplex uses its SMART-handle overhang to capture the indexed Tn5 end. The HY head supplies s5. Supplementary Table 1 lists <code>Block_tail_primer_oligo</code> (the bare 23-nt SMART handle) to &ldquo;block free annealing hybridization primers&rdquo;; the STAR Methods hybridize at 37 &deg;C for 90 min and then block in each well (30 min, 37 &deg;C) before pooling.'),
        '<h3>(3) Fill, amplify and add MGI ends</h3>',
        info('The STAR Methods treat the pooled nuclei with T4 polynucleotide kinase and ATP (37 &deg;C, 30 min), re-split them into the final wells, and add a gap-fill mix of T4 DNA polymerase, dNTP and Ampligase DNA ligase (37 &deg;C, 60 min). Proteinase K/SDS lysis (55 &deg;C, 15 min) and a Tween quench follow, then KAPA HiFi PCR in the same well: 72 &deg;C 5 min, 98 &deg;C 3 min, 12 cycles of 98/65/72 &deg;C. Which gaps and nicks the ligase seals is worked out from the oligo layout; the paper does not draw the junctions.'),
        panel([*duplex.rows(), *annotation_rows(lib)], cls="small", caption="Final DNBSEQ library assembled from Supplementary Table 1. The HY/Tn5 junction is the unique 23-nt overlap between the reverse complement of Barcoded_HY_oligo and the indexed Tn5 transfer strand, computed by ch_atac.hy_tn5_junction(); only the genomic insert is a placeholder."),
        sp.section(lib, C.SEQ_PRIMERS, heading="Sequencing",
                   intro='The PCR library is circularized with a VAHTS Circularization Kit for MGI, made into DNA nanoballs and sequenced on DNBSEQ-T7. The STAR Methods print &ldquo;100 cycles of read1 with dark reaction from 11-33, 44-62, 100 cycles of read2&rdquo;. Those windows fall on the shared handle (11&ndash;33) and on the mosaic end (44&ndash;62) of the modelled Read 1, after the HY barcode (1&ndash;10) and Tn5 barcode (34&ndash;43). The paper calls its primers &ldquo;custom TM (Tn5 modified) sequencing primers&rdquo; but does not print their bases. The standard Nextera primers are therefore placed here by computation, not transcription. The 10-cycle index read follows from the 10-nt i7, not from a printed cycle count.',
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
