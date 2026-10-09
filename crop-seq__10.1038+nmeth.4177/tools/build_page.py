#!/usr/bin/env python3
"""Build the diagram-centric page for the original CROP-seq protocol."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))
import crop_seq as C
import seqprimers as sp
from chemdraw import Construct, Scene, annotation_rows, oligo, panel, strand_row
from page import head, info, table
OUT = HERE.parent / "crop-seq.html"


def one_strand(con, caption: str) -> str:
    return panel([strand_row(con), *annotation_rows(con)], cls="small", caption=caption)


def vector_trick() -> str:
    viral, provirus = C.transfer_genome().viral_rna(), C.transfer_genome().provirus()
    cassette = Construct(list(C.guide_expression_cassette()), name="guide cassette")
    return f'''<h2>The vector trick</h2>
<h3>(1) Clone each guide between the hU6 promoter and scaffold</h3>
<seq>{oligo("pooled 74-nt Gibson oligo", C.pooled_guide_oligo())}</seq>
{one_strand(cassette, "The 20-nt guide is the only variable region. The 5' homology ends with the added Pol-III +1 G; the 3' homology is the start of the gRNA scaffold.")}
<h3>(2) Package a transfer RNA whose 3′ U3 contains that cassette</h3>
{one_strand(viral, "Packaged transfer RNA, shown as DNA letters. Its incomplete LTR ends are R-U5 and U3-R.")}
<h3>(3) Reverse transcription rebuilds both LTRs from the 3′ U3</h3>
{one_strand(provirus, "Integrated provirus. The model accepts the U3 cassette once and derives both copies; the second guide cassette is not hand-entered.")}
'''


def transcripts() -> str:
    return f'''<h2>Two RNAs from the integrated cassette</h2>
<h3>Pol III: the functional guide</h3>
{one_strand(C.pol3_guide_transcript(), "hU6 produces the non-polyadenylated guide RNA used by Cas9.")}
<h3>Pol II: the guide-identifying transcript</h3>
{one_strand(C.pol2_reporter_transcript(), "The same guide cassette lies inside the polyadenylated puromycin-resistance transcript, so oligo-dT single-cell RNA-seq can capture it.")}
'''


def dropseq_library() -> str:
    frag, lib = C.guide_bearing_cdna_fragment(), C.final_library()
    return f'''<h2>Drop-seq library</h2>
<h3>(4) Capture the polyadenylated reporter on a barcoded bead</h3>
{info("CROP-seq uses Drop-seq: bead oligo-dT reverse transcription, template switching, SMART PCR, Nextera XT tagmentation, then selective P5-SMART + N70x enrichment.")}
{one_strand(frag, "Representative bead-end cDNA whose random Tn5 cut leaves the guide within Read 2 reach. The cDNA is antisense; sequencing reports the guide in its original sense.")}
<h3>(5) Select the bead end and add flow-cell adapters</h3>
{panel([*Scene.duplex(list(lib), label="library").rows(), *annotation_rows(lib)], cls="small", caption="Final guide-bearing CROP-seq library. Cell barcode and UMI are on Read 1; guide identity is in Read 2.")}
'''


def sequencing() -> str:
    return sp.section(C.final_library(), C.sequencing_primers(), heading="Sequencing",
        intro="The Supplementary Protocol specifies 20 cycles of custom Read 1, 8 cycles of Index 1, and 64 cycles of Read 2. Read 2 also supplies transcriptome sequence.",
        required_roles=("Read 1", "Index 1 (i7)", "Read 2"),
        read_lengths={"Read 1": 20, "Index 1 (i7)": 8, "Read 2": 64})


def render() -> str:
    return "\n".join([head("CROP-seq library chemistry"), '<div class="wrap">',
        '<h1>CROP-seq &mdash; pooled CRISPR screening with single-cell RNA readout</h1>',
        '<p class="research-notes"><a href="01_crop-seq.html">Research notes</a></p>',
        info('Defining source: <a href="https://doi.org/10.1038/nmeth.4177">Datlinger et al., <i>Nature Methods</i> (2017)</a>. The original implementation uses Drop-seq.'),
        vector_trick(), transcripts(), dropseq_library(), sequencing(), '</div>'])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
