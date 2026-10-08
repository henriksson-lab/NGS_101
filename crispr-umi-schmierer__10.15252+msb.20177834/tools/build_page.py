#!/usr/bin/env python3
"""Build the paper-bound CRISPR-UMI (Schmierer) protocol schematic."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import crisprumi as C
import seqprimers as sp
from chemdraw import Construct, Row, Scene, Segment, annotation_rows, complement_segments, panel, strand_row
from page import caveat, head, info

OUT = HERE.parent / "crisprumi.html"


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def cassette() -> Construct:
    return Construct([
        seg("U6 3′ end", C.U6_OVERLAP[-16:], "r3"),
        seg("sgRNA spacer", "N" * C.SPACER_LEN, "cbc", placeholder=True),
        seg("AU-flip scaffold", C.SCAFFOLD_AU_FLIP, "r2"),
        seg("Pol III terminator", C.TERMINATOR, "me"),
        seg("Illumina Index 1 site", C.ILLUMINA_ADAPTER, "t7"),
        seg("6-nt RSL", "N" * C.RSL_LEN, "umi", placeholder=True),
        seg("vector flank", C.DOWNSTREAM, "w1"),
    ], name="integrated guide–RSL cassette")


def read1_scene() -> Scene:
    top = [seg("Read 1 site", C.READ1_SITE[-6:], "t7"),
           seg("U6 3' end", C.U6_OVERLAP[-23:], "r3"),
           seg("sgRNA spacer", "N" * C.SPACER_LEN, "cbc", placeholder=True),
           seg("scaffold", C.SCAFFOLD_AU_FLIP[:12], "r2")]
    sc = Scene()
    sc.strand("library", complement_segments(top), label="library template", rev=True)
    sc.anneal("CRIPSRSEQ", [seg("Read 1 tail", C.CUSTOM_SEQ_PRIMER[:6], "t7"),
                            seg("U6 anneal", C.CUSTOM_SEQ_PRIMER[6:], "r3")],
              to="library", pair=("U6 anneal", "U6 3' end'"))
    sc.arrow("CRIPSRSEQ", "20-cycle Read 1")
    sc.mark("library", "sgRNA spacer'", "cycle 1 is spacer base 1")
    return sc


def render() -> str:
    cas = cassette()
    lib = C.final_library()
    return "\n".join([
        head("CRISPR-UMI (Schmierer)"), '<div class="wrap">',
        '<h1>CRISPR-UMI (Schmierer)</h1>',
        '<p class="research-notes"><a href="01_crispr-umi-schmierer.html">Research notes</a></p>',
        info('Schmierer et al. 2017, <a href="https://doi.org/10.15252/msb.20177834">doi:10.15252/msb.20177834</a>. The paper calls the lineage barcode a Random Sequence Label (RSL).'),
        caveat('This page is only the Schmierer protocol. The unrelated Michlits method and an experimental padlock design remain in the research notes.'),
        '<h2>Build and deliver the guide–RSL library</h2>',
        '<h3>(1) Clone each guide together with a random 6-nt lineage label</h3>',
        panel([strand_row(cas), *annotation_rows(cas)], cls="long",
              caption="The AU-flip guide cassette contains a built-in Illumina Index 1 site followed immediately by the 6-nt RSL. The guide and RSL therefore travel in the same lentiviral construct."),
        '<h3>(2) Transduce the pooled library and apply the CRISPR screen</h3>',
        panel([Row(chunks=[("one virion: [guide]—[RSL]  ->  one integrated cassette  ->  one labelled cell lineage", None, False)])],
              cls="long", caption="Distinct RSLs attached to one guide distinguish independently transduced lineages."),
        '<h2>Recover the integrated cassette</h2>',
        '<h3>(3) PCR1 amplifies the genomic integration in parallel reactions</h3>',
        panel([Row(chunks=[("PCR1-F  --->  [integrated U6—guide—scaffold—Index 1—RSL cassette]  <---  PCR1-R", None, False)])],
              caption="The paper uses 40 parallel genomic-DNA PCR1 reactions for 14 cycles, then pools them."),
        '<h3>(4) PCR2 and indexed PCR3 complete the flow-cell library</h3>',
        panel([Row(chunks=[("PCR2: 19 cycles  ->  PCR3: P5—i5 + P7, 14 cycles  ->  288-bp product", None, False)])],
              caption="Nested amplification enriches the guide/RSL locus. PCR3 adds the P5 end, sample i5 and P7 end."),
        '<h2>Completed sequencing library</h2>',
        panel([strand_row(lib), *annotation_rows(lib)], cls="long",
              caption="The finished PCR3 product. The vector-templated Index 1 site places the RSL in the i7 read."),
        '<h3>Custom Read 1 primer</h3>',
        panel(read1_scene().rows(), cls="small",
              caption="The paper’s CRIPSRSEQ primer ends on the U6 +1 G, so Read 1 starts at guide base 1."),
        sp.section(lib, C.SEQ_PRIMERS,
                   intro="Read 1 reports the 20-nt guide, Index 1 reports the 6-nt RSL, and the forward-strand i5 read reports the sample index."),
        '</div>',
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
