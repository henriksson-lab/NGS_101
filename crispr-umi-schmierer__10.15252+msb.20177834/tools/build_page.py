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
from chemdraw import (Construct, Row, Scene, Segment, annotation_rows,
                      complement_segments, panel, protocol_graph_panels, strand_row)
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
    graph = C.protocol_graph()
    return "\n".join([
        head("CRISPR-UMI (Schmierer)"), '<div class="wrap" data-protocol-complete="true">',
        '<h1>CRISPR-UMI (Schmierer)</h1>',
        '<p class="research-notes"><a href="01_crispr-umi-schmierer.html">Research notes</a></p>',
        info('Schmierer et al. 2017, <a href="https://doi.org/10.15252/msb.20177834">doi:10.15252/msb.20177834</a>. The paper calls the lineage barcode a Random Sequence Label (RSL).'),
        caveat('This page is only the Schmierer protocol. The unrelated Michlits method and an experimental padlock design remain in the research notes.'),
        '<h2>Published guide-list discrepancy</h2>',
        '<p>The original <a href="https://pmc-oa-opendata.s3.amazonaws.com/PMC5658704.1/MSB-13-945-s002.csv">EV2 guide CSV</a> writes 11 NKX-family guide IDs with an underscore before the guide number, for example <code>NKX2-1_10</code>. The original <a href="https://pmc-oa-opendata.s3.amazonaws.com/PMC5658704.1/MSB-13-945-s003.zip">EV3 gene-set ZIP</a> writes the same 11 with a hyphen, for example <code>NKX2-1-10</code>. Their other 23,321 IDs match. This page uses the paper\'s 23,279 distinct-sequence count and draws a generic guide; it assigns no per-guide ID. For sequence lookup, use the EV2 ID and treat the EV3 form as an unresolved alias. Whether the EV3 spelling was intentional is unknown.</p>',
        '<h2>Guide–RSL cassette</h2>',
        '<p>The Schmierer lab later deposited <a href="https://www.addgene.org/196709/">Addgene #196709</a> as an empty <code>lentiGuide-Puro.3xBsmBI</code> backbone linked to this paper and parent #52963. A <a href="https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3003297">later original study</a> calls its cloning backbone <code>pLenti-Puro-AU-flip-3xBsmBI</code> (Addgene #196709). The guide and random label are cloned separately. Addgene lists a verified full map, but its sequence download requires login; this page has not checked the paper-derived reconstruction against that deposit.</p>',
        '<h3>(1) Clone each guide together with a random 6-nt lineage label</h3>',
        panel([strand_row(cas), *annotation_rows(cas)], cls="long",
              caption="The AU-flip guide cassette contains a built-in Illumina Index 1 site followed immediately by the 6-nt RSL. The guide and RSL therefore travel in the same lentiviral construct."),
        '<h2>Reaction workflow</h2>', protocol_graph_panels(graph, cls="long"),
        '<h2>Completed sequencing library</h2>',
        panel([strand_row(lib), *annotation_rows(lib)], cls="long",
              caption="The finished PCR3 product. The vector-templated Index 1 site places the RSL in the i7 read."),
        '<h3>Custom Read 1 primer</h3>',
        '<p>Naming discrepancy inside the defining paper: the Methods sentence on library preparation calls the custom Read 1 oligo <code>CRISPRSeq</code>, while the primer table printed immediately below it heads the same oligo <code>CRIPSRSEQ</code>. Only one sequence is given (<code>' + C.CUSTOM_SEQ_PRIMER + '</code>), so these are one primer under two spellings; this page uses the primer-table spelling. Which spelling the authors intended is unknown.</p>',
        panel(read1_scene().rows(), cls="small",
              caption="The paper’s CRIPSRSEQ primer ends on the U6 +1 G, so Read 1 starts at guide base 1."),
        sp.section(lib, C.SEQ_PRIMERS,
                   intro="Read 1 reports the 20-nt guide, Index 1 reports the 6-nt RSL, and the forward-strand i5 read reports the sample index.",
                   read_lengths={"Read 1": C.READ1_CYCLES,
                                 "Index 1 (i7)": C.I7_CYCLES,
                                 "Index 2 (i5)": C.I5_CYCLES},
                   show_read_content=False),
        sp.read_content_table(lib, C.RUN_PRIMERS, C.READ_LENGTHS),
        '</div>',
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
