#!/usr/bin/env python3
"""Shared renderer for individual inDrop v1 and v2 protocol pages."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import indrop as I
from chemdraw import (Construct, Scene, Segment, annotation_rows, complement_segments,
                      junction_row, oligo, panel, protocol_graph_panels, revcomp,
                      strand_row)
import illumina as il
from page import head, info, table
from protocol_graph import AnnotationSlot
from protocol_views import render_annotations
import seqprimers as sp

PAPERS = {
    I.V1: ("10.1016/j.cell.2015.04.044", "inDrop v1"),
    I.V2: ("10.1038/nprot.2016.154", "inDrop v2"),
}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def preamble(protocol: str) -> str:
    doi, title = PAPERS[protocol]
    notes = ("01_indrop.html" if protocol == I.V1 else
             "../indrop-v1__10.1016+j.cell.2015.04.044/01_indrop.html")
    caveat = ('''<div class="caveat"><b>Sources of the v1 sequences.</b> Every v1 oligo drawn here is printed in the Cell 2015 Supplemental Information, Table S5 (hydrogel primer, barcode templates, RNA ligation oligo, 2nd RT primer and both PCR enrichment primers). Table S5 numbers barcode 1 as ten positions, so a 10-nt barcode 1 is drawn; the 2017 v2 plate sheet instead lists 8&ndash;11 nt. Table S5 prints no sequencing primer: the Read 1 and Read 2 primers below are canonical Illumina primers derived from the library ends, not author-stated. The Supplemental Information 'DNA-Library Preparation' prints the post-RT step order and its modifications (ExoI with HinFI, SPRI clean-ups, no DNase after IVT, Shrimp Alkaline Phosphatase before ligation, 10&ndash;13 PCR cycles); the remaining enzymes and conditions follow the cited CEL-Seq/MARS-seq protocol and are not printed.</div>'''
              if protocol == I.V1 else "")
    # The v2 graph owns its bead-QC discrepancy; the v1 page only notes the shared reference.
    bead_qc = ('''<div class="caveat"><b>Published bead-QC size conflict.</b> Zilionis et al. 2017 Supplementary Method 2 lists the intermediate bead band as 102 nt in the gel-imaging instruction and 101 nt in Anticipated Results; both passages list the full and unextended bands as 135 and 65 nt. The bead primer shown here keeps the published variable 8–11 nt first barcode, so no single intermediate length is asserted. The v1 bead drawing instead follows Cell 2015 Table S5, which uses the same bead-primer and template sequences.</div>'''
               if protocol == I.V1 else "")
    return f'''<div class="wrap"><h1>{title} &mdash; droplet barcoding with photo-released hydrogel primers</h1>
<p class="research-notes"><a href="{notes}">Research notes</a></p>
{info(f'Defining source: <a href="https://doi.org/{doi}">doi:{doi}</a>.')}
{info('Each bead carries a split-pool cell barcode, a 6-nt UMI, oligo-dT and a T7 promoter. UV releases the primers into the droplet; pooled cDNA is amplified linearly by IVT.')}{caveat}{bead_qc}'''


def bead_oligos(protocol: str = I.V2) -> str:
    bc1 = I.V1_BC1_LEN if protocol == I.V1 else 8
    rows = [
        oligo("Acrydite bead primer", [seg("leader", I.LEADER), seg("T7 promoter", I.T7_PROMOTER, "tso"), seg("middle", I.MID), seg("PE1 partial", I.PE1[:24], "r2")], mods="Acrydite; photocleavable spacer"),
        oligo("P1 plate", [seg("W1*", I.W1_STAR, "r3"), seg("barcode 1", "A" * bc1, "cbc", placeholder=True), seg("PE1*", revcomp(I.PE1), "r2")]),
        oligo("P2 plate", [seg("BA19", "B" + "A" * 19, placeholder=True), seg("N6", "N" * 6, "umi", placeholder=True), seg("barcode 2", "B" * 8, "cbc", placeholder=True), seg("W1*", I.W1_STAR, "r3")]),
    ]
    return "<h2>Bead synthesis oligos</h2><seq>" + "".join(rows) + "</seq>"


def library_oligos(protocol: str) -> str:
    if protocol == I.V1:
        rows = [
            oligo("RNA ligation oligo (RLO)", [seg("RLO", I.RLO, "r2")], mods="5Phos; 3SpC3"),
            oligo("2nd RT primer", [seg("RT2", I.V1_RT2, "r2")]),
            oligo("PCR enrichment primer 1", [seg("P5", il.P5, "p5"), seg("PE1", I.V1_PCR1[len(il.P5):], "r1")]),
            oligo("PCR enrichment primer 2", [seg("P7", il.P7, "p7"), seg("RT2", I.V1_PCR2[len(il.P7):], "r2")]),
        ]
    else:
        rows = [
            oligo("PE2-N6", [seg("PE2", I.PE2, "r1"), seg("N6", "N" * 6, "umi", placeholder=True)]),
            oligo("PE2 PCR primer", [seg("P5", il.P5, "p5"), seg("GGTC", "GGTC"), seg("PE2 anneal", I.PE2[:18], "r1")]),
            oligo("PE1 indexed PCR primer", [seg("P7", il.P7, "p7"), seg("i7 as ordered", "I" * 6, "cbc", placeholder=True), seg("PE1 anneal", I.PE1[:16], "r2")]),
            oligo("Custom Read 1 / Index / Read 2", [seg("R1", I.V2_R1, "r1"), seg(" / ", "", placeholder=True), seg("Index", I.V2_INDEX, "r2"), seg(" / ", "", placeholder=True), seg("R2", I.V2_R2, "r2")]),
        ]
    return "<h2>Library oligos</h2><seq>" + "".join(rows) + "</seq>"


def bead_steps(protocol: str) -> str:
    bead = I.bead_primer(I.V1_BC1_LEN) if protocol == I.V1 else I.bead_primer()
    caption = "Finished bead oligo; barcode 1 is 8-11 nt, so one representative length is drawn."
    if protocol == I.V1:
        caption = "Finished v1 bead oligo as assembled in Cell 2015 Table S5 (barcode 1 numbered as ten positions there)."
    return f'''<h2>Bead and droplet steps</h2>
<h3>(1) Build the bead barcode in two 384-way split-pool extensions</h3>{panel([strand_row(bead), *annotation_rows(bead)], caption=caption)}
<h3>(2) Co-encapsulate one cell and one bead; release primers with UV</h3>{info('The photocleavable spacer frees many copies of one bead barcode into the droplet.')}
<h3>(3) Reverse-transcribe from T19V, then make the T7 promoter double-stranded</h3>{panel([strand_row(Construct([*bead.segments, seg("antisense cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True)]))], caption="First strand carries bead barcode and UMI upstream of antisense cDNA.")}
<h3>(4) T7 in-vitro transcription</h3>{info('IVT linearly amplifies the pooled material as antisense RNA carrying the UMI and both bead-barcode blocks.')}
'''


def v1_steps() -> str:
    lig = Construct([seg("aRNA fragment", "XXXXXXXX...XXXXXXXX", placeholder=True), seg("RLO", I.RLO, "r2")])
    return f'''<h2>v1 library conversion</h2>
<h3>(5) Fragment aRNA and ligate the blocked RLO</h3>{panel([strand_row(lig), junction_row(lig, "aRNA fragment", "RLO"), *annotation_rows(lig)], caption="Cell 2015 Table S5 RLO (/5Phos/ … /3SpC3/) ligated to the aRNA 3' end; ligase and conditions follow the cited CEL-Seq/MARS-seq protocol.")}
<h3>(6) Second RT from RLO, then P5/P7 enrichment PCR</h3>{final_panel(I.V1)}'''


def v2_steps() -> str:
    rt = Construct([seg("PE2", I.PE2, "r1"), seg("N6", "N" * 6, "umi", placeholder=True), seg("sense cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True), seg("bead-meta copy", "AAAA...AAAA", placeholder=True)])
    return f'''<h2>v2 library conversion</h2>
<h3>(5) Random-prime the aRNA with PE2-N6</h3>{panel([strand_row(rt), *annotation_rows(rt)], caption="PE2-N6 makes sense-oriented cDNA from the antisense aRNA; the opposite end retains the bead metadata.")}
<h3>(6) Add P5 on the PE2 end and indexed P7 on the PE1 end</h3>{final_panel(I.V2)}'''


def final_panel(protocol: str) -> str:
    lib = I.library(protocol)
    caption = "Final library; one representative 8-nt barcode-1 length is drawn."
    if protocol == I.V1:
        caption = "Final v1 library computed from the Cell 2015 Table S5 oligos; one 10-nt barcode 1 is drawn."
    return panel([strand_row(lib), strand_row(lib, "bottom"), *annotation_rows(lib)], caption=caption)


def sequencing(protocol: str) -> str:
    lib = I.library(protocol); primers = I.V1_PRIMERS if protocol == I.V1 else I.V2_PRIMERS
    if any(sp.locate(lib, p) is None for p in primers): raise ValueError("inDrop read site missing")
    if protocol == I.V1:
        return (sp.section(
            lib, primers, heading="Read layout",
            intro=("Read 1 reports the two bead-barcode blocks and UMI; Read 2 "
                   "reports cDNA (read assignment also stated in the authors' GEO GSE65525 "
                   "data processing; Supplemental Information 'Sequencing and Data Filtering'). "
                   "Table S5 prints no sequencing primer, so the canonical Illumina Read 1 and "
                   "PE Read 2 primers are assumed from the library ends; no per-read cycle count is stated."),
            required_roles=tuple(p.role for p in primers),
            read_lengths=I.V1_READ_LENGTHS, read_uses=I.V1_READ_USES,
            show_read_content=True) + '</div>')
    rows = [("Read 1", "not stated", "sense cDNA"), ("Index 1", "not stated (6-nt listed index)", "sample index"), ("Read 2", "not stated", "barcode 1 (8&ndash;11) &middot; W1 (22) &middot; barcode 2 (8) &middot; UMI (6) &middot; poly(T)")]
    specs = [
        sp.read_spec("Read 1", "not stated", "sense cDNA", "include"),
        sp.read_spec("Index 1 (i7)", "not stated", "6-nt sample index as listed in Supplementary Table", "demultiplex"),
        sp.read_spec("Read 2", "not stated", "barcode 1, W1, barcode 2, UMI, then poly(T)", "extract"),
    ]
    return ('<h2>Read layout</h2>' + sp.diagram(lib, primers)
            + table(("Read", "Cycles", "Content"), rows)
            + sp.read_spec_table(specs, "Acquired cycle counts are not stated in the inspected supplements (MOESM454/455/460); the main protocol text was not available. Read contents follow the oligo table and the published indrops.py read anatomy.")
            + '</div>')


def render_page(protocol: str, out: Path, graph=None) -> None:
    graph = graph if graph is not None else I.protocol_graph(protocol)
    version_steps = v1_steps() if protocol == I.V1 else v2_steps()
    html = "\n".join([head(f"{protocol} library chemistry"), preamble(protocol),
                       render_annotations(graph, AnnotationSlot.INTRO),
                       bead_oligos(protocol), library_oligos(protocol),
                       '<h2>Reaction workflow</h2>', protocol_graph_panels(graph, cls="long"),
                       '<h2>Molecular details</h2>', bead_steps(protocol), version_steps,
                       sequencing(protocol)])
    html = html.replace('<div class="wrap">',
                        '<div class="wrap" data-protocol-complete="true">', 1)
    out.write_text(html, encoding="utf-8")
    print(f"wrote {out}  ({out.stat().st_size:,} bytes)")
