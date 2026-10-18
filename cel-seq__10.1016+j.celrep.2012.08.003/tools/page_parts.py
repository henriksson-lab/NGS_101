"""Sparse, source-specific renderers for CEL-Seq and CEL-Seq2."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import celseq as C
import seqprimers as sp
from chemdraw import Construct, Segment, oligo, panel, protocol_graph_panels, strand_row
from page import head, info, table
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


PAPERS = {
    "CEL-Seq": ("10.1016/j.celrep.2012.08.003", "Cell Reports", "2012"),
}


def preamble(protocol: str, *, complete: bool = False,
             source: tuple[str, str, str] | None = None) -> str:
    doi, journal, year = source or PAPERS[protocol]
    supporting = ""
    if protocol == "CEL-Seq":
        supporting = (' Primer designs and bench steps: the 2012 Supplemental Information '
                      '&ldquo;CEL-Seq Protocol&rdquo; (all 40 RT primers); the same design is '
                      'restated in <a href="https://doi.org/10.1186/s13059-016-0938-8">CEL-Seq2 '
                      'Table S2</a>.')
    complete_attr = ' data-protocol-complete="true"' if complete else ""
    return f"""<div class="wrap"{complete_attr}>
<h1>{protocol} &mdash; barcoded 3' RNA-seq by T7 linear amplification</h1>
{info(f'Defining source: <a href="https://doi.org/{doi}">Hashimshony et al., '
      f'<i>{journal}</i> ({year})</a>.{supporting}')}
<div class="caveat"><b>Small-RNA library boundary.</b> The 2012 CEL-Seq protocol names
the Illumina TruSeq Small RNA kit reagents (RA3 3&prime; adapter, RTP, RP1, indexed RPIX) but
does not print their sequences (&ldquo;available from Illumina&rdquo;), and the kit revision
used in 2012 has not been read. Dotted outer arms are therefore placeholders, not a
TruSeq Small RNA reconstruction.</div>"""


def oligos(protocol: str, design: C.PrimerDesign | None = None) -> str:
    rows = [oligo(f"{protocol} RT primer", C.rt_primer_segments(protocol, design))]
    if protocol == "CEL-Seq2":
        if design is None:
            raise ValueError("CEL-Seq2 requires its source-specific primer design")
        rows.append(oligo("randomhexRT primer", [
            seg("library-RT tail", design.randomhex_tail, "r2"),
            seg("random hexamer", "N" * design.randomhex_nt, placeholder=True),
        ]))
    return f"<h2>Key oligos</h2><seq>{''.join(rows)}</seq>"


def cel1_readout() -> str:
    final = C.final_library("CEL-Seq")
    reads = [
        ("Read 1", "cell barcode, then poly(T); exact primer sequence unavailable"),
        ("Index", "library index; exact primer and index design unavailable"),
        ("Read 2", "transcript sequence from the 3' tag"),
    ]
    return f"""<h2>Read layout</h2>
{sp.unavailable_diagram(final, "exact library-primer sequences are unavailable")}
{table(("Read", "Reports"), reads)}
{sp.read_spec_table([
    sp.read_spec("Read 1", "at least 15", "8-nt cell barcode, then poly(T)", "extract"),
    sp.read_spec("Index 1 (i7)", None, "library index (RPIX), read when needed", "demultiplex"),
    sp.read_spec("Read 2", 50, "transcript sequence from the 3' tag", "include"),
], "Read roles and lengths are source-declared (HiSeq 2000 paired-end: at least 15 bases "
   "for Read 1, 50 for Read 2, plus the Illumina index when needed); exact primer "
   "sequences are unavailable.")}"""


def graph_for(protocol: str, **kw):
    if protocol != "CEL-Seq": raise ValueError("native graph is original CEL-seq only")
    from native import protocol_graph
    return protocol_graph()


def render_page(protocol: str, **kw):
    if protocol != "CEL-Seq": raise ValueError("native page is original CEL-seq only")
    import native
    from batch_page import render
    return render(native)
