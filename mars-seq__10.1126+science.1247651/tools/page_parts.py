"""Sparse individual page renderers for MARS-seq and MARS-seq2.0."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import marsseq as M
import seqprimers as sp
from chemdraw import Construct, Segment, oligo, panel, protocol_graph_panels, strand_row
from page import head, info, table


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


from mars_graph import PAPERS, graph_for as protocol_graph


def render_page(protocol: str) -> str:
    graph = protocol_graph(protocol)
    doi, journal, year = PAPERS[protocol]
    author = "Jaitin et al." if protocol == "MARS-seq" else "Keren-Shaul et al."
    exact = protocol == "MARS-seq2.0"
    caveat = ("The 2014 author supplement publishes RT1, the ligation adapters and both "
              "library-PCR primers; the final library below is assembled from those oligos."
              if not exact else
              "The supplements publish RT1 and the pool-barcoded ligation adapters, but not "
              "RT2 or the PCR-primer sequences; dotted outer arms stop at that boundary.")
    oligos = [oligo("RT1", M.rt1_segments(protocol))]
    if exact:
        oligos.append(oligo("representative pool ligation adapter",
                            list(M.mars2_ligation_adapter()), mods="/5Phos/ … /3SpC3/"))
    else:
        oligos.extend([
            oligo("representative pool ligation adapter",
                  [seg("pool barcode", M.POOL_ORDERED, "cbc"),
                   seg("random diversity", "N" * M.MARS1_POOL_RANDOM_NT,
                       placeholder=True),
                   seg("ligation/RT2 handle", M.LIG_CONSTANT, "r1")],
                  mods="/5Phos/ … /3SpC3/"),
            oligo("P5_Rd1 PCR", [seg("P5 + Read 1", M.MARS1_P5_PCR, "r1")]),
            oligo("P7_Rd2 PCR", [seg("P7 + Read 2", M.MARS1_P7_PCR, "r2")]),
        ])
    final = M.final_library(protocol)
    if exact:
        r1 = Construct([seg("random diversity", "N" * 5, placeholder=True),
                        seg("pool barcode", M.POOL_READ, "cbc"),
                        seg("cDNA", "XXXXXXXX...", placeholder=True)])
        r2 = Construct([seg("well barcode", "N" * 7, "cbc", placeholder=True),
                        seg("UMI", "N" * 8, "umi", placeholder=True),
                        seg("poly(T), ignored", "TTTTT")])
        # The Supplementary Information's read-design fields print the worked example
        # R1_design = 3I.4P.50M with R2_design = 7W.8R.5I.  Read 1's 3I matches the
        # MARS-seq1 N3 adapter, so 5I is what the 2.0 N5 adapter implies -- derived here,
        # never published, and labelled as such.
        read_caption = ("Published MARS-seq2.0 read design: 3I.4P.50M and 7W.8R.5I. "
                        "DERIVED — with this version's N5 ligation adapter Read 1 is "
                        "5I.4P instead; the 5I is derived from the adapter, not published.")
        read_rows = [("Read 1", "N5, 4-nt pool barcode, sense cDNA"),
                     ("Read 2", "7-nt well barcode, 8-nt UMI, poly(T)")]
    else:
        r1 = Construct([seg("random diversity", "N" * 3, placeholder=True),
                        seg("pool barcode", M.POOL_READ, "cbc"),
                        seg("cDNA", "XXXXXXXX...", placeholder=True)])
        r2 = Construct([seg("well barcode", "N" * 6, "cbc", placeholder=True),
                        seg("UMI", "N" * 4, "umi", placeholder=True),
                        seg("poly(T), ignored", "TTTTT")])
        read_caption = "Published MARS-seq read design: N3, 4-nt pool barcode, then cDNA; Read 2 reports cell barcode and UMI."
        read_rows = [("Read 1", "N3, 4-nt pool barcode, sense cDNA"),
                     ("Read 2", "6-nt well barcode, 4-nt UMI, poly(T)")]
    return "\n".join([
        head(f"{protocol} library chemistry"),
        '<div class="wrap" data-protocol-complete="true">',
        f'<h1>{protocol} &mdash; three-level barcoded RNA-seq with T7 amplification</h1>',
        info(f'Defining source: <a href="https://doi.org/{doi}">{author}, '
             f'<i>{journal}</i> ({year})</a>.'),
        f'<div class="caveat"><b>Sequence boundary.</b> {caveat}</div>',
        f'<h2>Key oligos</h2><seq>{"".join(oligos)}</seq>',
        '<h2>Library construction</h2>', protocol_graph_panels(graph),
        (sp.unavailable_diagram(final, "the defining sources do not print the library-primer sequences")
         if exact else sp.diagram(final, M.SEQ_PRIMERS)),
        '<h2>Read layout</h2>',
        panel([strand_row(r1, prefix="Read 1  ", suffix=""),
               strand_row(r2, prefix="Read 2  ", suffix="")], cls="small",
              caption=read_caption),
        table(("Read", "Reports"), read_rows),
        sp.read_spec_table([
            sp.read_spec("Read 1", None, read_rows[0][1], "extract_then_include"),
            sp.read_spec("Read 2", None, read_rows[1][1], "extract"),
        ], "Read contents are source-declared; the defining sources do not fix a complete run recipe."),
        '</div>',
    ])
