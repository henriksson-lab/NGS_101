"""Sparse individual page renderers for MARS-seq and MARS-seq2.0."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import marsseq as M
import seqprimers as sp
from chemdraw import (Construct, Scene, Segment, annotation_rows, junction_row, oligo, panel,
                      strand_row)
from page import head, info, table


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


PAPERS = {
    "MARS-seq": ("10.1126/science.1247651", "Science", "2014"),
    "MARS-seq2.0": ("10.1038/s41596-019-0164-4", "Nature Protocols", "2019"),
}


def render_page(protocol: str) -> str:
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
    ds = M.ds_cdna(protocol)
    lig = M.ligated_arna(protocol)
    ligation_adaptor = "ordered pool barcode"
    final = M.final_library(protocol)
    if exact:
        r1 = Construct([seg("random diversity", "N" * 5, placeholder=True),
                        seg("pool barcode", M.POOL_READ, "cbc"),
                        seg("cDNA", "XXXXXXXX...", placeholder=True)])
        r2 = Construct([seg("well barcode", "N" * 7, "cbc", placeholder=True),
                        seg("UMI", "N" * 8, "umi", placeholder=True),
                        seg("poly(T), ignored", "TTTTT")])
        read_caption = "Published MARS-seq2.0 read design: 5I.4P.M and 7W.8R.5I."
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
        head(f"{protocol} library chemistry"), '<div class="wrap">',
        f'<h1>{protocol} &mdash; three-level barcoded RNA-seq with T7 amplification</h1>',
        info(f'Defining source: <a href="https://doi.org/{doi}">{author}, '
             f'<i>{journal}</i> ({year})</a>.'),
        f'<div class="caveat"><b>Sequence boundary.</b> {caveat}</div>',
        f'<h2>Key oligos</h2><seq>{"".join(oligos)}</seq>',
        '<h2>Library construction</h2>',
        '<h3>(1) Well-barcoded RT1 primes at the poly(A) tail</h3>',
        panel(M.rt_scene(protocol).rows(), cls="long",
              caption="RT1 installs the T7 promoter, cell identity and molecular tag."),
        '<h3>(2) Second-strand synthesis makes the T7 promoter double-stranded</h3>',
        panel(Scene.duplex(list(ds)).rows(), cls="long"),
        '<h3>(3) T7 IVT linearly amplifies, then the aRNA is fragmented</h3>',
        panel([strand_row(Construct([seg("barcoded 5' aRNA fragment",
                                          "[T7 TAG / WELL BC / UMI]XXXXXXXX...",
                                          placeholder=True)]))], cls="long"),
        '<h3>(4) Ligate a pool-barcoded adapter to the fragmented aRNA</h3>',
        panel([strand_row(lig), junction_row(lig, "antisense insert", ligation_adaptor),
               *annotation_rows(lig)], cls="long",
              caption="The 5'-phosphorylated, 3'-blocked adapter adds pool identity."),
        '<h3>(5) RT2, then library PCR</h3>',
        panel([strand_row(final), *annotation_rows(final)], cls="long",
              caption=("INFERRED — supported inner order with unpublished PCR arms dotted."
                       if exact else "Final MARS-seq library; no index read is used.")),
        (sp.unavailable_diagram(final, "the defining sources do not print the library-primer sequences")
         if exact else sp.diagram(final, M.SEQ_PRIMERS)),
        '<h2>Read layout</h2>',
        panel([strand_row(r1, prefix="Read 1  ", suffix=""),
               strand_row(r2, prefix="Read 2  ", suffix="")], cls="small",
              caption=read_caption),
        table(("Read", "Reports"), read_rows), '</div>',
    ])
