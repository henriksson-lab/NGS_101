#!/usr/bin/env python3
"""Build the diagram-centric ASTAR-seq chemistry page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import astar as A
from chemdraw import Scene, Segment, annotation_rows, complement_segments, oligo, panel
import illumina as il
import nextera as nx
from page import head, info
import seqprimers as sp

OUT = HERE.parent / "astar-seq.html"


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def preamble() -> str:
    return f'''<div class="wrap">
<h1>ASTAR-seq &mdash; chromatin accessibility and transcriptome from one cell</h1>
{info('Defining source: <a href="https://doi.org/10.1101/829960">Xing et al., 2019</a>.')}
{info('Tagment genomic DNA first, then reverse-transcribe RNA in the same Fluidigm C1 chamber. Biotin separates the two products before their final library reactions.')}
'''


def oligos() -> str:
    rows = [
        oligo("Biotin oligo-dT (C1-P2-T31)",
              [seg("handle", A.C1_HANDLE, "tso"), seg("CG", A.C1_DT_SPACER),
               seg("T31", "T" * A.C1_DT_LEN)], mods=A.C1_DT_MOD),
        oligo("RNA template-switching oligo (C1-P2-RNA-TSO)",
              [seg("handle", A.C1_HANDLE, "tso"), seg("rGrGrG", "rGrGrG", placeholder=True)],
              mods="all-ribo"),
        oligo("Biotin C1 PCR primer (C1-P2-PCR-2)",
              [seg("handle", A.C1_PCR, "tso")], mods=A.C1_PCR_MOD),
        oligo("ATAC adaptor 1 / qPCR primer", [seg("s5", nx.S5, "s5"), seg("ME", nx.ME, "me")]),
        oligo("ATAC adaptor 2 / qPCR primer", [seg("s7", nx.S7, "s7"), seg("ME", nx.ME, "me")]),
        oligo("v2_Ad1.N i5 indexing primer",
              [seg("P5", il.P5, "p5"), seg("i5", "N" * A.INDEX_LEN, "cbc", placeholder=True),
               seg("s5", nx.S5, "s5"), seg("ME prefix", nx.ME[:10], "me")]),
        oligo("v2_Ad2.N i7 indexing primer",
              [seg("P7", il.P7, "p7"), seg("i7 as ordered", "N" * A.INDEX_LEN, "cbc", placeholder=True),
               seg("s7", nx.S7, "s7"), seg("ME prefix", nx.ME[:7], "me")]),
    ]
    return "<h2>Key oligos</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>"


BODY = "XXXXXXXX...XXXXXXXX"


def tagmentation_scene() -> Scene:
    """Representative s5/s7 product with the two 9-nt gaps still open."""
    g = nx.TAGMENTATION_GAP
    top = [seg("s5", nx.S5, "s5"), seg("ME-L", nx.ME, "me"),
           seg("top gap-facing", "X" * g, placeholder=True),
           seg("genomic core", BODY, placeholder=True)]
    bottom = [seg("s7", nx.S7, "s7"), seg("ME-R", nx.ME, "me"),
              seg("bottom gap-facing", "x" * g, placeholder=True),
              *complement_segments([top[-1]])]
    sc = Scene()
    sc.strand("top", top, label="DNA")
    sc.anneal("bottom", bottom, to="top", pair=("genomic core'", "genomic core"), label="DNA")
    sc.anneal("left ME'", complement_segments([top[1]], suffix=""), to="top",
              pair=("ME-L", "ME-L"), label="")
    sc.anneal("right ME'", complement_segments([bottom[1]], suffix=""), to="bottom",
              pair=("ME-R", "ME-R"), label="", above=True)
    sc.mark("top", "top gap-facing", "9-nt gap opposite this end")
    sc.mark("bottom", "bottom gap-facing", "9-nt gap opposite this end")
    return sc


def rt_scenes() -> tuple[Scene, Scene]:
    mrna = [seg("body", BODY, placeholder=True), seg("poly(A)", "A" * A.C1_DT_LEN)]
    primer = [seg("handle", A.C1_HANDLE, "tso"), seg("CG", A.C1_DT_SPACER),
              seg("dT", "T" * A.C1_DT_LEN)]
    first = [*primer, *complement_segments([mrna[0]]), seg("CCC", "CCC", "tso")]

    p = Scene()
    p.strand("mRNA", mrna)
    p.anneal("biotin oligo-dT", primer, to="mRNA", pair=("dT", "poly(A)"),
             mod5="biotin")
    p.arrow("biotin oligo-dT", "SuperScript IV")

    sw = Scene()
    sw.strand("mRNA", mrna)
    sw.anneal("first-strand cDNA", first, to="mRNA", pair=("dT", "poly(A)"),
              mod5="biotin")
    sw.anneal("all-RNA TSO", [seg("TSO handle", A.C1_HANDLE, "tso"),
                               seg("GGG", "GGG", "tso")],
              to="first-strand cDNA", pair=("GGG", "CCC"), above=True)
    sw.arrow("first-strand cDNA", "template switch")
    return p, sw


def cdna_pcr_scene() -> Scene:
    con = A.c1_cdna()
    sc = Scene()
    sc.strand("top", list(con), label="cDNA", mod5="biotin")
    sc.anneal("bottom", complement_segments(list(con)), to="top",
              pair=("C1 handle, oligo-dT end'", "C1 handle, oligo-dT end"),
              label="cDNA", mod5="biotin")
    sc.labels("top")
    return sc


def separated_scene() -> Scene:
    con = A.c1_cdna()
    sc = Scene()
    sc.strand("released strand", complement_segments(list(con)), label="RNA-side cDNA")
    sc.note("released strand", "heat releases the non-biotinylated strand; biotinylated strands remain bead-bound")
    return sc


def final_panel(lib, caption: str) -> str:
    sc = Scene.duplex(lib.segments)
    sc.strands["top"].label = sc.strands["bottom"].label = ""
    return panel([*sc.rows(), *annotation_rows(lib)], caption=caption)


def steps() -> str:
    prime, switch = rt_scenes()
    atac = A.atac_library()
    rna = A.rna_library()
    return f'''<h2>Library generation</h2>
<h3>(1) Lyse and tagment chromatin in each C1 chamber</h3>
{panel(tagmentation_scene().rows(), cls="small", caption="Tn5 inserts s5/ME or s7/ME ends at accessible genomic DNA. The representative amplifiable product is shown before gap fill.")}

<h3>(2) Stop Tn5, anneal the biotin oligo-dT, then reverse-transcribe</h3>
{panel(prime.rows(), cls="small", caption="EDTA stops Tn5. MgCl2 is then supplied with the RT mix; the unanchored T31 primer can anneal within the poly(A) tract.")}
{panel(switch.rows(), cls="small", caption="MMLV adds CCC at the RNA 5' end; the all-RNA TSO's GGG pairs there and contributes the same C1 handle as the oligo-dT primer.")}

<h3>(3) Five cycles of on-chip PCR add biotin to both cDNA strands</h3>
{panel(cdna_pcr_scene().rows(), caption="The single C1 handle primer amplifies between the two terminal handle copies. ATAC fragments lack this handle and are not amplified.")}

<h3>(4) Streptavidin separates the modalities</h3>
{panel(separated_scene().rows(), cls="small", caption="Supernatant: tagmented ATAC-DNA. Beads: cDNA; heating releases its non-biotinylated strand for further amplification.")}

<h3>(5a) ATAC-DNA: gap-fill and index PCR</h3>
{final_panel(atac, "One s5/s7 ATAC product after PCR with a per-cell v2_Ad1.N / v2_Ad2.N index pair.")}

<h3>(5b) cDNA: finish amplification, retagment, and Nextera XT index PCR</h3>
{final_panel(rna, "The separated full-length cDNA is tagmented anew. Only an s5/s7 fragment becomes a dual-indexed RNA library molecule.")}
'''


def read_scene(lib, primer, arrow: str) -> str:
    hit = sp.locate(lib, primer)
    if hit is None:
        raise ValueError(f"{primer.name} does not land")
    sc = Scene.duplex(list(lib), label="")
    if primer.role == "Read 1":
        pseg = [seg("s5 primer", nx.S5, "s5"), seg("ME primer", nx.ME, "me")]
        sc.anneal("primer", pseg, to="bottom", pair=("s5 primer", "s5'"))
    elif primer.role == "Read 2":
        pseg = [seg("s7 primer", nx.S7, "s7"), seg("ME primer", nx.ME, "me")]
        sc.anneal("primer", pseg, to="top", pair=("s7 primer", "s7'"), above=True)
    else:
        raise ValueError(f"no read drawing for {primer.role}")
    sc.arrow("primer", arrow)
    return panel(sc.rows(), cls="small")


def sequencing() -> str:
    atac, rna = A.atac_library(), A.rna_library()
    return f'''<h2>Sequencing</h2>
{info('ATAC: paired-end 2 &times; 50. RNA: paired-end 2 &times; 101. Both use standard Nextera sequencing primers; neither molecule carries a cell barcode or UMI. Cell identity is the C1 chamber\'s dual-index pair.')}
<h3>ATAC Read 1 and Read 2</h3>
{read_scene(atac, sp.NEXTERA["R1"], "Read 1 into accessible genomic DNA")}
{read_scene(atac, sp.NEXTERA["R2"], "Read 2 into the opposite end")}
<h3>RNA Read 1 and Read 2</h3>
{read_scene(rna, sp.NEXTERA["R1"], "Read 1 into cDNA")}
{read_scene(rna, sp.NEXTERA["R2"], "Read 2 into the opposite end")}
{sp.section(atac, A.SEQ_PRIMERS, heading="Sequencing-primer sites", intro="The ATAC and RNA libraries share the same adapter geometry; only the insert differs.")}
</div>'''


def main() -> None:
    OUT.write_text("\n".join([head("ASTAR-seq library chemistry"), preamble(), oligos(),
                              steps(), sequencing()]), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
