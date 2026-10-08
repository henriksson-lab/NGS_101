#!/usr/bin/env python3
"""Build the diagram-centric page for transcriptome-wide RNA FISSEQ."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))

import fisseq as F
from chemdraw import annotation_rows, circle_rows, oligo, panel, strand_row
from page import head, info, table

OUT = HERE.parent / "fisseq.html"


def one_strand(con, caption: str) -> str:
    return panel([strand_row(con), *annotation_rows(con)], cls="small", caption=caption)


def construction() -> str:
    cdna = F.linear_cdna()
    return f'''<h2>Make and immobilize adapter-tagged cDNA</h2>
<h3>(1) Random-prime reverse transcription in fixed cells or tissue</h3>
<seq>{oligo("5′-phosphorylated RT primer", list(cdna)[:2], mods="5Phos")}</seq>
{panel(F.rt_scene().rows(), cls="small", caption="The 5′ adapter does not pair to RNA. M-MuLV extends the random hexamer and incorporates aminoallyl-dUTP into cDNA.")}
<h3>(2) Cross-link cDNA, then remove the RNA</h3>
{one_strand(cdna, "BS(PEG)9 anchors the aminoallyl-bearing cDNA in place. RNase plus RNase H leaves a linear single-stranded cDNA with a 5′ phosphate and 3′ hydroxyl.")}
'''


def circle_and_rca() -> str:
    cdna, product = F.linear_cdna(), F.rolony()
    rows = [*circle_rows(cdna, "CircLigase II"), *annotation_rows(cdna, prefix_width=6)]
    return f'''<h2>Circularize and make a rolony</h2>
<h3>(3) CircLigase II closes each cDNA on itself</h3>
{panel(rows, cls="small", caption="The return path is the covalent circle. ** marks the exact last-to-first ligation: the cDNA 3′-OH joins the adapter 5′ phosphate.")}
<h3>(4) Anneal the protected RCA primer</h3>
<seq>{oligo("RCA primer", [F.seg("RCA primer", F.RCA_PRIMER, "r2")], three="-3′; * = each of the final two phosphorothioate bonds")}</seq>
{panel(F.rca_priming_scene().rows(), cls="small", caption="The circle is opened only for drawing. The RCA primer is the exact reverse complement of the adapter.")}
<h3>(5) Phi29 makes tandem complementary copies</h3>
{one_strand(product, "Three repeats are drawn; the physical product contains many. The repeat sequence and repeat boundaries are generated from the closed cDNA circle.")}
'''


def sequencing() -> str:
    rows = [(name, f"5′-p {seq}-3′", len(seq)) for name, seq in F.SEQUENCING_PRIMERS]
    return f'''<h2>Sequencing primers on the finished rolony</h2>
{panel(F.sequencing_scene().rows(), cls="long", caption="All five sliding primers are shown on the middle repeat. The three primer-N sites in this three-repeat drawing are found computationally; every physical repeat carries the same nested sites.")}
{info("SOLiD sequencing-by-ligation uses primer N, then four primers shortened successively at the 5′ end. The offsets make the five-base ligation stride interrogate all positions.")}
{table(("Primer", "Sequence", "Length"), rows, scroll=False)}
'''


def render() -> str:
    return "\n".join([head("FISSEQ library chemistry"), '<div class="wrap">',
        '<h1>FISSEQ &mdash; transcriptome sequencing inside fixed cells</h1>',
        '<p class="research-notes"><a href="01_fisseq.html">Research notes</a></p>',
        info('Defining source: <a href="https://doi.org/10.1126/science.1250212">Lee et al., <i>Science</i> (2014)</a>. Exact oligos and reaction order are from the authors\' open <a href="https://doi.org/10.1038/nprot.2014.191"><i>Nature Protocols</i> procedure</a>.'),
        construction(), circle_and_rca(), sequencing(), '</div>'])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
