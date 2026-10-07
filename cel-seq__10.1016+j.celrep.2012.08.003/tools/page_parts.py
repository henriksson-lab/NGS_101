"""Sparse, source-specific renderers for CEL-Seq and CEL-Seq2."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import celseq as C
from chemdraw import Construct, Scene, Segment, annotation_rows, oligo, panel, strand_row
from page import head, info, table


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


PAPERS = {
    "CEL-Seq": ("10.1016/j.celrep.2012.08.003", "Cell Reports", "2012"),
    "CEL-Seq2": ("10.1186/s13059-016-0938-8", "Genome Biology", "2016"),
}


def preamble(protocol: str) -> str:
    doi, journal, year = PAPERS[protocol]
    supporting = ""
    if protocol == "CEL-Seq":
        supporting = (' Supporting primer source: <a href="https://doi.org/'
                      '10.1186/s13059-016-0938-8">CEL-Seq2 Table S2</a>.')
    return f"""<div class="wrap">
<h1>{protocol} &mdash; barcoded 3' RNA-seq by T7 linear amplification</h1>
{info(f'Defining source: <a href="https://doi.org/{doi}">Hashimshony et al., '
      f'<i>{journal}</i> ({year})</a>.{supporting}')}
<div class="caveat"><b>Small-RNA library boundary.</b> The accessible primary sources do
not print the Illumina library-primer sequences. Dotted outer arms are unpublished
placeholders, not the secondary-source TruSeq Small RNA reconstruction.</div>"""


def oligos(protocol: str) -> str:
    rows = [oligo(f"{protocol} RT primer", C.rt_primer_segments(protocol))]
    if protocol == "CEL-Seq2":
        rows.append(oligo("randomhexRT primer", [
            seg("library-RT tail", C.RANDOMHEX_TAIL, "r2"),
            seg("random hexamer", "N" * C.RANDOMHEX_NT, placeholder=True),
        ]))
    return f"<h2>Key oligos</h2><seq>{''.join(rows)}</seq>"


def common_front(protocol: str) -> str:
    ds = C.double_stranded_cdna(protocol)
    arna = C.arna(protocol)
    return f"""<h2>Library construction</h2>
<h3>(1) Reverse transcribe from the anchored, cell-barcoded oligo-dT primer</h3>
{panel(C.rt_scene(protocol).rows(), cls="long",
       caption="The primer installs the T7 promoter and cell identity before amplification.")}
<h3>(2) Make the second strand</h3>
{panel(Scene.duplex(list(ds)).rows(), cls="long",
       caption="Second-strand synthesis makes the T7 promoter double-stranded; cells can now "
               "be pooled without losing their identities.")}
<h3>(3) Amplify linearly by T7 in-vitro transcription</h3>
{panel([strand_row(arna), *annotation_rows(arna)], cls="long",
       caption="Each cDNA template produces many antisense-RNA copies.")}"""


def cel1_back() -> str:
    ligated = C.cel1_ligated_arna()
    final = C.final_library("CEL-Seq")
    reads = [
        ("Read 1", "cell barcode, then poly(T); exact primer sequence unavailable"),
        ("Index", "library index; exact primer and index design unavailable"),
        ("Read 2", "transcript sequence from the 3' tag"),
    ]
    return f"""<h3>(4) Fragment aRNA and ligate the RNA 3' adaptor</h3>
{panel([strand_row(ligated), *annotation_rows(ligated)], cls="long",
       caption="INFERRED — the defining chemistry ligates a second adaptor to fragmented aRNA, "
               "but its exact bases are unavailable in the accessible primary sources.")}
<h3>(5) Reverse transcribe and library-PCR</h3>
{panel([strand_row(final), *annotation_rows(final)], cls="long",
       caption="INFERRED — supported inner order with unpublished library-PCR arms shown as "
               "dotted placeholders.")}
<h2>Read layout</h2>
{table(("Read", "Reports"), reads)}"""


def cel2_back() -> str:
    final = C.final_library("CEL-Seq2")
    r1 = [seg("UMI", "N" * C.CEL2_UMI_NT, "umi", placeholder=True),
          seg("cell barcode", "N" * C.CEL2_BARCODE_NT, "cbc", placeholder=True),
          seg("poly(T)", "TTT")]
    r2 = [seg("transcript", "XXXXXXXX...", placeholder=True)]
    reads = [
        ("Read 1", "15 nt: 6-nt UMI + 6-nt cell barcode + three T"),
        ("Index", "7 cycles: indexed library identity"),
        ("Read 2", "36 nt: transcript sequence from the random-primed end"),
    ]
    return f"""<h3>(4) Fragment aRNA and reverse transcribe with randomhexRT</h3>
{panel(C.cel2_random_rt_scene().rows(), cls="long",
       caption="Only a 5'-most aRNA fragment carries the UMI and cell barcode. The tailed "
               "random hexamer adds the opposite library handle without RNA ligation.")}
<h3>(5) Library PCR</h3>
{panel([strand_row(final), *annotation_rows(final)], cls="long",
       caption="INFERRED — supported inner order with unpublished library-PCR arms shown as "
               "dotted placeholders.")}
<h2>Read layout</h2>
{panel([strand_row(Construct(r1), prefix="Read 1  ", suffix=""),
        strand_row(Construct(r2), prefix="Read 2  ", suffix="")],
       cls="small", caption="Paired-end sequencing keeps cell identity in Read 1 and "
                            "transcript sequence in Read 2.")}
{table(("Read", "Published layout"), reads)}"""


def render_page(protocol: str) -> str:
    if protocol not in PAPERS:
        raise ValueError(f"unknown CEL-Seq protocol: {protocol!r}")
    back = cel1_back() if protocol == "CEL-Seq" else cel2_back()
    return "\n".join([head(f"{protocol} library chemistry"), preamble(protocol),
                      oligos(protocol), common_front(protocol), back, "</div>"])
