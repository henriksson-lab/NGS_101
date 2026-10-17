#!/usr/bin/env python3
"""Generate smallseq.html -- the Small-seq chemistry page.

Every ladder on the page is derived from the constants in smallseq.py. Nothing is
hand-aligned: every multi-strand drawing is a chemdraw.Scene (strands placed by which
segment pairs with which, every column checked), single-strand ladders come from
Construct, and the sizes in the table come from library_bp().
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import smallseq as sm
from chemdraw import (Row, Scene, Segment, annotation_rows, complement_segments,
                      oligo, panel, strand_row, tm)
import seqprimers as sp
from illumina import P5, P7
from page import caveat, head, info, legend, table
from protocol_views import assemble_page, WorkflowSection

OUT = HERE.parent / "smallseq.html"
PRE = 4          # width of the "5'- " strand prefix used by the Construct ladders


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


# ---------------------------------------------------------------- Scene helpers
# Every multi-strand drawing below is a chemdraw.Scene: strands are given 5'->3' and placed
# by naming the segments that pair, so no column is ever typed by hand. Two small local
# additions, both positioned from Scene spans rather than typed:

def rtp_segs() -> list[Segment]:
    """RTP 5'->3' as ordered: biotin-CCTTGGCACCCGAGAATTCC-rA (the last base is RNA)."""
    return [seg("RTP", sm.RTP[:-1], "r2"), seg("rA", sm.RTP[-1])]


def rp1_segs() -> list[Segment]:
    return [seg("P5", P5, "p5"), seg("RP1 foot", sm.RP1[len(P5):], "r1")]


def ligated_segs() -> list[Segment]:
    return list(sm.ligated_rna().segments)


INSERT = sm.insert_placeholder(sm.TYPICAL_MIRNA_NT)
UMI = sm.RA5_UMI_BASE * sm.RA5_UMI_LEN


# =============================================================================== preamble
def preamble() -> str:
    return f"""<div class="wrap">
<h1>Small-seq &mdash; single-cell small-RNA sequencing</h1>

{info("""<b>Small-seq</b> &mdash; Hagemann-Jensen M, Abdullayev I, Sandberg R, Faridani OR,
&ldquo;Small-seq for single-cell small-RNA sequencing&rdquo;, <i>Nat Protoc</i>
2018;<b>13</b>(10):2407&ndash;2424,
<a href="https://doi.org/10.1038/s41596-018-0049-y">doi:10.1038/s41596-018-0049-y</a>.
Original method: Faridani <i>et al.</i>, <i>Nat Biotechnol</i> 2016;<b>34</b>:1264.
Pipeline: <a href="https://github.com/eyay/smallseq">github.com/eyay/smallseq</a>.""")}

{info("""A miRNA is 22 nucleotides of single-stranded RNA with no poly(A) tail, no cap and
no handle. It cannot be tagmented, cannot be primed with oligo-dT, and offers a reverse
transcriptase nothing to hold on to. So Small-seq <b>ligates</b> a defined end onto it
&mdash; twice, in sequence, with no polymerase involved in either step.""")}

{legend("""<b>Four things here appear nowhere else in this repo.</b>
<b>(1)</b> TruSeq <i>Small RNA</i> adapters (<r1>RA5</r1>, <r2>RA3</r2>) rather than the
TruSeq DNA adapter, so the read-1 landing site is a different sequence entirely.
<b>(2)</b> <b>Sequential ligation</b>: 3' adapter, then a digestion step that destroys the
leftovers, then the 5' adapter.
<b>(3)</b> A <umi>UMI</umi> inside the ligated 5' adapter &mdash; the first thing
sequenced, before any biology.
<b>(4)</b> A <w1>masking oligonucleotide</w1> that depletes 5.8S rRNA by occupying its 3'
end, so no adapter can be ligated there. Written up once in
<code>ref/concepts/small-rna-ligation.md</code>.""")}

{info("""<b>Implementation shown: the 2018 detailed protocol.</b> Its Supplementary
Table 1 supplies the 192 eight-base SRX index primers; SR001 is the representative here.
The 2016 defining paper instead used 48 six-base INDEX-RPI primers (Supplementary Table 2)
and 43- or 51-base single reads, and pooled and purified the libraries with a Zymo DNA Clean
&amp; Concentrator-5 column without a described size selection (Online Methods, library
preparation). Its adapter, RT primer, RP1 and mask sequences match the 2018 Reagent setup;
it writes the mask's 3&prime; label as TEG-biotin. The 2018 protocol recommends at least 51 bases
(Steps 37–38). These are published versions, not contradictory source values.""")}
{caveat("""<b>Read 1 is represented by its derived entry site.</b> The displayed RA5
binding span puts the first read base on the UMI, as specified by the detailed protocol.
It does not identify the oligo loaded into an instrument. Instrument-specific primer
instructions remain outside this model's established evidence; this is not a claim that
they were unpublished.""")}
"""


# ================================================================================= oligos
def oligos() -> str:
    rows = [
        oligo("RA5 &mdash; 5' adapter (RNA)",
              [seg("", sm.RA5_HANDLE, "r1"), seg("", UMI, "umi", placeholder=True),
               seg("", sm.RA5_LINKER, "tso")],
              mods="NH2-", three="-3'OH"),
        oligo("RA3 &mdash; 3' adapter (DNA)", [seg("", sm.RA3, "r2")],
              mods="rApp-", three="-ddC"),
        oligo("RTP &mdash; RT primer, and reverse primer of PCR 1",
              [seg("", sm.RTP[:-1], "r2"), seg("", "rA", None, placeholder=True)],
              mods="biotin-"),
        oligo("RP1 &mdash; forward primer of both PCRs",
              [seg("", P5, "p5"), seg("", sm.RP1[len(P5):], "r1")]),
        oligo("SRX &mdash; index primer, one per well",
              [seg("", P7, "p7"), seg("", sm.SRX_001_PRIMER_INDEX, "cbc"),
               seg("", sm.SRX_COMMON_LINK, "r3"), seg("", sm.RTP, "r2")]),
        oligo("5.8S rRNA masking oligo (DNA)", [seg("", sm.MASK_58S, "w1")],
              three="-biotin"),
    ]
    return ("<h2>Oligos</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>\n"
            + info(f"""All six verbatim from <i>Reagent setup</i>, which states
            &ldquo;All oligonucleotides are listed in the 5&prime; to 3&prime;
            direction&rdquo;. <cbc>SRX</cbc> is from the detailed protocol's Supplementary
            Table 1; SR001 is shown and the complete 192-member family is published there.
            <code>r</code> marks a ribonucleotide; <code>rH</code> is rA, rU or rC;
            <code>rApp</code> is a chemically pre-adenylated 5' end.
            Lengths: <r1>RA5</r1> {len(sm.RA5_HANDLE)}&nbsp;+&nbsp;{sm.RA5_UMI_LEN}
            &nbsp;+&nbsp;{len(sm.RA5_LINKER)}&nbsp;=&nbsp;{len(sm.ra5_oligo())}&nbsp;nt,
            <r2>RA3</r2> {len(sm.RA3)}&nbsp;nt, RTP {len(sm.RTP)}&nbsp;nt,
            RP1 {len(sm.RP1)}&nbsp;nt, SRX {len(sm.SRX_SEQUENCE)}&nbsp;nt, mask
            {len(sm.MASK_58S)}&nbsp;nt."""))


def scene_rtp_ra3() -> Scene:
    """RTP on RA3: a blunt 21-bp duplex. Scene raises if any column fails to pair."""
    sc = Scene()
    sc.strand("RA3", [seg("RA3", sm.RA3, "r2")], mod5="rApp", mod3="ddC")
    sc.anneal("RTP", rtp_segs(), to="RA3", pair=("rA", "RA3"), mod5="biotin")
    sc.mark("RTP", "RTP", f"RTP == revcomp(RA3), {len(sm.RA3)} for {len(sm.RA3)}",
            through="rA")
    return sc


def ra5_split() -> list[Segment]:
    """RA5 as ordered, with the handle split where RP1's 3' end falls."""
    n = sm.RP1_ANNEAL_NT
    return [seg("RP1 site", sm.RA5_HANDLE[:n], "r1"), seg("CGATC", sm.RA5_HANDLE[n:], "r1"),
            seg("UMI", UMI, "umi", placeholder=True), seg("CA", sm.RA5_LINKER, "tso")]


def scene_rp1() -> Scene:
    """RP1 does not pair with RA5 -- it has RA5's sequence. It primes on the cDNA copy of
    RA5, so that is what it is annealed to here."""
    sc = Scene()
    ra5 = ra5_split()
    sc.strand("RA5 (RNA)", ra5, mod5="NH2")
    sc.anneal("cDNA", complement_segments(ra5), to="RA5 (RNA)",
              pair=("RP1 site'", "RP1 site"))
    sc.anneal("RP1", rp1_segs(), to="cDNA", pair=("RP1 foot", "RP1 site'"))
    sc.footer(f"^ RP1's 3' end stops here, {len(sm.RA5_HANDLE) - sm.RP1_ANNEAL_NT} nt short of "
         f"RA5's: {sm.RA5_HANDLE[sm.RP1_ANNEAL_NT:]}",
              "cDNA", "CGATC'")
    return sc


def interlock() -> str:
    foot = sm.RP1[len(P5):]
    tail_nt = sm.RA5_HANDLE[len(foot):]
    p = panel([*scene_rtp_ra3().rows(), Row(), *scene_rp1().rows()],
              cls="long", caption="How the four published oligos interlock. RTP pairs with "
              "RA3 directly. RP1 has RA5's own sequence, so it pairs with RA5's cDNA copy "
              "(drawn under the RA5 it was copied from). Every column is placed by pairing "
              "and checked; every relationship is also asserted in the self-test.")
    rows = [
        ["<code>RTP == revcomp(RA3)</code>", "&#9989;",
         "exactly, 21 for 21 &mdash; this is the whole digestion step: RTP anneals to "
         "leftover RA3 and makes lambda exonuclease a blunt duplex to chew on"],
        ["<code>RP1 == P5 + RA5[:21]</code>", "&#9989;",
         f"P5 taken from <code>lib/illumina.py</code> (a vendor document), not from this "
         f"paper. RP1 is the small-RNA P5 arm."],
        [f"RP1 stops {len(tail_nt)} nt short of RA5's 3' end", "&#9989;",
         f"the missing bases are <code>{tail_nt}</code> &mdash; which is why RP1 cannot be "
         f"the read-1 primer"],
        ["<code>RA3</code> begins <code>TGG</code>", "&#9989;",
         "hence the pipeline's precursor filter <code>remove_reads_with_genomic_TGG.py</code>"],
    ]
    thermo = table(
        ["Duplex", "Tm (nearest-neighbour)", "Published temperature it explains"],
        [["RTP &middot; RA3, 21 bp", f"<b>{tm(sm.RTP):.1f}&nbsp;&deg;C</b>",
          f"PCR 1 anneals at <b>{sm.PCR1_ANNEAL_C}&nbsp;&deg;C</b> &mdash; within "
          f"{abs(sm.PCR1_ANNEAL_C - tm(sm.RTP)):.1f}&nbsp;&deg;C"],
         ["RP1, full 50 nt", f"{tm(sm.RP1):.1f}&nbsp;&deg;C",
          f"PCR 2 anneals at {sm.PCR2_ANNEAL_C}&nbsp;&deg;C, below it"],
         [f"RP1's {len(foot)}-nt annealing foot", f"{tm(foot):.1f}&nbsp;&deg;C",
          "the first cycle of PCR 1 is the hard one; after that RP1 binds full-length"],
         [f"masking oligo, {len(sm.MASK_58S)} nt", f"{tm(sm.MASK_58S):.1f}&nbsp;&deg;C",
          f"lysis is {sm.LYSIS_MIN}&nbsp;min at {sm.LYSIS_TEMP_C}&nbsp;&deg;C &mdash; "
          f"the mask stays annealed throughout"]])
    return ("<h2>The oligos interlock exactly</h2>\n" + p
            + table(["Identity", "", "Why it matters"], rows)
            + "<h3>And the annealing temperatures fall out of the sequences</h3>\n" + thermo
            + info("""The cleanest internal consistency in the paper. In PCR 1 the reverse
            primer is the <i>leftover RT primer</i> &mdash; a bare 21-mer &mdash; and the
            annealing temperature is its Tm to within a fifth of a degree. In PCR 2 both
            primers carry full-length tails, and the annealing temperature goes up
            7&nbsp;&deg;C."""))


# ================================================================================== steps
def scene_mask() -> Scene:
    """The mask on the 3' end of 5.8S rRNA. The rRNA is drawn as revcomp(mask) -- what the
    oligo must be annealing to -- plus a placeholder for the rest of the rRNA upstream."""
    up = sm.RRNA_58S_NT - len(sm.MASK_58S)
    rrna = [seg(f"5.8S rRNA, upstream ~{up} nt", "XXXXXX...XXX", None, placeholder=True,
                inferred=True),
            seg("5.8S rRNA, 3'-terminal 76 nt", sm.MASK_58S_TARGET, "w1", inferred=True)]
    sc = Scene()
    sc.strand("5.8S rRNA", rrna, mod3="OH")
    sc.anneal("mask", [seg("mask", sm.MASK_58S, "w1")], to="5.8S rRNA",
              pair=("mask", "5.8S rRNA, 3'-terminal 76 nt"), mod3="biotin")
    sc.note("mask", "3'-biotin blocks extension; the mask's 5' end is flush with the "
                    "rRNA's 3'-OH, which is therefore inside the duplex")
    sc.labels("5.8S rRNA")
    return sc


def step_mask() -> str:
    rows = scene_mask().rows()
    return ("<h2>Step by step</h2>\n"
            + f"<h3>(1) Lysis, {sm.LYSIS_MIN} min at {sm.LYSIS_TEMP_C} &deg;C &mdash; "
              f"and the mask anneals</h3>\n"
            + panel(rows, cls="long",
                    caption="The bottom strand IS the masking oligo, drawn 3'->5' so it pairs "
                            "column for column; its 3'-biotin therefore sits at the left. "
                            "The rRNA's 76 3'-terminal nt are revcomp(mask) -- what the "
                            "oligo must be annealing to. INFERRED: the paper never prints "
                            "the rRNA.")
            + info(f"""2&nbsp;&micro;M masking oligo goes into the 3&nbsp;&micro;l lysis
            buffer with 0.13% Triton X-100 and RNase inhibitor, <i>before</i> the cell is
            sorted in. &ldquo;During this incubation, the RNA molecules are unfolded, and
            the rRNA-masking oligonucleotide binds to the 3&prime; end of the highly
            abundant 5.8S rRNA, preventing ligation of the 3&prime; adapter to these
            (interfering) rRNA molecules.&rdquo; An rRNA whose 3' end is inside a
            {len(sm.MASK_58S)}-bp duplex never gets <r2>RA3</r2>; with no <r2>RA3</r2>
            there is no RT primer site; so it never enters the library. Nothing is pulled
            down and nothing is destroyed &mdash;
            {tm(sm.MASK_58S):.1f}&nbsp;&deg;C of duplex is the entire mechanism.""")
            + caveat("""<b>Unverified.</b> The reverse complement above is computed from the
            published oligo, but checking it base-by-base against 5.8S rRNA needs a
            reference sequence &mdash; which is data, and data is never committed to this
            repo. It is drawn here so the next person can BLAST it rather than retype it.
            The mask is species-specific: the published one is for human and mouse."""))


def scene_lig3() -> Scene:
    """Three separate single strands -- nothing pairs, so each starts at column 0."""
    sc = Scene()
    sc.strand("small RNA", [seg("small RNA", INSERT, None, placeholder=True)],
              mod5="p", mod3="OH")
    sc.strand("+ RA3", [seg("RA3", sm.RA3, "r2")], mod5="rApp", mod3="ddC")
    sc.strand("= ligated", [seg("small RNA", INSERT, None, placeholder=True),
                            seg("RA3", sm.RA3, "r2")], mod5="p", mod3="ddC")
    sc.junction("= ligated", "small RNA", "RA3")
    return sc


def step_lig3() -> str:
    p = panel(scene_lig3().rows(), cls="long",
        caption="(3) T4 RNA ligase 2 truncated KQ, 8% PEG8000, 30 C for 6 h then 4 C "
                "overnight. No ATP in the reaction. The adapter's rApp is the donor; the "
                "small RNA's 3'-OH is the acceptor.")
    return ("<h3>(2) Ligation of the 3' adapter &mdash; 17 h</h3>\n" + p
            + info("""<b>Two consequences of no ATP.</b> The ligase cannot adenylate
            anything, so the only possible donor in the tube is the pre-adenylated adapter
            &mdash; small RNAs cannot be ligated to each other. And the adapter's own 3' end
            is blocked (<code>ddC</code>), so it cannot extend onward.
            &ldquo;It is necessary to order the oligo with the 5&prime; adenylation
            modification. Enzymatic adenylation of the oligo after synthesis may reduce the
            efficiency of the method.&rdquo; A bonus: ligation works through 2&prime;-O-methyl
            3' ends, so piRNAs and plant miRNAs are captured too."""))


def scene_digest() -> Scene:
    """RTP anneals to EVERY RA3 -- free or ligated. Only the free one has a DNA 5'-phosphate
    at a duplex end, so only it is digested."""
    sc = Scene()
    sc.strand("free RA3", [seg("RA3", sm.RA3, "r2")], mod5="p", mod3="ddC")
    sc.anneal("RTP", rtp_segs(), to="free RA3", pair=("rA", "RA3"), mod5="biotin")
    sc.mark("RTP", "RTP", "lambda exonuclease digests the 5'-phosphorylated strand "
                          "-- the RA3 above, not the biotin-blocked RTP", through="rA")
    lig = [seg("small RNA", INSERT, None, placeholder=True), seg("RA3", sm.RA3, "r2")]
    sc.strand("ligated", lig, mod5="p", mod3="ddC")
    sc.anneal("RTP on ligated", rtp_segs(), to="ligated", pair=("rA", "RA3"),
              mod5="biotin", label="RTP")
    sc.mark("RTP on ligated", "RTP", "protected: RA3's 5'-phosphate is internal now, and "
                                     "RTP stays on as the RT primer", through="rA")
    return sc


def step_digest() -> str:
    p = panel(scene_digest().rows(), cls="long",
        caption="(4) 5' deadenylase strips rApp to leave 5'-P; RTP anneals; lambda "
                "exonuclease digests the 5'-phosphorylated strand of the duplex. 30 C "
                "15 min, then 37 C 15 min. RTP's 3'-terminal base is an rA.")
    return ("<h3>(3) Digestion of the unligated 3' adapter &mdash; the step that makes "
            "the protocol work</h3>\n" + p
            + info("""Three enzymes, one target. <b>5&prime; deadenylase</b> converts
            <code>rApp</code> to a plain 5'-phosphate. <b>RTP</b> &mdash; which is the exact
            reverse complement of <r2>RA3</r2> &mdash; anneals to it, giving
            &ldquo;a double-stranded DNA, a preferred substrate for lambda
            exonuclease&rdquo;. <b>Lambda exonuclease</b> then eats the phosphorylated
            strand. RTP survives because &ldquo;The RT primer is blocked at the 5&prime; end,
            which protects it against the exonuclease.&rdquo; Adapter already ligated to a
            small RNA is protected too: its 5'-phosphate is an internal phosphodiester now.
            So the reaction removes exactly the failure mode and nothing else.""")
            + caveat("""<b>Skip this and you have nothing.</b> &ldquo;If this procedure is
            not implemented, the adapter dimer fraction will dominate the sequence
            reads.&rdquo; RTP is left in the tube afterwards and reused twice: to prime RT,
            and as the reverse primer of PCR&nbsp;1."""))


def scene_lig5() -> Scene:
    sc = Scene()
    ra5 = [seg("RA5", sm.RA5_HANDLE, "r1"), seg("UMI", UMI, "umi", placeholder=True),
           seg("CA", sm.RA5_LINKER, "tso")]
    acc = [seg("small RNA", INSERT, None, placeholder=True), seg("RA3", sm.RA3, "r2")]
    sc.strand("RA5", ra5, mod5="NH2", mod3="OH")
    sc.strand("+ acceptor", acc, mod5="p", mod3="ddC")
    sc.strand("= ligated", ra5 + acc, mod5="NH2", mod3="ddC")
    sc.junction("= ligated", "CA", "small RNA")
    sc.labels("= ligated")
    return sc


def step_lig5() -> str:
    p = panel(scene_lig5().rows(), cls="long",
        caption="(5) T4 RNA ligase 1 with ATP, 37 C for 1 h. RA5's 5' end is blocked, its "
                "3'-OH is free, and the acceptor needs a 5'-phosphate.")
    return ("<h3>(4) Ligation of the 5' adapter &mdash; and with it the UMI</h3>\n" + p
            + info(f"""<b>The 5' block is what stops 5'-adapter concatemers</b>, and the
            5'-phosphate requirement is a selectivity feature rather than an accident:
            &ldquo;The 5&prime; adapter does not bind to mRNAs, as these species are capped
            at their 5&prime; end and, therefore, cannot be ligated.&rdquo; mRNA is excluded
            by chemistry, not by size.""")
            + info(f"""<b>The <umi>UMI</umi> is {sm.RA5_UMI_LEN} &times;
            <code>rH</code></b> &mdash; rA, rU or rC, never rG:
            &ldquo;The H mix was chosen over the &lsquo;N&rsquo; mix &hellip; because it
            results in less unwanted primer annealing.&rdquo; That gives
            3<sup>{sm.RA5_UMI_LEN}</sup>&nbsp;=&nbsp;<b>{sm.RA5_UMI_DIVERSITY:,}</b> UMIs,
            the paper's own figure, against 4<sup>8</sup>&nbsp;=&nbsp;{4**8:,} for an 8-mer
            of N. The <tso>CA</tso> that follows &ldquo;serves as a linker between the small
            RNA insert and the UMI&rdquo;."""))


def scene_rt() -> Scene:
    sc = Scene()
    sc.strand("RNA", ligated_segs(), mod5="NH2", mod3="ddC")
    sc.anneal("RTP", rtp_segs(), to="RNA", pair=("rA", "RA3"), mod5="biotin")
    sc.arrow("RTP", "reverse transcription")
    sc.labels("RNA")
    return sc


def scene_pcr1() -> tuple[Scene, Scene]:
    """Cycle 1: RP1 primes on the cDNA's RA5' end. Later cycles: RTP primes on that copy."""
    a = Scene()
    a.strand("cDNA", complement_segments(ligated_segs()), rev=True, mod5="biotin")
    a.anneal("RP1", rp1_segs(), to="cDNA", pair=("RP1 foot", "RA5'"))
    a.arrow("RP1", "cycle 1: copies the cDNA, adding P5")
    b = Scene()
    top = [seg("P5", P5, "p5"), *ligated_segs()]
    b.strand("RP1 copy", top)
    b.anneal("RTP", rtp_segs(), to="RP1 copy", pair=("rA", "RA3"), mod5="biotin")
    b.arrow("RTP", "leftover RT primer = reverse primer")
    b.labels("RP1 copy")
    return a, b


def srx_segs() -> list[Segment]:
    """SRX 5'->3' as the library forces it: P7 + i7 + link (all unpublished) + an
    RTP-like 3' end that anneals over RA3."""
    lib = list(sm.library().segments)
    i = [s.name for s in lib].index("RA3")
    return complement_segments(lib[i:])


def scene_pcr2() -> tuple[Scene, Scene]:
    """Each strand of the PCR-1 product with the primer that copies it."""
    top = [seg("P5", P5, "p5"), *ligated_segs()]
    a = Scene()
    a.strand("top", top)
    a.anneal("SRX", srx_segs(), to="top", pair=("RA3'", "RA3"))
    a.arrow("SRX", "SRX: adds i7 + P7")
    b = Scene()
    b.strand("bottom", complement_segments(top), rev=True)
    b.anneal("RP1", rp1_segs(), to="bottom", pair=("P5", "P5'"))
    b.arrow("RP1", "RP1, now full-length")
    b.labels("bottom")
    return a, b


def step_rt_pcr() -> str:
    p_rt = panel(scene_rt().rows(), cls="long",
        caption="(6) Reverse transcription. RTP was already added at the digestion step; "
                "SuperScript II, 42 C for 1 h, in Taq buffer rather than RT buffer because "
                "carried-over MgCl2 is already in excess.")

    a, b = scene_pcr1()
    p_pcr1 = panel([*a.rows(), Row(), *b.rows()], cls="long",
        caption="(7) PCR 1: only RP1 is added -- the leftover RT primer is the reverse "
                "primer. 13 cycles, annealing at 60 C, which is RTP's Tm. RP1 primes on the "
                "cDNA's 3' end (RA5'), adds P5, and the copy then carries all 26 nt of "
                "RA5 because the cDNA template did.")

    a, b = scene_pcr2()
    p_pcr2 = panel([*a.rows(), Row(), *b.rows()], cls="long",
        caption="(8) PCR 2: the same forward primer plus one of 192 published SRX index "
                "primers per well. 13 cycles, annealing at 67 C. SR001 is shown: its "
                "RTP-like 3' foot anneals over RA3 and its P7-index-common tail overhangs.")
    return ("<h3>(5) Reverse transcription</h3>\n" + p_rt
            + "<h3>(6&ndash;7) Two rounds of PCR</h3>\n" + p_pcr1 + p_pcr2
            + info("""&ldquo;In this reaction, only the forward primer is added to the PCR
            mix, and the leftover RT primer serves as the reverse primer.&rdquo; The second
            round uses &ldquo;the same forward primer &hellip; together with index reverse
            primers&rdquo;, one per well, so all 192 libraries can be pooled."""))


# ================================================================================ library
def final_library() -> str:
    lib = sm.library()
    rows = [strand_row(lib, "top"), strand_row(lib, "bottom")] + annotation_rows(
        lib, prefix_width=PRE)
    pub = f"""<pre>
<align class="long">
<p5>{P5}</p5> <r1>{sm.RA5_HANDLE}</r1> <umi>{UMI}</umi> <tso>{sm.RA5_LINKER}</tso> &hellip; <r2>{sm.RA3}</r2>
</align>
</pre>"""
    sizes = table(
        ["Insert", "Library", "Published size it has to match"],
        [["0 &mdash; adapter dimer", f"<b>{sm.library_bp(0)}&nbsp;bp</b>", "&mdash;"],
         [f"{sm.SMALLRNA_MIN_NT}&nbsp;nt &mdash; shortest analysed",
          f"{sm.library_bp(sm.SMALLRNA_MIN_NT)}&nbsp;bp", "&mdash;"],
         [f"{sm.TYPICAL_MIRNA_NT}&nbsp;nt &mdash; typical miRNA",
          f"<b>{sm.library_bp(sm.TYPICAL_MIRNA_NT)}&nbsp;bp</b>",
          f"cut target {sm.SIZE_SELECT_TARGET_BP[0]}&ndash;"
          f"{sm.SIZE_SELECT_TARGET_BP[1]}&nbsp;bp &#9989;"],
         [f"{sm.SMALLRNA_MAX_NT}&nbsp;nt &mdash; longest &lsquo;small RNA&rsquo;",
          f"{sm.library_bp(sm.SMALLRNA_MAX_NT)}&nbsp;bp",
          f"profile {sm.LIBRARY_PROFILE_BP[0]}&ndash;{sm.LIBRARY_PROFILE_BP[1]}"
          f"&nbsp;bp &#9989;"],
         [f"{sm.RRNA_58S_NT}&nbsp;nt &mdash; 5.8S rRNA, mask omitted",
          f"<b>{sm.library_bp(sm.RRNA_58S_NT)}&nbsp;bp</b>",
          f"peak {sm.RRNA_PEAK_BP[0]}&ndash;{sm.RRNA_PEAK_BP[1]}&nbsp;bp &#9989;"]])
    bracket = table(
        ["Published observation", "Size", "Implied right arm"],
        [[f"miRNA cut window, insert &asymp; {sm.TYPICAL_MIRNA_NT}&nbsp;nt",
          f"{sm.SIZE_SELECT_TARGET_BP[0]}&ndash;{sm.SIZE_SELECT_TARGET_BP[1]}&nbsp;bp",
          f"{sm.SIZE_SELECT_TARGET_BP[0] - sm.LEFT_ARM_NT - sm.TYPICAL_MIRNA_NT}&ndash;"
          f"{sm.SIZE_SELECT_TARGET_BP[1] - sm.LEFT_ARM_NT - sm.TYPICAL_MIRNA_NT}&nbsp;nt"],
         [f"5.8S peak with the mask left out, insert &asymp; {sm.RRNA_58S_NT}&nbsp;nt",
          f"{sm.RRNA_PEAK_BP[0]}&ndash;{sm.RRNA_PEAK_BP[1]}&nbsp;bp",
          f"{sm.RRNA_PEAK_BP[0] - sm.LEFT_ARM_NT - sm.RRNA_58S_NT}&ndash;"
          f"{sm.RRNA_PEAK_BP[1] - sm.LEFT_ARM_NT - sm.RRNA_58S_NT}&nbsp;nt"]])
    return ("<h2>The final library</h2>\n" + panel(rows, cls="long")
            + info(f"""{len(lib)}&nbsp;bp for a {sm.TYPICAL_MIRNA_NT}-nt miRNA. The
            published left arm is {sm.LEFT_ARM_NT}&nbsp;bp and every base of it is
            printed in the paper:""")
            + pub
            + "<h3>Sizes cross-check the published SRX arm</h3>\n"
            + info(f"""<cbc>SRX</cbc>'s 3' end anneals where RTP anneals &mdash; over
            <r2>RA3</r2> &mdash; so what it adds is P7 + index + its 12-nt common linker.
            Two independent size observations provide a separate cross-check:""")
            + bracket
            + info(f"""The windows overlap, which they did not have to: <b>37&ndash;47
            nt</b>. The published arm &mdash; P7&nbsp;(24) + index&nbsp;(8) + common
            link&nbsp;({sm.SRX_LINKER_NT}) &mdash; is
            <b>{sm.SRX_RIGHT_ARM_NT}&nbsp;nt</b>, inside both. The sizes below now derive
            from the printed SR001 sequence rather than a fitted placeholder.""")
            + sizes
            + caveat(f"""<b>This explains the paper's own complaint.</b> It lists as a
            limitation &ldquo;the abundance of adapter dimers in the sequencing reads,
            despite the effects of the digestion step that was designed to reduce the amount
            of these dimers&rdquo;. The dimer is
            <b>{sm.ADAPTER_DIMER_BP}&nbsp;bp</b> &mdash; and the Pippin window
            <i>starts</i> at {sm.PIPPIN_WINDOW_BP[0]}&nbsp;bp, while the manual gel cut is
            {sm.GEL_CUT_BP[0]}&ndash;{sm.GEL_CUT_BP[1]}&nbsp;bp, which contains it outright.
            The shortest real library is only
            {sm.library_bp(sm.SMALLRNA_MIN_NT) - sm.ADAPTER_DIMER_BP}&nbsp;bp away. On a 3%
            agarose cassette that is not a separation, so size selection can
            <i>reduce</i> dimers and never remove them &mdash; exactly what the paper
            reports."""))


# ============================================================================= sequencing
def scene_read1() -> Scene:
    lib = list(sm.library().segments)
    sc = Scene()
    sc.strand("top", lib, label="")
    sc.anneal("bottom", complement_segments(lib), to="top", pair=("P5'", "P5"), label="")
    sc.anneal("read 1 primer", [seg("RA5", sm.READ1_SITE, "r1")], to="bottom",
              pair=("RA5", "RA5'"), label="")
    sc.arrow("read 1 primer", "read 1")
    sc.footer("|" + "-" * (sm.READ_LEN - 2) + f"| {sm.READ_LEN} bp single read", "bottom", "UMI'")
    sc.labels("top")
    return sc


def sequencing() -> str:
    lib = sm.library()
    rows = scene_read1().rows()
    return ("<h2>Sequencing</h2>\n"
            + panel(rows, cls="long", caption="Derived Read 1 entry site: the first called base is the UMI. The drawn RA5 span is a geometry query, not an identified loaded primer.")
            + info(f"""The 2018 detailed protocol's adapter-ligation section places the
            UMI before CA and the RNA insert. Steps 37–38 refer instrument operation to
            the manufacturer's protocol and recommend a single read of at least
            {sm.READ_LEN} bases. Supplementary Table 1 prints the complete SRX index
            primers and the reported eight-base indexes; the Index 1 length used below is
            that barcode length, not a printed acquired index-cycle count. The binding calculations below
            establish compatibility with the declared library; they do not establish which
            sequencing-primer mixture the authors loaded.""")
            + sp.section(lib, sm.SEQ_PRIMERS, intro=f"""Single-read, single-indexed
            2018 implementation. The Read 1 entry is the derived RA5 site ending immediately
            before the UMI. RP1 and the TruSeq DNA entries are binding comparisons.
            Read 1 reports {sm.RA5_UMI_LEN} UMI bases and {len(sm.RA5_LINKER)} CA bases
            before the insert.""",
            required_roles=("Read 1", "Index 1 (i7)"),
            content_primers=(sm.SEQ_PRIMERS[0], sm.SEQ_PRIMERS[3]),
            read_lengths={"Read 1": sm.READ_LEN,
                          "Index 1 (i7)": sm.SRX_INDEX_NT})
            + "<h3>Read processing</h3>\n"
            + table(["Step", "What happens"],
                    [[f"{sm.READ_LEN}-bp single read",
                      f"&ldquo;a single read of {sm.READ_LEN} bp or longer&rdquo;"],
                     [f"trim {sm.TRIM_5P} nt",
                      f"{sm.RA5_UMI_LEN}-nt UMI moved to the read name "
                      f"(<code>--bc-pattern=NNNNNNNN</code>), then "
                      f"<code>cutadapt -u {len(sm.RA5_LINKER)}</code> removes the CA"],
                     ["trim the 3' adapter",
                      f"<code>cutadapt -a file:cutadapt_3prime.fa --minimum-length "
                      f"{sm.SMALLRNA_MIN_NT}</code>"],
                     ["map", "STAR to hg38; reads soft-clipped at the 5' end, or by more "
                             "than 3 nt at the 3' end, are discarded"],
                     [f"call a small RNA",
                      f"{sm.SMALLRNA_MIN_NT}&ndash;{sm.SMALLRNA_MAX_NT} nt; reads at the "
                      f"maximum length ({sm.PRECURSOR_MIN_NT} bp for a {sm.READ_LEN}-bp "
                      f"run) are set aside as possible precursors"],
                     ["deduplicate",
                      "UMI-tools adjacency method, Hamming distance &ge; 2"],
                     ["assign",
                      "miRBase &rarr; GtRNAdb &rarr; Gencode, hierarchically; multimappers "
                      "split evenly across annotated loci"]])
            + info(f"""&#9989; The arithmetic checks out both ways: a {sm.READ_LEN}-bp read
            less {sm.TRIM_5P}&nbsp;nt leaves {sm.READ_LEN - sm.TRIM_5P}&nbsp;bp, and the
            paper's own worked example of a 76-bp read leaves {76 - sm.TRIM_5P}&nbsp;bp.""")
            + "</div>")


def limits() -> str:
    return ("<h2>What this chemistry cannot do</h2>\n"
            + table(["Limitation", "The paper's words"],
                    [["Ligation bias at both ends",
                      "&ldquo;Small-seq is, however, biased toward certain sequence "
                      "compositions at the ends of small RNAs&rdquo; &mdash; though "
                      "&ldquo;using UMIs with random sequences at the 5&prime; adapter may "
                      "reduce the bias at one end&rdquo;"],
                     ["Adapter dimers survive size selection",
                      "&ldquo;the abundance of adapter dimers in the sequencing reads, "
                      "despite the effects of the digestion step&rdquo;"],
                     ["Degradation products look like small RNAs",
                      "&ldquo;the products of RNA degradation can also be captured, mapped, "
                      "and annotated alongside the &lsquo;genuine&rsquo; small RNAs&rdquo;"],
                     ["It sequences the reagents",
                      "&ldquo;this approach is capable of sequencing the small amounts of "
                      "RNAs that are likely to contaminate the enzymes and reagents&rdquo; "
                      "&mdash; a no-cell control is mandatory, and its Bioanalyzer trace "
                      "&ldquo;looked almost identical&rdquo; to a single cell's"],
                     ["Plate scale only",
                      "192 samples per index set, one cell per well; the authors name the "
                      "fix themselves &mdash; &ldquo;adapt dual-sample indexing and "
                      "nano-well plates&rdquo;"],
                     [f"{sm.RA5_UMI_DIVERSITY:,} UMIs is not many",
                      "&ldquo;in principle, high enough to count all the molecules "
                      "transcribed from a single locus&rdquo; &mdash; an assumption, not a "
                      "result; saturation would bias counts downward"]]))


def main() -> None:
    graph = sm.protocol_graph()
    parts = [head("Small-seq Chemistry"), preamble().replace(
             '<div class="wrap">','<div class="wrap" data-protocol-complete="true">'),
             oligos(), interlock(), WorkflowSection("Reaction workflow"),
             final_library(), limits(), sequencing()]
    OUT.write_text(assemble_page(graph, parts), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
