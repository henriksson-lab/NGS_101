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
from illumina import P5
from page import caveat, head, info, legend, table

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

{caveat("""<b>One oligo in this chemistry is not published.</b> The <cbc>SRX</cbc> index
primers live in Supplementary Table 1, which is not in the extracted text, so everything
right of <r2>RA3</r2> is drawn as a placeholder and marked inferred (dotted underline).
Its <i>length</i> is bracketed to 37&ndash;47&nbsp;nt by two independent published sizes
&mdash; see <i>Sizes</i> below. Nothing on this page invents a base of it.""")}
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
              [seg("", "[not published &mdash; Supplementary Table 1]", None,
                   placeholder=True, inferred=True)]),
        oligo("5.8S rRNA masking oligo (DNA)", [seg("", sm.MASK_58S, "w1")],
              three="-biotin"),
    ]
    return ("<h2>Oligos</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>\n"
            + info(f"""All six verbatim from <i>Reagent setup</i>, which states
            &ldquo;All oligonucleotides are listed in the 5&prime; to 3&prime;
            direction&rdquo; &mdash; except <cbc>SRX</cbc>, which is not there at all.
            <code>r</code> marks a ribonucleotide; <code>rH</code> is rA, rU or rC;
            <code>rApp</code> is a chemically pre-adenylated 5' end.
            Lengths: <r1>RA5</r1> {len(sm.RA5_HANDLE)}&nbsp;+&nbsp;{sm.RA5_UMI_LEN}
            &nbsp;+&nbsp;{len(sm.RA5_LINKER)}&nbsp;=&nbsp;{len(sm.ra5_oligo())}&nbsp;nt,
            <r2>RA3</r2> {len(sm.RA3)}&nbsp;nt, RTP {len(sm.RTP)}&nbsp;nt,
            RP1 {len(sm.RP1)}&nbsp;nt, mask {len(sm.MASK_58S)}&nbsp;nt."""))


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
        caption="(8) PCR 2: the same forward primer plus one SRX index primer per well. "
                "13 cycles, annealing at 67 C. INFERRED -- SRX is not published: its 3' end "
                "must anneal over RA3 where RTP does, and its 5' tail (link + i7 + P7) is "
                "drawn as placeholder overhang.")
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
            + "<h3>Sizes, and how far the unpublished arm can be pinned down</h3>\n"
            + info(f"""<cbc>SRX</cbc>'s 3' end must anneal where RTP anneals &mdash; over
            <r2>RA3</r2> &mdash; so what it <i>adds</i> is P7 + index + some linker. Two
            independent published sizes bracket that length:""")
            + bracket
            + info(f"""The windows overlap, which they did not have to: <b>37&ndash;47
            nt</b>. A TruSeq-small-RNA-shaped arm &mdash; P7&nbsp;(24) + index&nbsp;(8) +
            a 12-nt link &mdash; is <b>{sm.SRX_RIGHT_ARM_NT}&nbsp;nt</b>, inside both.
            Taking that figure, <code>SRX_LINKER_NT</code> in
            <code>tools/smallseq.py</code>, every size below re-derives. Flip it and they
            all move together.""")
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
    quotes = [
        ["Overview of the protocol",
         "&ldquo;Finally, the libraries are pooled and purified, <b>ready for sequencing on "
         "an Illumina sequencer such as a HiSeq or NextSeq platform</b>.&rdquo;"],
        ["Experimental design",
         "&ldquo;<b>Sequencing is performed using Illumina platforms.</b> We normally "
         "sequence 1&ndash;2 million reads per cell.&rdquo;"],
        ["Experimental design",
         "&ldquo;<b>Sequencing starts from the UMI</b>, which is immediately followed by a "
         "&lsquo;-CA-&rsquo; couple, which separates the UMI from the small RNA "
         "sequences.&rdquo;"],
        ["Equipment",
         "&ldquo;Illumina DNA sequencer (e.g., HiSeq, MiSeq, NovaSeq, or NextSeq).&rdquo;"],
        ["Step 37",
         "&ldquo;Perform a pilot sequencing experiment with a pool of 48 samples, using a "
         "MiSeq or a NextSeq instrument. &hellip; <b>Carry out sequencing according to the "
         "manufacturer&rsquo;s protocol.</b>&rdquo;"],
        ["Step 38",
         "&ldquo;Perform sequencing for a <b>single read of 51 bp or longer</b>, depending "
         "on the small RNA molecules of interest. Sequencing of 51 bp is enough to detect "
         "miRNAs.&rdquo;"],
        ["Reagent setup, SRX",
         "&ldquo;SRX DNA index primers are <b>modified from standard Illumina TruSeq small "
         "RNA index primers</b> to have barcodes of 8 bp for 192 samples. For sequence "
         "information, see Supplementary Table 1.&rdquo;"],
        ["Fig. 1 legend",
         "&ldquo;Sequences required for Illumina cluster generation and sample indexing are "
         "added through two rounds of PCR.&rdquo;"],
    ]
    answers = table(
        ["Question", "Answer"],
        [["Does the paper name a sequencing primer?",
          "&#128308; <b>No.</b> The words &ldquo;sequencing primer&rdquo;, &ldquo;custom "
          "primer&rdquo;, &ldquo;primer mix&rdquo;, &ldquo;Read 1&rdquo;, &ldquo;flow "
          "cell&rdquo;, &ldquo;spike-in&rdquo; and &ldquo;PhiX&rdquo; occur <b>zero</b> "
          "times in the article."],
         ["Does it say a custom primer must be supplied?", "&#128308; <b>No.</b>"],
         ["Does it say the standard instrument primers suffice?",
          "&#128308; <b>No.</b> The closest it comes is &ldquo;ready for sequencing on an "
          "Illumina sequencer such as a HiSeq or NextSeq platform&rdquo;, &ldquo;relies on "
          "standard reagents and instruments&rdquo;, and &ldquo;Carry out sequencing "
          "according to the manufacturer&rsquo;s protocol&rdquo;."],
         ["Can a TruSeq <b>DNA</b> read-1 primer, or RP1 itself?",
          "&#9989; <b>No</b> &mdash; computed in <i>Sequencing primers</i> below."],
         ["Can a TruSeq <b>Small RNA</b> read-1 primer?",
          "&#128993; <b>Yes</b>, on the structure: it anneals to the same 26-nt "
          "<r1>RA5</r1> region and simply emits 10 extra bases &mdash; the "
          "<umi>UMI</umi> + <tso>CA</tso> &mdash; before the insert, which is exactly what "
          "the pipeline trims."],
         ["Which primer did the authors actually use?",
          "&#128308; <b>Unknown.</b> Possibly in Supplementary Table 1, which is not in "
          "the extracted text."]])
    return ("<h2>Sequencing &mdash; and which primer reads this library</h2>\n"
            + panel(rows, cls="long",
                    caption="Read 1. Its 3' end is fixed by the construct: sequencing "
                            "starts at the UMI, so the primer must reach the last base of "
                            "RA5 and no further.")
            + info(f"""<b>The construct fixes the primer's 3' end even though the paper
            does not name it.</b> The <umi>UMI</umi> starts at offset
            <b>{lib.offset("UMI")}</b>, immediately after P5 + RA5, so any usable read-1
            primer is a <i>suffix of <p5>P5</p5>&nbsp;+&nbsp;<r1>RA5</r1></i> ending on
            RA5's last base. Everything upstream of base {lib.offset("UMI") + 1} is those
            two sequences and nothing else.""")
            + caveat("""<b>There is no TruSeq DNA adapter anywhere in this library</b>
            (computed in <i>Sequencing primers</i> below), so a run configured only with the
            standard TruSeq DNA primers would yield nothing from it. That is true of any
            TruSeq Small RNA library, not a Small-seq peculiarity &mdash; but the paper never
            mentions it.""")
            + "<h3>Everything the paper says about sequencing, verbatim</h3>\n"
            + table(["Where", "What it says"], quotes)
            + "<h3>The answer, stated precisely</h3>\n" + answers
            + info(f"""<b>So:</b> Small-seq does not change the read-1 landing site relative
            to a stock TruSeq Small RNA library &mdash; it only pushes the insert
            {sm.TRIM_5P}&nbsp;nt further into the read. Whatever primer arrangement a
            facility uses for TruSeq Small RNA libraries applies unchanged. The paper
            supports the construct half of that claim and is silent on the primer half.
            To close it: get Supplementary Table 1, and add a verbatim TruSeq Small RNA
            block to <code>lib/illumina.py</code> so the equality can be asserted rather
            than inferred.""")
            + sp.section(lib, sm.SEQ_PRIMERS, intro=f"""Small-seq is a single-read,
            single-indexed library. The TruSeq <i>DNA</i> set is listed because it is what a
            stock run loads; each has no site, for the reason given. The read-1 entry is the
            site the construct requires, not a vendor oligo &mdash; the TruSeq Small RNA
            primers are not yet in <code>lib/</code>. Read 1 reports <umi>UMI</umi> +
            <tso>CA</tso> before the insert, the {sm.TRIM_5P}&nbsp;nt the pipeline trims.""")
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
    parts = [head("Small-seq Chemistry"), preamble(), oligos(), interlock(),
             step_mask(), step_lig3(), step_digest(), step_lig5(), step_rt_pcr(),
             final_library(), limits(), sequencing()]
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
