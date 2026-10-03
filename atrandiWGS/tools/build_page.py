#!/usr/bin/env python3
"""Generate atrandi_wgs.html -- the scg_lib_structs-style chemistry page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import atrandi as A
from chemdraw import (Construct, Row, Segment, annotation_rows, complement, oligo,
                      panel, revcomp, strand_row, tm)

OUT = HERE.parent / "atrandi_wgs.html"


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


def chunks(*cs):
    """(text, tag, inferred) tuples, tolerating 1- and 2-tuples."""
    out = []
    for c in cs:
        t, g, i = (list(c) + [None, False])[:3]
        out.append((t, g, i))
    return out


def row(*cs, indent=0, prefix="", suffix=""):
    return Row(chunks=chunks(*cs), indent=indent, prefix=prefix, suffix=suffix)


# ============================================================== shared sequence objects
LIG_TOP = [seg("", A.LIGADAPT_STEM, "s5"), seg("", A.LIGADAPT_ARM, "s5")]
LIG_BOT = [seg("", A.LIGADAPT_BOT, "s5")]
I7_PRIMER = [seg("", A.P7_SEQ, "p7"), seg("", A.ATRANDI_I7_INDEX),
             seg("", A.ATRANDI_I7_ANNEAL, "t7")]
# NB: P5_SEQ and the Read-1 annealing portion overlap by "ACAC" -- do not concatenate them,
# split the real primer instead, so the rendered line always equals ATRANDI_P5 exactly.
P5_PRIMER = [seg("", A.P5_SEQ, "p5"), seg("", A.ATRANDI_P5[len(A.P5_SEQ):], "s5")]

lib = A.final_library()
PRE = 4  # width of the "5'- " prefix


# =========================================================================== page parts
def head() -> str:
    return """<title>Atrandi PTA Library Chemistry</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@400;500;600&display=swap">
<style>
/* ---- tokens: complete light palette on bare :root ------------------------ */
:root {
  --bg:#fbfcfd; --surface:#f2f6fa; --surface-2:#e9eff5;
  --ink:#16202b; --ink-muted:#566879; --rule:#dce4ed;
  --link:#0b5fa5; --accent:#0b5fa5;
  --note-bg:#fff6e3; --note-rule:#d79a16;
  --key-bg:#eef3f8; --key-rule:#7f95aa;
  --c-p5:#08519c; --c-p7:#a50f15; --c-s5:#3d87bd; --c-s7:#e2663f;
  --c-me:#7d7d7d; --c-t7:#1f4fd8; --c-cbc:#d6417c; --c-umi:#6f58b5;
  --c-r1:#4a1486; --c-r2:#6a51a3; --c-r3:#807dba; --c-tso:#1f8a4d; --c-w1:#d6341a;
  --inf:#8a98a6;
  --mono:'IBM Plex Mono','DejaVu Sans Mono',Menlo,Consolas,monospace;
  --sans:'IBM Plex Sans','Helvetica Neue',Helvetica,Arial,sans-serif;
}
/* ---- same tokens, dark values, for system-dark (un-stamped) -------------- */
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg:#0e1419; --surface:#161e26; --surface-2:#1d272f;
    --ink:#e2e9ef; --ink-muted:#93a4b4; --rule:#28333d;
    --link:#78b4ee; --accent:#78b4ee;
    --note-bg:#2a2214; --note-rule:#b98a1f;
    --key-bg:#17202a; --key-rule:#5a7389;
    --c-p5:#6ba6e2; --c-p7:#f0837f; --c-s5:#79bde4; --c-s7:#f59b78;
    --c-me:#9ba3ab; --c-t7:#8aa5ff; --c-cbc:#ff8ab8; --c-umi:#b0a2e8;
    --c-r1:#bda4ec; --c-r2:#b9a6e6; --c-r3:#a99ee0; --c-tso:#4fc98a; --c-w1:#ff8360;
    --inf:#78868f;
  }
}
/* ---- and again for an explicit dark choice ------------------------------ */
:root[data-theme="dark"] {
  --bg:#0e1419; --surface:#161e26; --surface-2:#1d272f;
  --ink:#e2e9ef; --ink-muted:#93a4b4; --rule:#28333d;
  --link:#78b4ee; --accent:#78b4ee;
  --note-bg:#2a2214; --note-rule:#b98a1f;
  --key-bg:#17202a; --key-rule:#5a7389;
  --c-p5:#6ba6e2; --c-p7:#f0837f; --c-s5:#79bde4; --c-s7:#f59b78;
  --c-me:#9ba3ab; --c-t7:#8aa5ff; --c-cbc:#ff8ab8; --c-umi:#b0a2e8;
  --c-r1:#bda4ec; --c-r2:#b9a6e6; --c-r3:#a99ee0; --c-tso:#4fc98a; --c-w1:#ff8360;
  --inf:#78868f;
}

body { background:var(--bg); color:var(--ink); font-family:var(--sans);
       margin:0; padding-block:32px 96px; padding-inline:16px; line-height:1.5; }
.wrap { max-width:1180px; margin:0 auto; display:flex; flex-direction:column; gap:4px; }

h1 { font-size:clamp(1.5rem,1.1rem + 1.6vw,2.1rem); font-weight:600; line-height:1.2;
     text-wrap:balance; margin:0 0 .2em; letter-spacing:-0.01em; }
h2 { font-size:1.3rem; font-weight:600; text-wrap:balance; margin:2.2em 0 .2em;
     padding-bottom:.3em; border-bottom:2px solid var(--rule); }
h3 { font-size:1rem; font-weight:600; text-wrap:balance; margin:2em 0 .3em;
     color:var(--ink); }
h3::before { content:''; display:block; width:28px; height:2px; background:var(--accent);
             margin-bottom:.6em; opacity:.55; }
p { margin:.5em 0; }
info { display:block; max-width:68ch; color:var(--ink); }
a { color:var(--link); text-underline-offset:2px; }
a:focus-visible, :focus-visible { outline:2px solid var(--accent); outline-offset:2px; }
code { font-family:var(--mono); font-size:.9em; background:var(--surface-2);
       padding:.08em .3em; border-radius:3px; }

seq { font-family:var(--mono); font-size:.86rem; display:block; margin:.4em 0 0; }
seq p { margin:.3em 0; text-indent:-2.9em; padding-left:2.9em; }

pre { overflow-x:auto; margin:.6em 0 .2em; background:var(--surface);
      border-left:3px solid var(--rule); padding:12px 14px; }
align { font-family:var(--mono); display:block; white-space:pre; }
align.long  { font-size:.74rem; line-height:1.25; }
align.small { font-size:.84rem; line-height:1.3; }
align i { color:var(--ink-muted); display:block; font-style:italic; line-height:1.5;
          margin-bottom:.7em; max-width:96ch; white-space:pre-wrap; }

.legend, .caveat { padding:12px 16px; margin:1em 0; max-width:88ch; font-size:.95rem; }
.legend { background:var(--key-bg); border-left:3px solid var(--key-rule); }
.caveat { background:var(--note-bg); border-left:3px solid var(--note-rule); }

.tw { overflow-x:auto; }
table { border-collapse:collapse; margin:.8em 0; font-size:.92rem;
        font-variant-numeric:tabular-nums; }
th, td { border:1px solid var(--rule); padding:5px 12px; text-align:left; }
th { background:var(--surface); font-weight:600; }

/* region colours: what a stretch of sequence IS */
p5{color:var(--c-p5)} p7{color:var(--c-p7)} s5{color:var(--c-s5)} s7{color:var(--c-s7)}
me{color:var(--c-me)} t7{color:var(--c-t7)} cbc{color:var(--c-cbc)} umi{color:var(--c-umi)}
r1{color:var(--c-r1)} r2{color:var(--c-r2)} r3{color:var(--c-r3)} tso{color:var(--c-tso)}
w1{color:var(--c-w1)}
/* and how well it is KNOWN -- orthogonal, nests freely inside the above */
inf { border-bottom:1px dotted var(--inf); }
</style>"""


def preamble() -> str:
    return f"""<div class="wrap">
<h1>Atrandi semi-permeable capsules + PTA &mdash; single-microbe whole-genome sequencing</h1>

<p><info>This page documents the molecular structure of the single-cell microbial WGS libraries
produced by the protocol in <a href="https://doi.org/10.1101/2025.06.20.660799">Gourl&eacute; et al.,
"Scalable single-cell metagenomic analysis with Bascet and Zorn", bioRxiv 2025.06.20.660799</a>.
The workflow is built on the <b>Atrandi Biosciences Single-Microbe DNA Barcoding Kit</b>
(P/N CKP-BARK1) and its <i>Library Prep for Sequencing</i> guide, with two deliberate changes:
whole-genome amplification is <b>PTA</b> (BioSkryb ResolveDNA, #100954) rather than the kit's MDA,
and the <b>debranching step is therefore omitted</b>.</info></p>

<p><info>Cells are encapsulated in <b>semi-permeable capsules</b> (SPCs) &mdash; a picolitre aqueous
core inside a hydrogel shell that lets enzymes, detergents and salts diffuse freely while retaining
genomic DNA. That is what makes harsh lysis, repeated buffer exchange and four rounds of split-pool
ligation possible without ever isolating a single cell. Four barcode rounds of 24 variants each give
24<sup>4</sup> = <b>331,776</b> combinations.</info></p>

<div class="legend">
<b>How to read this page.</b> Colour says what a region <i>is</i>; a
<inf>dotted underline</inf> says how well it is <i>known</i>. Dotted regions are
<b>inferred</b>, not documented, and each inferred step carries a caption saying what the inference
rests on. Barcodes and linkers are drawn as placeholders
(<cbc>AAAAAAAA</cbc>&hellip;<cbc>DDDDDDDD</cbc>, <r2>LLLL</r2>) because the real sequences do not
affect the architecture.
</div>

<div class="caveat">
<b>What is <i>not</i> documented anywhere.</b> Atrandi publish the sequences of their ligation
adapter and indexing primers, but <b>not the barcode cassette</b> &mdash; no strand design, no
overhangs, no linker sequences. Since the ligation adapter demonstrably installs only the
Read&nbsp;1 / P5 side (step 8), the Read&nbsp;2 / P7 arm <i>must</i> come from the barcode cassette,
and steps 4&ndash;5 below are a <b>model</b> consistent with everything we can check &mdash; not
published chemistry. The two sequences it predicts are marked.
</div>
"""


def oligos() -> str:
    tmp5 = tm(A.ATRANDI_P5_ANNEAL)
    tmp7 = tm(A.ATRANDI_I7_ANNEAL)
    rows = [
        oligo("Exo-resistant random primer (PTA)",
              [seg("", "NpNpNpNpsNpsN", "tso", placeholder=True)],
              mods=""),
        oligo("Terminator (PTA)", [seg("", "alpha-thio-ddNTP", "w1", placeholder=True)], five="", three=""),
        oligo("Ligation Adapter, top strand (CRP-LGA1)", LIG_TOP, mods="/5Phos/"),
        oligo("Ligation Adapter, bottom strand", LIG_BOT, mods="/5AmMC6/"),
        oligo("Barcode cassette, rounds A-D (CRP-PLT1)",
              [seg("", "[8-bp barcode]", "cbc", placeholder=True),
               seg("", "[4-nt linker]", "r2", placeholder=True, inferred=True)],
              five="", three=""),
        oligo("Indexed Library PCR Primer 1 (i7)", I7_PRIMER),
        oligo("Library PCR Primer 2 (P5)", P5_PRIMER),
        oligo("Read 1 sequencing primer", [seg("", A.TRUSEQ_R1, "s5")]),
        oligo("Index 1 sequencing primer", [seg("", "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "t7")]),
        oligo("Read 2 sequencing primer", [seg("", A.TRUSEQ_R2, "t7")]),
    ]
    return f"""<h2>Adapter and primer sequences</h2>
<seq>
{chr(10).join(rows)}
</seq>
<p><info>The <code>*</code> in Atrandi's published primers is a 3'-terminal phosphorothioate,
protecting against Q5's proofreading exonuclease. Note both PCR primers are <b>truncated</b>
relative to canonical Illumina: the P5 primer's 3' annealing portion is exactly
{len(A.ATRANDI_P5_ANNEAL)}&nbsp;nt, the same length as the ligation adapter's single-stranded arm,
with zero slack. That sizing is what sets the 54&nbsp;&deg;C annealing temperature
(T<sub>m</sub> {tmp5:.1f}&nbsp;&deg;C for P5 vs {tmp7:.1f}&nbsp;&deg;C for i7) &mdash; and why NEB's
own primers cannot be substituted.</info></p>
"""


# ------------------------------------------------------------------------------ steps
def step_pta() -> str:
    amp_top = chunks(("NpNpNpNpsNpsN", "tso"), ("XXXXXXXXXXXX...XXXXXXXXXXXX",), ("ddN", "w1"))
    amp_bot = chunks(("Ndd", "w1"), ("XXXXXXXXXXXX...XXXXXXXXXXXX",), ("NsNsNpNpNpN", "tso"))
    rows = [
        row(("XXXXXXXXXXXXXXXXXXXXXXXX...XXXXXXXXXXXXXXXXXXXXXXXX",), prefix="5'- ", suffix=" -3'"),
        row(("XXXXXXXXXXXXXXXXXXXXXXXX...XXXXXXXXXXXXXXXXXXXXXXXX",), prefix="3'- ", suffix=" -5'"),
        Row(chunks=chunks(("    ",), ("NNNNNN", "tso"), ("------------>",)), indent=4),
        Row(chunks=chunks(("        phi29 extends, displacing downstream strand,",)), indent=4),
        Row(chunks=chunks(("        until an alpha-thio-ddN is incorporated",)), indent=4),
    ]
    p1 = panel(rows, caption="Random primers anneal across the genome; phi29 extends with strand "
                             "displacement.")
    p2 = panel([Row(chunks=amp_top, prefix="5'- ", suffix=" -3'"),
                Row(chunks=amp_bot, prefix="3'- ", suffix=" -5'")],
               caption="The resulting amplicon. Both termini are chemically distinctive and both "
                       "matter downstream.")
    return f"""<h3>(2) Whole-genome amplification by PTA (BioSkryb ResolveDNA #100954; 30&nbsp;&deg;C
2.5&nbsp;h &rarr; 65&nbsp;&deg;C 5&nbsp;min)</h3>
{p1}
{p2}
<p><info><b>5' end</b> &mdash; the exo-resistant random primer, <tso>5'-NpNpNpNpsNpsN-3'</tso>: six
random bases whose last two internucleotide linkages are phosphorothioate (Dean <i>et al.</i> 2002;
6&ndash;9mer in the commercial kit). Synthetic oligos carry a <b>5'-hydroxyl</b>, so these ends are
<b>not ligatable until kinased</b> &mdash; which is why end-prep cannot be skipped.</info></p>
<p><info><b>3' end</b> &mdash; an <w1>alpha-thio-dideoxynucleotide</w1>, i.e. a
2',3'-dideoxyribonucleoside 5'-O-(1-thiotriphosphate). Two separate features at two different atoms,
and both are needed: the <b>2',3'-dideoxy sugar</b> blocks extension, and the <b>alpha-thio group</b>
makes the new linkage a phosphorothioate that resists phi29's 3'&rarr;5' proofreading exonuclease.
With plain ddNTPs the polymerase simply removes them and repriming generates chimeras &mdash; mapping
rates were 15.0&nbsp;&plusmn;&nbsp;2&nbsp;%; with the alpha-thio group, 97.9&nbsp;&plusmn;&nbsp;0.6&nbsp;%
(<a href="https://doi.org/10.1073/pnas.2024176118">Gonzalez-Pena <i>et al.</i> 2021</a>).</info></p>
<div class="caveat">
<b>No debranching step.</b> Terminating extension after ~250&ndash;2000&nbsp;bp means PTA never
builds MDA's hyperbranched network, so the Atrandi kit's debranching enzyme (37&nbsp;&deg;C, 1&nbsp;h)
is dropped. <a href="https://doi.org/10.1101/2025.09.10.675331">Negreira <i>et al.</i> 2025</a>,
running PTA in the same Atrandi SPC workflow, state it directly: the PTA samples
<i>"did not require debranching and were directly submitted to barcoding instead."</i>
End-prep is still required.
<br><br>
<b>Amplicon size vs capsule retention.</b> PTA products centre near 1.3&nbsp;kb in human cells and
<b>~900&nbsp;bp in bacteria</b>, against an SPC retention cutoff of &gt;500&nbsp;bp.
<a href="https://doi.org/10.1101/2025.03.14.643253">Mullaney <i>et al.</i> 2025</a> observed
150&ndash;1500&nbsp;bp amplified DNA leaking into the supernatant during a 30&nbsp;&deg;C PTA
reaction, and entering empty capsules. This is the central unquantified risk of PTA-in-capsules.
</div>
"""


def step_endprep() -> str:
    rows = [
        row(("XXXXXXXXXXXX...XXXXXXXXXXXX",), ("A",), prefix="5'-p ", suffix=" -3'"),
        row(("A",), ("XXXXXXXXXXXX...XXXXXXXXXXXX",), prefix=" 3'- ", suffix=" p-5'"),
    ]
    dead = [
        row(("XXXXXXXXXXXX...XXXXXXXXXXXX",), ("ddN", "w1"), prefix="5'- ", suffix=" -3'  <-- cannot be A-tailed"),
        row(("Ndd", "w1"), ("XXXXXXXXXXXX...XXXXXXXXXXXX",), prefix="3'- ", suffix=" -5'"),
    ]
    return f"""<h3>(3) End preparation (20&nbsp;&deg;C 30&nbsp;min &rarr; 65&nbsp;&deg;C
30&nbsp;min): blunting, 5'-phosphorylation and dA-tailing</h3>
{panel(rows, caption="Ends carrying a normal 3'-OH are blunted, kinased and dA-tailed, ready for "
                     "T-overhang ligation.")}
{panel(dead, caption="Ends terminating in the alpha-thio-dideoxynucleotide are not: no 3'-OH to "
                     "extend, and the phosphorothioate linkage resists the exonuclease that would "
                     "otherwise trim it back.")}
<div class="caveat">
<b>An open question, shown rather than smoothed over.</b> The PTA patents list "removing the
terminator" as a workflow step (US11643682B2 Fig.&nbsp;1D; WO2019148119A1 claim&nbsp;80) but
<b>name no enzyme</b> &mdash; and the same specification states elsewhere that irreversible
terminators are <i>"not capable of substantial removal by an exonuclease"</i>. The PNAS paper drops
the clause entirely and says amplicons <i>"undergo direct ligation of adapters"</i>.
<br><br>
<i>Most likely resolution (inferred):</i> barcoding happens at the ends that <b>do</b> carry a normal
3'-OH &mdash; from polymerase dissociation before terminator incorporation, from nicks, and from
fill-in at recessed ends. Each ~1&nbsp;kb amplicon has two ends and only needs one.
</div>
"""


def step_barcode_a() -> str:
    cassette = [
        row(("    ",), ("LLLL", "r2", True), ("AAAAAAAA", "cbc"), ("T",),
            prefix="5'- ", suffix=" -3'"),
        row(("        ",), ("aaaaaaaa", "cbc"),
            prefix="3'- ", suffix="  p-5'"),
        Row(),
        row(("    ^^^^        ^",), indent=4),
        row(("    |           3'-T overhang, joins the dA-tailed amplicon",), indent=4),
        row(("    4-nt overhang left for round B",), indent=4),
    ]
    joined = [
        row(("LLLL", "r2", True), ("AAAAAAAA", "cbc"), ("T",), ("XXXXXXXX...XXXXXXXX",),
            prefix="5'- ", suffix=" -3'"),
        row(("    ",), ("aaaaaaaa", "cbc"), ("A",), ("xxxxxxxx...xxxxxxxx",),
            prefix="3'- ", suffix=" -5'"),
    ]
    return f"""<h3>(4) Round A: split into 24 wells and ligate the first barcode</h3>
{panel(cassette, cls="small",
       caption="INFERRED -- Atrandi do not publish the barcode cassette design. To join the "
               "dA-tailed amplicon it must present a 3'-T overhang, and to accept round B it must "
               "leave a fresh overhang behind. Barcodes and linkers are placeholders.")}
{panel(joined, cls="small",
       caption="After ligation. The T/A junction base is the '+1' the demultiplexer trims "
               "(8+4+8+4+8+4+8+1 = 45). Lowercase marks the complement of a placeholder.")}
"""


def step_barcode_bcd() -> str:
    cassette = [
        row(("    ",), ("LLLL", "r2", True), ("BBBBBBBB", "cbc"),
            prefix="5'- ", suffix=" -3'"),
        row(("        ",), ("bbbbbbbb", "cbc"), ("llll", "r2", True),
            prefix="3'- ", suffix=" p-5'"),
        Row(),
        row(("    ^^^^        ^^^^",), indent=4),
        row(("    |           anneals to the growing molecule",), indent=4),
        row(("    overhang for the next round",), indent=4),
    ]
    blk = A.final_library()
    block = Construct(blk.segments[blk.segments.index(blk.get("BC-D")):-2])
    done = [strand_row(block, "top"), strand_row(block, "bottom")]
    return f"""<h3>(5) Rounds B, C and D: three more 24-way splits, each adding a barcode</h3>
{panel(cassette, cls="small",
       caption="INFERRED -- rounds B, C and D join through a 4-nt cohesive overhang and leave "
               "another behind. All 24 barcodes within a round share one linker, which is why four "
               "barcodes are separated by exactly three linkers.")}
{panel(done, cls="small", caption="After all four rounds.")}
<p><info>Ligation order is <b>A &rarr; B &rarr; C &rarr; D</b>, so D ends up outermost and is read
<b>first</b> in Read&nbsp;2. The round-D cassette must also carry the <t7>Read&nbsp;2 / P7 arm</t7>,
because &mdash; as step&nbsp;8 shows &mdash; the ligation adapter does not. That end must be
ligation-dead, or it would pick up a Read&nbsp;1 arm during library prep and the molecule would
carry P5 at both ends.</info></p>
<div class="caveat">
<b>The model reproduces an independent fact.</b> "Four 8-nt barcodes, three 4-nt cohesive junctions,
one TA ligation base" predicts a Read&nbsp;2 layout of <b>8+4+8+4+8+4+8+1</b>, with barcodes anchored
at offsets <b>0, 12, 24, 36</b> and <b>45&nbsp;nt</b> trimmed. Those are exactly the constants
hard-coded in <a href="https://github.com/henriksson-lab/bascet">Bascet</a>, derived independently
from sequencing data. This page's diagrams are generated from the segment table, and a self-test
asserts the agreement.
</div>
"""


def step_release_frag() -> str:
    frag = [
        row(("LLLL", "r2", True), ("AAAAAAAA", "cbc"), ("T",), ("XXXXXXXX...XXXXXXXX",),
            ("  |  ",), ("XXXXXX...XXXXXX",), prefix="5'- ", suffix=" -3'"),
        row(("    ",), ("aaaaaaaa", "cbc"), ("A",), ("xxxxxxxx...xxxxxxxx",),
            ("  |  ",), ("xxxxxx...xxxxxx",), prefix="3'- ", suffix=" -5'"),
        Row(),
        row(("^ the FS enzyme cuts here; only the barcode-proximal fragment",), indent=38),
        row(("  becomes a library molecule",), indent=38),
    ]
    return f"""<h3>(6) Release from the capsules, then (7) enzymatic fragmentation
(NEBNext Ultra II FS, #E7805S; 10&nbsp;min at 37&nbsp;&deg;C &rarr; 300&ndash;700&nbsp;bp)</h3>
{panel(frag, cls="small",
       caption="The Release reagent dissolves the hydrogel shell; the FS enzyme mix then "
               "fragments, end-repairs and dA-tails in one tube.")}
"""


def step_adapter() -> str:
    from chemdraw import complement
    ad = [
        row((A.LIGADAPT_STEM, "s5"), (A.LIGADAPT_ARM, "s5"),
            prefix="5'-/5Phos/ ", suffix=" -3'"),
        row(("T",), (complement(A.LIGADAPT_STEM), "s5"),
            prefix="      3'- ", suffix=" /5AmMC6/-5'"),
        Row(),
        row(("^ 3'-T overhang (the ligation end)",), indent=10),
        row(("^" + " " * 11 + "^",), indent=11),
        row(("12-bp stem  20-nt single-stranded arm = the P5 primer landing site",), indent=11),
    ]
    both = A.adapter_both_ends()
    return f"""<h3>(8) Ligate the Atrandi adapter (20&nbsp;&deg;C, 15&nbsp;min &mdash; no USER step)</h3>
{panel(ad, cls="small",
       caption="The Atrandi ligation adapter. This is NOT a forked/Y adapter: the 13-nt bottom "
               "strand pairs over all 12 of its 5'-proximal bases, leaving one arm, not two.")}
<p><info>Verified base-by-base: <code>revcomp(GCTCTTCCGATC) = GATCGGAAGAGC</code>, so the stem is
exactly 12&nbsp;bp; <code>revcomp(GTCGTGTAGGGAAAGAGTGT)</code> is the <b>first 20&nbsp;nt of the
TruSeq Read&nbsp;1 primer</b>; and the bottom strand is the <b>last 13&nbsp;nt</b> of that same
primer. After ligation the strand reads <code>AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT</code> &mdash;
character for character, Illumina's canonical Read&nbsp;2 adapter-trimming sequence, the leading A
being the dA tail.</info></p>
<p><info>The <code>/5AmMC6/</code> amino block sits at the fork point and makes that 5' end
permanently ligation-incompetent. Unlike a bare 5'-OH it cannot be rescued by the kinase activity
carried over from the FS end-repair mix, so adapter&ndash;adapter dimers cannot form.</info></p>
{panel([strand_row(both, "top"), strand_row(both, "bottom")], cls="small",
       caption="A fragment that received the adapter but NO barcode: the same Read-1 / P5 site at "
               "both ends.")}
<div class="caveat">
<b>Why only barcoded molecules amplify.</b> The adapter installs <b>only</b> the Read&nbsp;1 / P5
side, at every end. A fragment lacking the barcode cassette therefore has a P5 landing site at both
ends; after one P5 extension its new 3' end is <code>AGATCGGAAGAGC</code>, for which no primer
exists, so it amplifies <b>linearly</b>. Only molecules carrying the barcode cassette &mdash; and
hence the i7 landing site &mdash; go exponential. This is a suppression PCR.
</div>
"""


def step_pcr() -> str:
    pp = A.pre_pcr()
    IND = 30                       # room for the i7 primer's P7 + index tail
    top = strand_row(pp, "top", indent=IND)
    bot = strand_row(pp, "bottom", indent=IND)
    i7 = Row(chunks=chunks(*[(s.top, s.tag, s.inferred) for s in I7_PRIMER], ("------->",)),
             prefix="5'- ")
    arm_col = pp.offset("TruSeq Read 1") + len(A.LIGADAPT_STEM) + IND + PRE
    p5_rev = A.ATRANDI_P5[::-1]
    p5 = Row(chunks=chunks(("<------",), (p5_rev[:len(A.LIGADAPT_ARM)], "s5"),
                           (p5_rev[len(A.LIGADAPT_ARM):], "p5")),
             indent=arm_col - 7, suffix=" -5'")
    return f"""<h3>(9) Indexing PCR (98&nbsp;&deg;C 20&nbsp;s / 54&nbsp;&deg;C 30&nbsp;s /
72&nbsp;&deg;C 20&nbsp;s, 8&ndash;14 cycles)</h3>
{panel([i7, top, bot, p5],
       caption="The i7 primer lands on the Read-2 arm supplied by the round-D cassette, its P7 and "
               "index tail unpaired. The P5 primer lands on the adapter's 20-nt arm with zero "
               "slack -- its 25-nt P5 tail hangs off the end and is copied in on the next cycle.")}
<p><info>Both primers are truncated to match their landing sites exactly, which is why
Atrandi's own primers are required and NEB's cannot be substituted: a full-length Illumina P5 primer
would have 13 unpaired 3' bases here and could not prime at all. The 20-nt P5 arm also sets the
annealing temperature &mdash; T<sub>m</sub> {tm(A.ATRANDI_P5_ANNEAL):.1f}&nbsp;&deg;C, hence
54&nbsp;&deg;C rather than NEB's 65&nbsp;&deg;C.</info></p>
"""


def step_final() -> str:
    rows = [strand_row(lib, "top"), strand_row(lib, "bottom")] + annotation_rows(lib, prefix_width=PRE)
    return f"""<h3>(10) Final library structure</h3>
{panel(rows)}
<p><info>Total {len(lib)}&nbsp;bp excluding the insert. There is <b>no UMI</b> in this chemistry, and
the i5 position carries no index &mdash; Atrandi's P5 primer is the plain universal primer, so i5 is
optional.</info></p>
<div class="caveat">
<b>One inferred length.</b> The round-D Read&nbsp;2 arm is drawn as
{len(A.ATRANDI_I7_ANNEAL)}&nbsp;nt, the exact length of the i7 primer's annealing portion, following
the design rule proven on the P5 side. It may instead be the full 34&nbsp;nt canonical arm, with the
primer landing 7&nbsp;nt further out &mdash; the primer-binding site is identical either way, so only
those 7&nbsp;nt (<t7>CCGATCT</t7>) are uncertain.
</div>
"""


def sequencing() -> str:
    from chemdraw import complement
    lib_top = strand_row(lib, "top")
    lib_bot = strand_row(lib, "bottom")
    r1_start = lib.offset("TruSeq Read 1")
    r2_start = lib.offset("BC-D")
    i7_start = lib.offset("TruSeq Read 2")

    r1_seq = A.TRUSEQ_R1[::-1]                       # 3'->5' left to right
    r1 = Row(chunks=chunks(("<--------------",), (r1_seq, "s5")),
             indent=r1_start - 1 + PRE - 15, suffix=" -5'")
    i7_seq = complement(A.ATRANDI_I7_ANNEAL)         # pairs with the 27-nt arm
    i7 = Row(chunks=chunks(("<------",), (i7_seq, "t7")),
             indent=i7_start + PRE - 7, suffix=" -5'")
    r2 = Row(chunks=chunks((A.TRUSEQ_R2, "t7"), ("------------------->",)),
             indent=max(0, r2_start + PRE - len(A.TRUSEQ_R2) - 4), prefix="5'- ")
    return f"""<h2>Library sequencing</h2>
<p><info>Paired-end on an Illumina NovaSeq&nbsp;X. Atrandi's guide specifies Read&nbsp;1
128&nbsp;bp (genomic insert) and Read&nbsp;2 172&nbsp;bp (cell barcode + insert); the read lengths
actually used in the preprint are not stated, but Read&nbsp;2 must exceed 45&nbsp;nt to contain all
four barcodes.</info></p>

<h3>(1) Read 1 &mdash; genomic insert, primed from the P5 side (top strand as template)</h3>
{panel([lib_top, lib_bot, r1])}

<h3>(2) Index 1 read &mdash; the 6-bp i7 sample index (top strand as template)</h3>
{panel([lib_top, lib_bot, i7])}
<p><info>Only the 27&nbsp;nt that pair are drawn; the standard index primer's 5'-terminal
<code>GATCGG</code> has nothing to anneal to and hangs off. The index is reported as the reverse
complement of the bases in the PCR primer: the oligo carries <code>{A.ATRANDI_I7_INDEX}</code>, so
the index read gives <code>{revcomp(A.ATRANDI_I7_INDEX)}</code>.</info></p>

<h3>(3) Read 2 &mdash; cell barcode then insert (bottom strand as template)</h3>
{panel([r2, lib_top, lib_bot])}
<p><info>Read&nbsp;2 reports the top-strand sequence: barcode <cbc>D</cbc> first, then
<cbc>C</cbc>, <cbc>B</cbc>, <cbc>A</cbc>, then the insert. The primer's 5' end extends past the
27-nt arm and is unpaired there.</info></p>

<h2>Read 2 layout</h2>
<div class="tw">
<table>
<tr><th>Offset (0-based)</th><th>Length</th><th>Content</th></tr>
<tr><td>0&ndash;7</td><td>8</td><td>Barcode D (ligated last, read first)</td></tr>
<tr><td>8&ndash;11</td><td>4</td><td>linker <i>(inferred)</i></td></tr>
<tr><td>12&ndash;19</td><td>8</td><td>Barcode C</td></tr>
<tr><td>20&ndash;23</td><td>4</td><td>linker <i>(inferred)</i></td></tr>
<tr><td>24&ndash;31</td><td>8</td><td>Barcode B</td></tr>
<tr><td>32&ndash;35</td><td>4</td><td>linker <i>(inferred)</i></td></tr>
<tr><td>36&ndash;43</td><td>8</td><td>Barcode A (ligated first, read last)</td></tr>
<tr><td>44</td><td>1</td><td>dA/dT ligation junction <i>(inferred)</i></td></tr>
<tr><td>45&ndash;</td><td>&mdash;</td><td>genomic insert</td></tr>
</table>
</div>
<p><info>Demultiplexing uses a cascading match (D, then C, then B, then A) allowing one mismatch per
barcode and four in total, with the anchors re-derived from the first 10,000 reads.
Barcode&nbsp;A is checked first, being the most likely to be lost to over-fragmentation.</info></p>
</div>
"""


def main() -> None:
    parts = [head(), preamble(), oligos(),
             "<h2>Step-by-step library generation</h2>",
             "<h3>(1) Encapsulate single cells in SPCs and lyse</h3>",
             "<p><info>Cells are diluted to a Poisson loading of &lambda;&nbsp;&le;&nbsp;0.1 and "
             "encapsulated on the Atrandi ONYX. Lysis is a four-enzyme cocktail (lysostaphin, "
             "mutanolysin, lysozyme, achromopeptidase) overnight at 37&nbsp;&deg;C, then "
             "proteinase&nbsp;K at 55&nbsp;&deg;C, then 4&nbsp;M guanidine thiocyanate &mdash; all "
             "performed through the capsule shell, which retains the genomic DNA throughout.</info></p>",
             step_pta(), step_endprep(), step_barcode_a(), step_barcode_bcd(),
             step_release_frag(), step_adapter(), step_pcr(), step_final(), sequencing()]
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
