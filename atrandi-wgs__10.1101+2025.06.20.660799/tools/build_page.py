#!/usr/bin/env python3
"""Generate atrandi_wgs.html -- the scg_lib_structs-style chemistry page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import atrandi as A
import illumina as il
import seqprimers as sp
from chemdraw import (Construct, Scene, Segment, annotation_rows, complement_segments, oligo,
                      panel, strand_row)
from page import head

OUT = HERE.parent / "atrandi_wgs.html"


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


# ============================================================== shared sequence objects
LIG_TOP = [seg("", A.LIGADAPT_STEM, "s5"), seg("", A.LIGADAPT_ARM, "s5")]
LIG_BOT = [seg("", A.LIGADAPT_BOT, "s5")]
I7_PRIMER = [seg("P7", A.P7_SEQ, "p7"), seg("i7 index", A.ATRANDI_I7_INDEX),
             seg("i7 anneal", A.ATRANDI_I7_ANNEAL, "t7")]
# NB: P5_SEQ and the Read-1 annealing portion overlap by "ACAC" -- do not concatenate them,
# split the real primer instead, so the rendered line always equals ATRANDI_P5 exactly.
P5_PRIMER = [seg("", A.P5_SEQ, "p5"), seg("", A.ATRANDI_P5[len(A.P5_SEQ):], "s5")]
# The same primer split where it matters for annealing: 25-nt tail + 20-nt arm-pairing 3' end.
P5_PRIMER_ANNEAL = [seg("P5 tail", A.ATRANDI_P5_TAIL, "p5"),
                    seg("P5 anneal", A.ATRANDI_P5_ANNEAL, "s5")]
INDEX1_SEQ = il.INDEX1_PRIMER                        # standard TruSeq Index 1 read primer
_I1_TAIL = len(INDEX1_SEQ) - len(A.ATRANDI_I7_ANNEAL)  # its 5' GATCGG
INDEX1_PRIMER = [seg("I1 tail", INDEX1_SEQ[:_I1_TAIL], "t7"),
                 seg("I1 anneal", INDEX1_SEQ[_I1_TAIL:], "t7")]
_R2_ARM = len(A.ATRANDI_I7_ANNEAL)
READ2_PRIMER = [seg("R2 5'", A.TRUSEQ_R2[:_R2_ARM], "t7"), seg("R2 3'", A.TRUSEQ_R2[_R2_ARM:], "t7")]

lib = A.final_library()
PRE = 4  # width of the "5'- " prefix


# =========================================================================== page parts
def preamble() -> str:
    return f"""<div class="wrap">
<h1>Atrandi semi-permeable capsules + PTA &mdash; single-microbe whole-genome sequencing</h1>

<p><info>Protocol: <a href="https://doi.org/10.1101/2025.06.20.660799">Gourl&eacute; et al.</a>, using
the Atrandi Single-Microbe DNA Barcoding Kit with PTA amplification and no debranching step.</info></p>

<div class="caveat">
<b>Inferred structure.</b> Atrandi does not publish the barcode-cassette strands, linkers or
Read&nbsp;2 arm. Dotted regions and steps labelled <b>INFERRED</b> show the minimal model used here.
</div>
"""


def oligos() -> str:
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
    ]
    our_rows = [oligo(f"{k} i7 + P7", [seg("", A.P7_SEQ, "p7"), seg("", v),
                                       seg("", A.ATRANDI_I7_ANNEAL, "t7")])
                for k, v in A.OUR_I7_INDEX.items()]
    our_rows.append(oligo("Generic P5", P5_PRIMER))
    return f"""<h2>Adapter and primer sequences</h2>
<seq>
{chr(10).join(rows)}
</seq>

<h3>Indexing primers</h3>
<seq>
{chr(10).join(our_rows)}
</seq>
"""


# ------------------------------------------------------------------------------ steps
def step_pta() -> str:
    # The template is given 5'->3' and drawn 3'->5', so the primer reads 5'->3' left to right
    # and phi29 extends to the right, towards the template's 5' end.
    ph = dict(placeholder=True)
    tmpl = [seg("ahead", "XXXXXXXXXXXXXXXXXXXXXXXX...XXXXXXXX", **ph), seg("site", "XXXXXX", **ph),
            seg("behind", "XXXXXXXX...XXXXXXXX", **ph)]
    prime = Scene()
    prime.strand("template", tmpl, label="genomic DNA", rev=True)
    prime.anneal("primer", [seg("primer", "N" * A.PTA_PRIMER_LEN, "tso", **ph)], to="template",
                 pair=("primer", "site"), label="primer")
    prime.arrow("primer", "phi29 extends, displacing any downstream strand,", length=12)
    prime.note("primer", "...until an alpha-thio-ddN is incorporated")
    amp = Scene()
    amp.strand("amplicon", [seg("primer", "N" * A.PTA_PRIMER_LEN, "tso", **ph),
                            seg("body", "XXXXXXXXXXXX...XXXXXXXXXXXX", **ph),
                            seg("ddN", "N", "w1", **ph)], label="")
    amp.mark("amplicon", "primer", "random primer: 5'-OH, NpNpNpNpsNpsN (two 3' phosphorothioates)")
    amp.mark("amplicon", "ddN", "alpha-thio-ddN: no 3'-OH")
    p1 = panel(prime.rows(), caption="Random primers anneal across the denatured genome; phi29 "
                                     "extends with strand displacement.")
    p2 = panel(amp.rows(), caption="One amplicon strand as synthesised (displaced, so single-"
                                   "stranded). Both termini are chemically distinctive and both "
                                   "matter downstream.")
    return f"""<h3>(2) Whole-genome amplification by PTA (BioSkryb ResolveDNA #100954; 30&nbsp;&deg;C
2.5&nbsp;h &rarr; 65&nbsp;&deg;C 5&nbsp;min)</h3>
{p1}
{p2}
"""


def step_endprep() -> str:
    body = seg("body", "XXXXXXXXXXXX...XXXXXXXXXXXX", placeholder=True)
    good = Scene()
    good.strand("top", [body, seg("dA", "A")], label="", mod5="p")
    good.anneal("bottom", [*complement_segments([body]), seg("dA", "A")], to="top",
                pair=("body'", "body"), label="", mod5="p")
    good.mark("top", "dA", "dA")
    good.mark("bottom", "dA", "dA")
    rows = good.rows()
    dead = Scene()
    dead.strand("top", [seg("5' end", "N"), body, seg("ddN", "N", "w1")], label="")
    dead.anneal("bottom", [seg("5' end", "N"), *complement_segments([body]),
                           seg("ddN", "N", "w1")], to="top", pair=("ddN", "5' end"), label="")
    dead.mark("top", "ddN", "alpha-thio-ddN: no 3'-OH, cannot be A-tailed")
    dead.mark("bottom", "ddN", "alpha-thio-ddN")
    dead = dead.rows()
    return f"""<h3>(3) End preparation (20&nbsp;&deg;C 30&nbsp;min &rarr; 65&nbsp;&deg;C
30&nbsp;min): blunting, 5'-phosphorylation and dA-tailing</h3>
{panel(rows, caption="Ends carrying a normal 3'-OH are blunted, kinased and dA-tailed, ready for "
                     "T-overhang ligation.")}
{panel(dead, caption="Ends terminating in the alpha-thio-dideoxynucleotide are not: no 3'-OH to "
                     "extend, and the phosphorothioate linkage resists the exonuclease that would "
                     "otherwise trim it back.")}
"""


def _bc(letter):
    return seg(f"BC-{letter}", letter * A.BARCODE_LEN, "cbc", placeholder=True)


def _linker(name="L"):
    return seg(name, "L" * A.LINKER_LEN, "r2", placeholder=True, inferred=True)


INSERT = seg("insert", "XXXXXXXX...XXXXXXXX", placeholder=True)


def step_barcode_a() -> str:
    cas_top = [_linker(), _bc("A"), seg("T", "T")]
    cas = Scene()
    cas.strand("top", cas_top, label="", mod5="p")
    cas.anneal("bottom", complement_segments([_bc("A")]), to="top", pair=("BC-A'", "BC-A"), label="",
               mod5="p")
    cas.mark("top", "L", "4-nt 5' overhang left for round B")
    cas.mark("top", "T", "3'-T overhang, joins the dA-tailed amplicon")
    cassette = cas.rows()
    # after ligation: the cassette's T sits opposite the amplicon's dA
    j_top = [_linker(), _bc("A"), seg("T", "T"), INSERT]
    j_bot = complement_segments([INSERT]) + [seg("dA", "A")] + complement_segments([_bc("A")])
    joined = Scene.duplex(j_top, on="BC-A", bottom=j_bot)
    joined.mark("bottom", "dA", "the amplicon's dA tail")
    joined = joined.rows()
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
    # 5'->3' bottom strand: the 5'-p overhang that pairs with the growing molecule's LLLL,
    # then the barcode complement.
    cas = Scene()
    cas.strand("top", [_linker(), _bc("B")], label="", mod5="p")
    cas.anneal("bottom", [*complement_segments([_linker("L join")]), *complement_segments([_bc("B")])], to="top",
               pair=("BC-B'", "BC-B"), label="", mod5="p")
    cas.mark("top", "L", "5' overhang for the next round")
    cas.mark("bottom", "L join'", "5' overhang: anneals to the growing molecule's LLLL")
    cassette = cas.rows()
    blk = A.final_library()
    block = Construct(blk.segments[blk.segments.index(blk.get("BC-D")):
                                   blk.segments.index(blk.get("insert")) + 1])
    done = [strand_row(block, "top"), strand_row(block, "bottom")]
    return f"""<h3>(5) Rounds B, C and D: three more 24-way splits, each adding a barcode</h3>
{panel(cassette, cls="small",
       caption="INFERRED -- rounds B, C and D join through a 4-nt cohesive overhang and leave "
               "another behind. All 24 barcodes within a round share one linker, which is why four "
               "barcodes are separated by exactly three linkers. The two overhangs of one "
               "cassette must differ in sequence (one linker per junction), or cassettes of a "
               "round would concatemerise; the placeholder letters do not show that.")}
{panel(done, cls="small", caption="After all four rounds.")}
<p><info>Ligation order is <b>A &rarr; B &rarr; C &rarr; D</b>, so D ends up outermost and is read
<b>first</b> in Read&nbsp;2.</info></p>
"""


def step_release_frag() -> str:
    near = seg("near", "XXXXXXXX...XXXXXXXX", placeholder=True)
    far = seg("far", "XXXXXX...XXXXXX", placeholder=True)
    top = [_linker(), _bc("A"), seg("T", "T"), near, far]
    bot = complement_segments([near, far]) + [seg("dA", "A")] + complement_segments([_bc("A")])
    sc = Scene.duplex(top, on="BC-A", bottom=bot)
    sc.mark("bottom", "far'", "FS cuts somewhere in here; only the barcode-proximal "
                              "fragment becomes a library molecule")
    frag = sc.rows()
    return f"""<h3>(6) Release from the capsules, then (7) enzymatic fragmentation
(NEBNext Ultra II FS, #E7805S; 10&nbsp;min at 37&nbsp;&deg;C &rarr; 300&ndash;700&nbsp;bp)</h3>
{panel(frag, cls="small",
       caption="The Release reagent dissolves the hydrogel shell; the FS enzyme mix then "
               "fragments, end-repairs and dA-tails in one tube.")}
"""


def step_adapter() -> str:
    top = [seg("stem", A.LIGADAPT_STEM, "s5"), seg("arm", A.LIGADAPT_ARM, "s5")]
    bot = [seg("stem'", A.LIGADAPT_BOT[:-1], "s5"), seg("T", A.LIGADAPT_BOT[-1:])]
    sc = Scene()
    sc.strand("top", top, label="", mod5="/5Phos/")
    sc.anneal("bottom", bot, to="top", pair=("stem'", "stem"), label="", mod5="/5AmMC6/")
    sc.mark("bottom", "T", "3'-T overhang (the ligation end)")
    sc.mark("top", "stem", "12-bp stem")
    sc.mark("top", "arm", "20-nt single-stranded arm = the P5 primer landing site")
    ad = sc.rows()
    both = A.adapter_both_ends()
    return f"""<h3>(8) Ligate the Atrandi adapter (20&nbsp;&deg;C, 15&nbsp;min &mdash; no USER step)</h3>
{panel(ad, cls="small",
       caption="The Atrandi ligation adapter. This is NOT a forked/Y adapter: the 13-nt bottom "
               "strand pairs over all 12 of its 5'-proximal bases, leaving one arm, not two.")}
{panel([strand_row(both, "top"), strand_row(both, "bottom")], cls="small",
       caption="A fragment that received the adapter but NO barcode: the same Read-1 / P5 site at "
               "both ends.")}
"""


def pre_pcr_scene() -> Scene:
    """Barcoded, adapter-ligated molecule. The bottom strand ends at the adapter's 13-nt
    amino-blocked oligo, so the 20-nt arm is single-stranded."""
    pp = A.pre_pcr()
    assert pp.segments[-1].name == "TruSeq Read 1"
    top = pp.segments[:-1] + [seg("stem", A.LIGADAPT_STEM, "s5"), seg("arm", A.LIGADAPT_ARM, "s5")]
    return Scene.duplex(top, bottom=complement_segments(top[:-1]), mod5="/5AmMC6/")


def step_pcr() -> str:
    sc = pre_pcr_scene()
    # i7 primer = top-strand sense: it anneals to the BOTTOM strand; P5 primer's 3' end is the
    # revcomp of the arm: it anneals to the TOP strand's single-stranded arm.
    sc.anneal("i7 primer", I7_PRIMER, to="bottom", pair=("i7 anneal", "TruSeq Read 2'"))
    sc.arrow("i7 primer", "")
    sc.anneal("P5 primer", P5_PRIMER_ANNEAL, to="top", pair=("P5 anneal", "arm"), above=True)
    sc.arrow("P5 primer", "")
    sc.strands["top"].label, sc.strands["bottom"].label = "top", "bottom"
    return f"""<h3>(9) Indexing PCR (98&nbsp;&deg;C 20&nbsp;s / 54&nbsp;&deg;C 30&nbsp;s /
72&nbsp;&deg;C 20&nbsp;s, 8&ndash;14 cycles)</h3>
{panel(sc.rows(),
       caption="The P5 primer lands on the top strand's 20-nt single-stranded adapter arm with "
               "zero slack -- its 25-nt P5 tail hangs off the end. The i7 primer lands on the "
               "bottom strand's copy of the Read-2 arm supplied by the round-D cassette, its P7 "
               "and index tail unpaired. Each tail is copied in on the next cycle.")}
"""


def step_final() -> str:
    rows = [strand_row(lib, "top"), strand_row(lib, "bottom")] + annotation_rows(lib, prefix_width=PRE)
    return f"""<h3>(10) Final library structure</h3>
{panel(rows, caption="INFERRED -- the barcode cassette and Read-2 arm are not published; dotted "
                     "segments show the modelled portions.")}
<p><info>Total {len(lib) - len(lib.get("insert"))}&nbsp;bp excluding the insert. There is <b>no UMI</b> in this chemistry, and
the i5 position carries no index &mdash; Atrandi's P5 primer is the plain universal primer, so i5 is
optional.</info></p>
"""


def lib_scene() -> Scene:
    sc = Scene.duplex(lib.segments)
    sc.strands["top"].label, sc.strands["bottom"].label = "top", "bottom"
    return sc


_R2_CONTENT = {"BC-D": "Barcode D (ligated last, read first)",
               "BC-A": "Barcode A (ligated first, read last)"}


def read2_layout_rows() -> str:
    """Read-2 offsets, computed from A.r2_read() (the same construct the selftest pins to
    Bascet)."""
    rows, pos = [], 0
    for s in A.r2_read():
        if s.name == "insert":
            rows.append(f"<tr><td>{pos}&ndash;</td><td>&mdash;</td><td>genomic insert</td></tr>")
            break
        what = (_R2_CONTENT.get(s.name, f"Barcode {s.name[3:]}") if s.name.startswith("BC-")
                else "linker" if s.tag == "r2" else "dA/dT ligation junction")
        if s.inferred:
            what += " <i>(inferred)</i>"
        span = f"{pos}" if len(s) == 1 else f"{pos}&ndash;{pos + len(s) - 1}"
        rows.append(f"<tr><td>{span}</td><td>{len(s)}</td><td>{what}</td></tr>")
        pos += len(s)
    return "\n".join(rows)


def sequencing() -> str:
    # Read 1 and Index 1 primers are bottom-strand sense: they anneal to the TOP strand, drawn
    # above it. The Read 2 primer is top-strand sense: it anneals to the BOTTOM strand.
    r1 = lib_scene()
    r1.anneal("Read 1 primer", [seg("R1", A.TRUSEQ_R1, "s5")], to="top", pair=("R1", "dA"),
              above=True)
    r1.arrow("Read 1 primer", "Read 1: insert")
    full_arm = A.D_ARM_INCLUDES_CCGATCT
    i1 = lib_scene()
    i1.anneal("Index 1 primer", INDEX1_PRIMER, to="top", pair=("I1 anneal", "TruSeq Read 2"),
              above=True, unpaired=[] if full_arm else ["I1 tail"])
    i1.arrow("Index 1 primer", "Index 1: i7", length=4)
    if not full_arm:
        i1.mark("Index 1 primer", "I1 tail",
                f"5' {INDEX1_SEQ[:_I1_TAIL]}: no partner on a 27-nt arm, hangs off")
    r2 = lib_scene()
    r2.anneal("Read 2 primer", READ2_PRIMER, to="bottom", pair=("R2 5'", "TruSeq Read 2'"),
              unpaired=[] if full_arm else ["R2 3'"])
    if full_arm:
        r2.arrow("Read 2 primer", "Read 2: D, C, B, A, insert")
    else:
        r2.mark("Read 2 primer", "R2 3'", f"3'-terminal {A._D_ARM_TAIL} has no partner on a 27-nt "
                                          "arm: a stock Read 2 primer cannot extend")
    return f"""<h2>Library sequencing</h2>
<p><info>Paired-end: Read&nbsp;1 is genomic insert; Read&nbsp;2 begins with the four cell barcodes.</info></p>

<h3>(1) Read 1 &mdash; genomic insert, primed from the P5 side (top strand as template)</h3>
{panel(r1.rows())}

<h3>(2) Index 1 read &mdash; the 6-bp i7 sample index (top strand as template)</h3>
{panel(i1.rows())}

<h3>(3) Read 2 &mdash; cell barcode then insert (bottom strand as template)</h3>
{panel(r2.rows())}

<div id="seqprimers"></div>
{sp.section(lib, A.SEQ_PRIMERS)}

<h2>Read 2 layout</h2>
<div class="tw">
<table>
<tr><th>Offset (0-based)</th><th>Length</th><th>Content</th></tr>
{read2_layout_rows()}
</table>
</div>
</div>
"""


def main() -> None:
    parts = [head("Atrandi PTA Library Chemistry"), preamble(), oligos(),
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
