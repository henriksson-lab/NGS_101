#!/usr/bin/env python3
"""Generate crisprmip.html -- the CRISPR-MIP chemistry page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import crisprmip as cm
import illumina as il
import seqprimers as sp
from chemdraw import (Construct, Row, Scene, Segment, annotation_rows, bridge,
                      complement_segments, oligo, panel, revcomp, strand_row)
from crispr import SCAFFOLD_V1, clone_guide
from padlock import capture
from page import caveat, head, info, legend, table
from plasmid import amplify, read_genbank

OUT = HERE.parent / "crisprmip.html"
PRE = 4
SPACER = cm.EXAMPLE_SPACER
VEC = HERE.parents[1] / "lenticrispr-v1-screening__10.1126+science.1247005" / "ref" / "plasmids" / "addgene-52961_lentiCRISPRv2.gb"

# The page is computed from the real vector map, so it cannot be built without it. That map
# is third-party material and is not committed; fail with the one thing the reader needs.
if not VEC.exists():
    print(f"cannot build {OUT.name}: missing {VEC.relative_to(HERE.parents[1])}"
          " -- the lentiCRISPRv2 map (Addgene #52961). It is third-party material and is not"
          " committed; see lenticrispr-v1-screening__10.1126+science.1247005/ref/plasmids/"
          "MANIFEST.md for how to re-obtain it.", file=sys.stderr)
    sys.exit(1)

probe = cm.probe()
vector = clone_guide(read_genbank(str(VEC)), SPACER)
cap = capture(vector, probe)[0]
lib_seq = amplify(cap.circle, cm.P5_TRACR_FWD, cm.P7_TRACR_REV, min_anneal=15)[0]
FILL = cap.fill


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


# ------------------------------------------------- drawing helpers (all computed)
# Every multi-strand drawing below is a chemdraw.Scene: strands are given 5'->3' and placed
# by which segment pairs with which, and the Scene refuses any column that does not pair.
# The helpers here only add what Scene lacks -- a padlock's loop -- and take every column
# from Scene/Construct spans, never from a typed number.

def _cells_to_chunks(cells: dict[int, tuple[str, str | None]]) -> list[tuple]:
    out, col = [], 0
    for c in sorted(cells):
        ch, tag = cells[c]
        if c > col:
            out.append((" " * (c - col), None, False))
        if out and out[-1][1] == tag and c == col:
            out[-1] = (out[-1][0] + ch, tag, False)
        else:
            out.append((ch, tag, False))
        col = c + 1
    return out


def padlock_rows(sc: Scene, arms: list[str], inner, label: str) -> list[Row]:
    """A Scene whose `arms` (annealed above their template, left to right) are the two ends
    of ONE oligo, drawn on one row with the backbone arching between their outer ends.

    The Scene has already checked every arm column against the template. Here the arm
    strands' separate rows are merged into one: outer end labels are dropped (the backbone
    continues there, into the bridge's risers) and only the ends facing the gap keep
    theirs (3'-OH, 5'-p). With one arm strand (the closed circle) no end is free at all.
    """
    origin, gutter = sc.origin()
    sts = [sc.strands[a] for a in arms]
    cells: dict[int, tuple[str, str | None]] = {}

    def put(col, text, tag=None):
        for i, ch in enumerate(text):
            if col + i in cells:
                raise ValueError(f"padlock arms overlap at column {col + i}")
            cells[col + i] = (ch, tag)

    for st in sts:
        pos = origin + st.col
        for s in st.drawn():
            put(pos, s.top, s.tag)
            pos += len(s)
    if len(sts) > 1:
        put(origin + sts[0].end(), sts[0].ends()[1])            # left arm: its 3' end faces the gap
        put(origin + sts[-1].col - len(sts[-1].ends()[0]), sts[-1].ends()[0])
    probe_row = Row(chunks=[(sts[0].label.ljust(gutter), None, False)]
                    + _cells_to_chunks({c - gutter: v for c, v in cells.items()}))
    left = origin + sts[0].col
    right = origin + sts[-1].end() - 1
    loop = bridge(left, right, inner, label=label, riser=1)
    rest = sc.rows(omit=arms)
    return loop + [probe_row] + rest


def band(width: int, text: str) -> str:
    """'<---- text ---->' exactly `width` columns wide."""
    lab = f" {text} "
    if len(lab) + 2 > width:
        raise ValueError(f"band label {text!r} does not fit in {width} columns")
    pad = width - 2 - len(lab)
    return "<" + "-" * (pad // 2) + lab + "-" * (pad - pad // 2) + ">"


def span_bands(con: Construct, items, lead: int) -> Row:
    """One row of bands under a Construct row, each spanning segments first..last.
    `lead` is the width of whatever precedes the construct's first base on its row."""
    cells: dict[int, tuple[str, str | None]] = {}
    for first, last, text in items:
        s, e = con.span(first)[0], con.span(last)[1]
        for i, ch in enumerate(band(e - s, text)):
            if lead + s + i in cells:
                raise ValueError(f"bands overlap at column {lead + s + i}")
            cells[lead + s + i] = (ch, None)
    return Row(chunks=_cells_to_chunks(cells))


BACKBONE_INNER = [[("[13-nt UMI]", "umi", False), ("  ", None, False),
                   ("[Read 2 site]", "t7", False), ("  [i7]  ", None, False),
                   ("[Read 1 site]", "s5", False)]]


# ------------------------------------------------- the final library, as segments
def library() -> Construct:
    s = lib_seq.seq
    return Construct([
        seg("Illumina P5", s[0:len(il.P5)], "p5"),
        seg("vector 3'", s[len(il.P5):48], "r2"),     # captured: past the scaffold
        seg("lig arm", cm.LIG_ARM, "r1"),
        seg("UMI", "N" * 13, "umi", placeholder=True),
        seg("Read 2 site", cm.READ2_SITE, "t7"),
        seg("i7", cm.PROBE_INDICES[0], None),
        seg("Read 1 site", cm.READ1_SITE, "s5"),
        seg("ext arm", cm.EXT_ARM, "r3"),
        seg("", "G", None),                           # the Pol III +1
        seg("sgRNA spacer", "N" * 20, "cbc", placeholder=True),
        seg("scaffold", s[195:245], "r2"),
        seg("Illumina P7", s[245:], "p7"),
    ], name="CRISPR-MIP library")


LIB = library()


def preamble() -> str:
    return f"""<div class="wrap">
<h1>CRISPR-MIP &mdash; counting sgRNAs with a padlock probe</h1>

{info("""A pooled CRISPR screen is normally read out by PCR off genomic DNA, and counted as
reads. That conflates two things: how many cells carried a guide, and how well that guide's
amplicon happened to amplify. <b>CRISPR-MIP replaces the readout PCR with a padlock
capture carrying a UMI</b>, so every captured genomic molecule is labelled <i>before</i>
any amplification and the readout becomes a count of molecules rather than of reads.""")}

{info("""Selinger M, Yakovenko I, Nazir I, Henriksson J. <i>CRISPR-MIP replaces PCR and
reveals GC and oversampling bias in pooled CRISPR screens.</i> bioRxiv 2024.03.28.587082.
<a href="https://doi.org/10.1101/2024.03.28.587082">doi:10.1101/2024.03.28.587082</a>.
Sequences are from Table S2; every length, position and product below is recomputed from
the real lentiCRISPRv2 map (Addgene #52961) each time this page is built.""")}

{legend("""<b>The idea in one line.</b> A padlock probe finds its target with <i>two</i>
arms instead of one primer, the gap between them is filled and sealed into a circle, and an
exonuclease then destroys everything that is still linear. The failure modes are removed by
enzymology rather than by a cleanup &mdash; and the probe's UMI ends up covalently attached
to the molecule it captured.""")}

{caveat(f"""<b>The probe is {cm.PROBE_LEN}&nbsp;nt, not the {cm.PROBE_LEN_IN_TEXT} stated
in the paper's Methods text.</b> Table S2's sequences are {cm.PROBE_LEN}&nbsp;nt and
decompose exactly as drawn below ({len(cm.LIG_ARM)} + {cm.UMI_LEN} + {len(cm.READ2_SITE)}
+ {cm.INDEX_LEN} + {len(cm.READ1_SITE)} + {len(cm.EXT_ARM)} = {cm.PROBE_LEN}). This page
follows Table S2.""")}
"""


def oligos() -> str:
    rows = [
        oligo("CRISPR-MIP probe (MIP_probe_1 of 9)",
              [seg("", cm.LIG_ARM, "r1"), seg("", "[13-bp UMI]", "umi", placeholder=True),
               seg("", cm.READ2_SITE, "t7"), seg("", f"[{cm.INDEX_LEN}-bp i7]", None, placeholder=True),
               seg("", cm.READ1_SITE, "s5"), seg("", cm.EXT_ARM, "r3")], mods="/5Phos/"),
        oligo("P5_tracrRNA_fwd", [seg("", il.P5, "p5"), seg("", cm.P5_TRACR_TAIL, "r2")]),
        oligo("P7_tracrRNA_rev", [seg("", il.P7, "p7"), seg("", cm.P7_TRACR_TAIL, "r2")]),
    ]
    idx = table(["probe", "i7 index"],
                [[f"MIP_probe_{i+1}", f"<code>{x}</code>"]
                 for i, x in enumerate(cm.PROBE_INDICES)])
    return f"""<h2>Oligos</h2>
<seq>
{chr(10).join(rows)}
</seq>
{info(f"""The probe is one ssDNA molecule. Its <r1>ligation arm</r1> and
<r3>extension arm</r3> are the only parts complementary to the target; everything between
them is backbone. <b>P5 and P7 appear nowhere in the probe</b> &mdash; they arrive on two
separate primers whose 3' tails anneal inside the <i>captured</i> sequence, which is what
makes an un-extended probe unamplifiable.""")}
{info(f"""The <t7>Read 2</t7> and <s5>Read 1</s5> sites are TruSeq sequence (where each
sequencing primer lands is computed under <i>Sequencing primers</i> below). The nine probes
differ <b>only</b> in an {cm.INDEX_LEN}-nt i7 index:""")}
{idx}
"""


# ------------------------------------------------- the target and the capture
# Top-strand-sense pieces of the captured region. The 112-nt gap is drawn with its two ends
# real and its middle elided, so every drawn base can be checked against its partner.
GAP_HEAD = 10                                   # scaffold bases shown after the spacer
GAP_TAIL = 8                                    # bases shown before the ligation arm


def gap_segments(suffix: str = "") -> list[Segment]:
    """The gap-fill, top-strand sense: +1 G, spacer, scaffold head ... last bases."""
    assert FILL[0] == "G" and len(FILL) == cap.gap
    return [seg("+1" + suffix, FILL[0]),
            seg("spacer" + suffix, "N" * 20, "cbc", placeholder=True),
            seg("scaffold" + suffix, FILL[21:21 + GAP_HEAD], "r2"),
            seg("elided" + suffix, "...", placeholder=True),
            seg("gap end" + suffix, FILL[-GAP_TAIL:], "r2")]


def target_segments() -> list[Segment]:
    """Vector top strand across the capture: ext-arm site, gap, lig-arm site."""
    return [seg("ext site", cm.EXT_ARM, "r3"), *gap_segments(" site"),
            seg("lig site", cm.LIG_ARM, "r1")]


def scene_target() -> Scene:
    top = [seg("U6 3' end", cm.EXT_ARM, "r3"), seg("+1", "G"),
           seg("spacer", "N" * 20, "cbc", placeholder=True),
           seg("scaffold", SCAFFOLD_V1[:24], "r2")]
    sc = Scene()
    sc.strand("top", top, label="")
    sc.anneal("bottom", complement_segments(top), to="top", pair=("spacer'", "spacer"),
              label="")
    sc.mark("bottom", "U6 3' end'", "U6 3' end")
    sc.mark("bottom", "+1'", "+1 (Pol III start)")
    sc.mark("bottom", "spacer'", "sgRNA spacer")
    sc.mark("bottom", "scaffold'", "scaffold")
    return sc


def scene_capture(filled: bool) -> tuple[Scene, list[str]]:
    """The probe's two arms on the template strand -- the strand OPPOSITE the one the arm
    sequences are written against. Unfilled: two arm strands with a gap between them.
    Filled: the extension arm, the new strand and the ligation arm are one continuous
    strand, so it is annealed as one."""
    sc = Scene()
    sc.strand("template", complement_segments(target_segments()), rev=True)
    if filled:
        sc.anneal("probe", [seg("ext arm", cm.EXT_ARM, "r3"), *gap_segments(),
                            seg("lig arm", cm.LIG_ARM, "r1")],
                  to="template", pair=("ext arm", "ext site'"), above=True)
        sc.mark("template", "+1 site'", f"copied: {cap.gap} nt filled in above it",
                through="gap end site'")
        return sc, ["probe"]
    sc.anneal("ext arm", [seg("ext arm", cm.EXT_ARM, "r3")], to="template",
              pair=("ext arm", "ext site'"), above=True, label="probe", mod3="OH")
    sc.anneal("lig arm", [seg("lig arm", cm.LIG_ARM, "r1")], to="template",
              pair=("lig arm", "lig site'"), above=True, label="probe", mod5="p")
    sc.mark("template", "+1 site'", f"gap, {cap.gap} nt", through="gap end site'")
    return sc, ["ext arm", "lig arm"]


def step_target() -> str:
    return f"""<h2>Step-by-step</h2>
<h3>(1) The target: one integrated provirus per cell</h3>
{panel(scene_target().rows(), cls="small",
       caption="The guide cassette of lentiCRISPRv2, carrying whichever 20-nt spacer that "
               "cell received. This is what the probe has to find, in a genome.")}
"""


def step_hybridise() -> str:
    sc, arms = scene_capture(filled=False)
    rows = padlock_rows(sc, arms, BACKBONE_INNER, "probe backbone")
    return f"""<h3>(2) Denature and hybridise &mdash; 94&nbsp;&deg;C, ramp to 60&nbsp;&deg;C at
&minus;0.1&nbsp;&deg;C/s, then overnight at 60&nbsp;&deg;C</h3>
{panel(rows, cls="small",
       caption="This is the shape the molecule is named for. Both arms anneal to the same "
               "strand and the backbone arches over the gap between them, leaving the "
               "probe draped across the target like a padlock waiting to be closed. "
               "The middle of the gap is elided (...).")}
{info(f"""Two arms, not one primer. A PCR primer needs one match; a padlock needs
<b>both</b> arms, in order, on the same strand, a fixed distance apart &mdash; which in
the whole 13&nbsp;kb vector happens at exactly one place.""")}
{info(f"""The arm T<sub>m</sub>s are deliberately unequal: <r1>ligation arm</r1>
{probe.arm_tms()[1]:.1f}&nbsp;&deg;C vs <r3>extension arm</r3>
{probe.arm_tms()[0]:.1f}&nbsp;&deg;C. The ligation arm must be the <b>hotter</b> of the
two, because the polymerase travels toward it &mdash; if it were the weaker duplex the
polymerase would displace it and the circle would never close.""")}
{info("""Note which strand is which: the arm sequences are written in the sense of the
vector's top strand, so the probe anneals to the strand <i>opposite</i> them &mdash; the
template shown beneath.""")}
"""


def flat_circle() -> Construct:
    """The closed circle, cut open at the extension arm's 5' end and laid flat."""
    return Construct([seg("ext arm", cm.EXT_ARM, "r3"), *gap_segments(),
                      seg("lig arm", cm.LIG_ARM, "r1"),
                      seg("UMI", "N" * cm.UMI_LEN, "umi", placeholder=True),
                      seg("Read 2", "[Read2]", "t7", placeholder=True),
                      seg("i7", "[i7]", None, placeholder=True),
                      seg("Read 1", "[Read1]", "s5", placeholder=True)], name="circle")


def step_fill() -> str:
    sc, arms = scene_capture(filled=True)
    filled = padlock_rows(sc, arms, BACKBONE_INNER, "probe backbone")
    con = flat_circle()
    lead = "  .-> "
    circ = strand_row(con, "top", prefix=lead, suffix=" -.")
    a, b = lead.index("."), len(circ.plain()) - 1
    close = Row(chunks=[(" " * a + "'" + "-" * (b - a - 1) + "'", None, False)])
    bands = span_bands(con, [("ext arm", "ext arm", "probe"),
                             ("+1", "gap end", "captured target"),
                             ("lig arm", "Read 1", "probe backbone")], len(lead))
    return f"""<h3>(3) Extend and ligate &mdash; Ampligase + Phusion HF, 60&nbsp;&deg;C, 1&nbsp;h</h3>
{panel(filled, cls="small",
       caption="The polymerase fills the gap from the extension arm's 3'-OH; the ligase "
               "then joins that new 3' end to the probe's 5' phosphate. Arm, new strand "
               "and arm are now one continuous strand, and the loop is a covalently "
               "closed circle threading the target.")}
{panel([circ, close, Row(), bands], cls="small",
       caption="The same %d-nt circle, cut open and laid flat so the order of its parts is "
               "visible (gap middle and backbone sites abbreviated). The captured target "
               "is now flanked by probe sequence -- which is what 'molecular inversion' "
               "means." % len(cap.circle))}
{info(f"""The gap is <b>{cap.gap}&nbsp;nt</b>: the Pol III <b>+1 G</b>, the complete
<cbc>20-nt spacer</cbc>, the <r2>76-nt scaffold</r2> and 15&nbsp;nt beyond it. A spacer
that already begins with G gives 111 instead, because the +1 is then supplied by the spacer
rather than appended &mdash; Brunello spacers are not G-initiated, so 112 is the working
number.""")}
"""


EXO_FATES = (
    ("closed circle (probe + captured target)", "", "SURVIVES"),
    ("unreacted probe", "5'-p ----------------- 3'", "destroyed by Exo I"),
    ("annealed but never ligated", "5'-p ------/ /-------- 3'", "destroyed by Exo I"),
    ("genomic DNA", "5'- ------------------ 3'", "destroyed by Exo III"),
)


def step_exo() -> str:
    w1 = max(len(n) for n, g, _ in EXO_FATES if g) + 2
    w2 = max(len(n) if not g else w1 + len(g) for n, g, _ in EXO_FATES) + 3
    rows = [Row(chunks=[((n.ljust(w1) + g if g else n).ljust(w2) + v, None, False)])
            for n, g, v in EXO_FATES]
    return f"""<h3>(4) Exonuclease I and III &mdash; 37&nbsp;&deg;C, 45&nbsp;min</h3>
{panel(rows, cls="small",
       caption="Everything with a free end is degraded. Only correctly closed circles "
               "remain.")}
{caveat("""<b>This is the selection, and there are four of them in series</b> &mdash; none
of which is a size selection or a bead cleanup:
<br>1. both arms must anneal to the same molecule, a fixed distance apart;
<br>2. the polymerase must cross the gap <i>and</i> the ligase must close the circle;
<br>3. the exonucleases destroy everything still linear;
<br>4. the PCR primer sites lie inside the <i>captured</i> region, so an un-extended probe
has nowhere for a primer to land.""")}
"""


# ------------------------------------------------- inverse PCR
P5_ADAPTER, P5_TAIL = il.P5, cm.P5_TRACR_TAIL
P7_ADAPTER, P7_TAIL = il.P7, cm.P7_TRACR_TAIL


def pcr_circle() -> list[Segment]:
    """The whole circle, once, 5'->3' (it is single-stranded: probe + fill), opened just
    after the P7 site so that the PCR product reads left to right without wrapping."""
    p7 = FILL.find(revcomp(P7_TAIL))
    p5 = FILL.find(P5_TAIL)
    assert 0 < p7 < p5 and FILL.count(P5_TAIL) == 1
    p7e, p5e = p7 + len(P7_TAIL), p5 + len(P5_TAIL)
    return [seg("between sites", FILL[p7e:p5], "r2"),
            seg("P5 site", FILL[p5:p5e], "r2"),
            seg("fill end", FILL[p5e:], "r2"),
            seg("lig arm", cm.LIG_ARM, "r1"),
            seg("UMI", "N" * cm.UMI_LEN, "umi", placeholder=True),
            seg("Read 2 site", cm.READ2_SITE, "t7"),
            seg("i7", cm.PROBE_INDICES[0]),
            seg("Read 1 site", cm.READ1_SITE, "s5"),
            seg("ext arm", cm.EXT_ARM, "r3"),
            seg("+1", FILL[0]),
            seg("spacer", "N" * 20, "cbc", placeholder=True),
            seg("scaffold", FILL[21:p7], "r2"),
            seg("P7 site", FILL[p7:p7e], "r2")]


def first_strand() -> list[Segment]:
    """What the P7 primer makes on the circle, 5'->3': its adapter, then the copy."""
    return [seg("P7 adapter", P7_ADAPTER, "p7"), *complement_segments(pcr_circle())]


def scene_pcr_cycle1() -> Scene:
    sc = Scene()
    sc.strand("circle", pcr_circle())
    sc.anneal("P7 primer", [seg("P7 adapter", P7_ADAPTER, "p7"), seg("P7 tail", P7_TAIL, "r2")],
              to="circle", pair=("P7 tail", "P7 site"), unpaired=["P7 adapter"])
    sc.arrow("P7 primer", "copies the circle, through backbone and arms", length=12)
    return sc


def scene_pcr_cycle2() -> Scene:
    sc = Scene()
    sc.strand("first strand", first_strand(), rev=True)
    sc.anneal("P5 primer", [seg("P5 adapter", P5_ADAPTER, "p5"), seg("P5 tail", P5_TAIL, "r2")],
              to="first strand", pair=("P5 tail", "P5 site'"), unpaired=["P5 adapter"])
    sc.arrow("P5 primer", "to the P7 adapter: the library", length=12)
    return sc


def step_pcr() -> str:
    return f"""<h3>(5) PCR off the circle</h3>
{panel(scene_pcr_cycle1().rows(), cls="long",
       caption="Inverse PCR, cycle 1. The circle is single-stranded and carries the "
               "P7 primer's site in its own sense, so only the P7 primer can anneal to it. "
               "The circle is drawn opened just past that site: on the circle the two "
               "primer sites are %d nt apart and face away from each other, so the product "
               "runs the long way round, through the probe backbone."
               % len(pcr_circle()[0]))}
{panel(scene_pcr_cycle2().rows(), cls="long",
       caption="Cycle 2. The P5 primer anneals to that first-strand copy, at a site that "
               "also lies inside the captured sequence, and extends to the P7 adapter's "
               "end. From here on it is ordinary PCR.")}
{info("""Both 3' tails lie in sequence that exists only if the gap was actually filled.
A probe that annealed but was never extended carries no primer site at all &mdash; the
fourth and last selection step.""")}
"""


def final() -> str:
    rows = [strand_row(LIB, "top"), strand_row(LIB, "bottom")] + annotation_rows(LIB, prefix_width=PRE)
    return f"""<h3>(6) Final library</h3>
{panel(rows, cls="long")}
{info(f"""<b>{len(LIB)}&nbsp;bp.</b> Note the order: the molecule begins with
<p5>P5</p5> and a piece of <r2>captured vector</r2> from just past the scaffold, then crosses the
<r1>ligation junction</r1> into the probe backbone &mdash; <umi>UMI</umi>,
<t7>Read&nbsp;2 site</t7>, i7, <s5>Read&nbsp;1 site</s5> &mdash; then crosses the
<r3>extension arm</r3> back into captured target, where the <cbc>sgRNA spacer</cbc> sits.
The inversion is why the UMI and the spacer end up in the same short amplicon.""")}
"""


def library_scene() -> Scene:
    sc = Scene()
    segs = list(LIB)
    sc.strand("top", segs, label="")
    sc.anneal("bottom", complement_segments(segs), to="top", pair=("UMI'", "UMI"), label="")
    return sc


def primer(role: str) -> sp.SeqPrimer:
    """The one declared primer for `role` that has a site on the library."""
    hits = [p for p in cm.SEQ_PRIMERS if p.role == role and sp.locate(LIB, p)]
    assert len(hits) == 1, role
    return hits[0]


def scene_primer(role: str, label: str, tag: str, text: str, length: int) -> Scene:
    """The declared `role` primer, placed where seqprimers.locate puts it on LIB. A primer
    in the top strand's sense anneals to the bottom strand; one in the opposite sense
    anneals to the top strand (the strand left after the paired-end turnaround). Any
    5'-terminal bases without a partner are drawn as unpaired."""
    p, hit = primer(role), sp.locate(LIB, primer(role))
    site = hit.covers[-1] if hit.strand == "bottom" else hit.covers[0]
    s0 = LIB.span(site)[0]
    if hit.strand == "bottom":
        sc = library_scene()
        segs = [seg(label, p.seq[hit.free5:], tag)]
        if hit.free5:
            segs.insert(0, seg("unpaired 5'", p.seq[:hit.free5], tag))
        sc.anneal(label, segs, to="bottom", pair=(label, site + "'"), shift=hit.start - s0,
                  unpaired=["unpaired 5'"])
    else:
        assert hit.free5 == 0, "draw the unpaired 5' end of a top-strand primer"
        sc = Scene()
        sc.strand("top", list(LIB), label="")
        sc.anneal(label, complement_segments([seg(label, revcomp(p.seq), tag)], suffix=""),
                  to="top", pair=(label, site), shift=hit.start - s0)
    sc.arrow(label, text, length=length)
    return sc


def scene_read1() -> Scene:
    return scene_primer("Read 1", "Read 1 primer", "s5",
                        f"Read 1, {cm.READ1_CYCLES} cycles: ext arm, +1 G, spacer", 14)


def scene_index1() -> Scene:
    return scene_primer("Index 1 (i7)", "Index 1 primer", "t7",
                        f"Index 1, {cm.INDEX_LEN} cycles: the i7", 6)


def scene_read2() -> Scene:
    return scene_primer("Read 2", "Read 2 primer", "t7",
                        f"Read 2, {cm.READ2_CYCLES} cycles: the UMI", 14)


def sequencing() -> str:
    return f"""<h2>Sequencing</h2>
{info("""Read 1, Index 1 and Read 2 prime on standard TruSeq sites carried in the probe
backbone, so no custom sequencing primer is needed. There is no i5: P5 runs straight into
captured vector sequence.""")}

<h3>Read 1 &mdash; {cm.READ1_CYCLES} cycles, the sgRNA</h3>
{panel(scene_read1().rows(), cls="long")}

<h3>Index 1 and Read 2 &mdash; {cm.INDEX_LEN} and {cm.READ2_CYCLES} cycles, the i7 and the UMI</h3>
{panel(scene_index1().rows(), cls="long",
       caption="Index 1 primes on the Read 2 site in its own sense and reads the i7 "
               "immediately 3' of it.")}
{panel(scene_read2().rows(), cls="long",
       caption="Read 2 primes on the same site from the other strand, reading the "
               "opposite way: straight into the UMI.")}

{sp.section(LIB, cm.SEQ_PRIMERS)}
{info("""<b>Why the UMI is the point.</b> It was attached during the capture, before any
amplification, so two reads sharing a UMI came from the same original genomic molecule.
Deduplicating on it converts read counts into molecule counts &mdash; which removes the
sequencing-depth bias, and is how the paper is able to show that PCR-based readouts carry a
GC% bias of their own.""")}
</div>
"""


def main() -> None:
    parts = [head("CRISPR-MIP Chemistry"), preamble(), oligos(), step_target(),
             step_hybridise(), step_fill(), step_exo(), step_pcr(), final(), sequencing()]
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
