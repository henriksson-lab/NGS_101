#!/usr/bin/env python3
"""Generate crisprscreen.html -- the pooled CRISPR screening chemistry page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import crisprscreen as cs
import illumina as il
import seqprimers as sp
from chemdraw import (Construct, Row, Scene, Segment, annotation_rows,
                      complement_segments, oligo, panel, revcomp, strand_row)
from crispr import (FILLER_LEN, OVERHANG_DOWNSTREAM, OVERHANG_UPSTREAM, SCAFFOLD_V1,
                    SPACER_LEN, U6_3PRIME, U6_PLUS1, clone_guide, guide_oligos)
from page import caveat, head, info, legend, table
from plasmid import BsmBI, Plasmid, amplify, digest, find_both, read_genbank

OUT = HERE.parent / "crisprscreen.html"
PL = HERE.parents[0] / "ref" / "plasmids"
PRE = 4
SPACER = cs.DEMO_SPACER                   # not G-initiated: every figure carries the +1 G
# vector top strand up to the upstream cut: the GPP site minus the CACC + G the insert
# supplies.
U6_TAIL = cs.SITE_U6_SHORT[:-(len(OVERHANG_UPSTREAM) + len(U6_PLUS1))]
SCAF_HEAD = SCAFFOLD_V1[:24]              # vector top strand, from the downstream cut

# the array-oligo flanks, from the Zhang lab's own design_library.py -- the same U6 3' end
# and scaffold head as everything else on this page, so they are cut from those.
ARRAY_5 = "TTTCTTGG" + cs.SITE_U6_LONG + U6_PLUS1
ARRAY_3 = SCAFFOLD_V1[:43]

# Every size and segment on this page is computed from these maps, so there is no page to
# build without them. They are third-party files the repo does not ship (see
# ref/plasmids/MANIFEST.md), so say exactly which one is missing and where it comes from
# rather than dying in read_genbank with a traceback.
MAP_FILES = (
    ("addgene-49535_lentiCRISPRv1.gb", "lentiCRISPR v1", "Addgene #49535"),
    ("addgene-52961_lentiCRISPRv2.gb", "lentiCRISPRv2", "Addgene #52961"),
    ("addgene-52963_lentiGuide-Puro.gb", "lentiGuide-Puro", "Addgene #52963"),
)


def require_maps() -> dict:
    missing = [(f, o) for f, _n, o in MAP_FILES if not (PL / f).exists()]
    if missing:
        root = HERE.parents[1]
        for f, origin in missing:
            try:
                where = (PL / f).relative_to(root)
            except ValueError:
                where = PL / f
            print(f"cannot build {OUT.name}: missing {where} -- download it from {origin} "
                  f"(https://www.addgene.org/{origin.split('#')[-1]}/); "
                  "see ref/plasmids/MANIFEST.md", file=sys.stderr)
        sys.exit(1)
    return {n: read_genbank(str(PL / f)) for f, n, _o in MAP_FILES}


VEC = require_maps()
CLONED = {n: clone_guide(p, SPACER) for n, p in VEC.items()}
# the same three, cloned with a G-initiated spacer: published sizes are quoted for one of
# those, and the whole amplicon shifts by a base between the two cases.
CLONED_G = {n: clone_guide(p, cs.DEMO_SPACER_G) for n, p in VEC.items()}


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


def ph(name, top, tag=None):
    """A schematic/stand-in region: brackets, ellipses, N-runs."""
    return Segment(name=name, top=top, tag=tag, placeholder=True)


ELLIPSIS = "..."                          # drawn on both strands; pairs with itself


# -------------------------------------------------- amplicons as segment tables
def amplicon(vector: str, fwd: str, rev: str, label: str) -> Construct:
    return amplicon_from(amplify(CLONED[vector], fwd, rev)[0].seq, label)


def amplicon_from(s: str, label: str) -> Construct:
    """Segment a finished Illumina amplicon, whichever route produced it."""
    u6 = s.find(U6_3PRIME)
    sp = u6 + len(U6_3PRIME) + 1
    sc = sp + SPACER_LEN
    r2 = s.find(revcomp(cs.READ2_IDX))
    p7 = s.find(revcomp(cs.P7))
    scaffold_end = min(sc + len(SCAFFOLD_V1), r2)
    segs = [seg("Illumina P5", s[0:len(il.P5)], "p5"),
            seg("Read 1 site", s[len(il.P5):u6 - 10], "s5"),
            seg("U6 promoter", s[u6 - 10:u6], "me"),
            seg("U6 3' end", U6_3PRIME, "r3"),
            seg("", U6_PLUS1, None),
            seg("sgRNA spacer", "N" * SPACER_LEN, "cbc", placeholder=True),
            seg("scaffold", s[sc:scaffold_end], "r2")]
    if scaffold_end < r2:
        segs.append(seg("vector", s[scaffold_end:r2], "w1"))
    segs += [seg("Read 2 site", s[r2:r2 + len(cs.READ2_IDX)], "t7"),
             seg("i7", s[r2 + len(cs.READ2_IDX):p7], None),
             seg("Illumina P7", s[p7:], "p7")]
    return Construct(segs, name=label)


GPP_F = cs.gpp_p5("")
KERMIT = cs.gpp_p7(cs.GPP_P7_INDEX_A01, cs.SITE_CPPT)
BEAKER = cs.gpp_p7(cs.GPP_P7_INDEX_A01, cs.SITE_EFS)
JF, JR = cs.joung_fwd(cs.JOUNG_STAGGERS[0]), cs.joung_rev(cs.JOUNG_REV_INDEX_1)

V2A_F, V2A_R = cs.V2ADAPTOR_F, cs.V2ADAPTOR_R
GRAFT = cs.V2ADAPTOR_GRAFT                # the reverse primer's non-templated 5' tail
ANNEAL = cs.V2ADAPTOR_ANNEAL              # its 3' half, which actually anneals


def step1_construct(vector: str) -> Construct:
    """The v2-adaptor step-1 product, segmented from the real amplicon.

    The point of the figure is the last segment: it is in the product but not in the
    template, because it rode in on the reverse primer's 5' tail.
    """
    s = amplify(CLONED[vector], V2A_F, V2A_R)[0].seq
    u6 = s.find(U6_3PRIME)
    sp = u6 + len(U6_3PRIME) + 1
    sc = sp + SPACER_LEN
    graft = len(s) - len(GRAFT)
    return Construct([
        seg("U6 promoter (Shalem F1 anneals here)", s[:u6], "s5"),
        seg("U6 3' end", U6_3PRIME, "r3"),
        seg("", U6_PLUS1, None),
        seg("sgRNA spacer", "N" * SPACER_LEN, "cbc", placeholder=True),
        seg("scaffold", s[sc:sc + len(SCAFFOLD_V1)], "r2"),
        seg("vector-specific", s[sc + len(SCAFFOLD_V1):graft], "w1"),
        seg("cPPT -- GRAFTED, absent from the template", s[graft:], "t7"),
    ], name=f"step 1 on {vector}")


def step2_construct(vector: str) -> Construct:
    """Step 2: the ordinary GPP pair, run on the step-1 product."""
    lin = Plasmid(f"step1-{vector}", amplify(CLONED[vector], V2A_F, V2A_R)[0].seq,
                  False, [])
    return amplicon_from(amplify(lin, GPP_F, KERMIT)[0].seq, f"step 2 on {vector}")


AMPS = [("GPP ARGON + KERMIT, on lentiGuide-Puro",
         amplicon("lentiGuide-Puro", GPP_F, KERMIT, "KERMIT")),
        ("GPP ARGON + BEAKER, on lentiCRISPRv2",
         amplicon("lentiCRISPRv2", GPP_F, BEAKER, "BEAKER")),
        ("Zhang / Joung one-step, on any of the three",
         amplicon("lentiCRISPRv2", JF, JR, "Joung"))]


# ------------------------------------------------------------------- sections
def preamble() -> str:
    return f"""<div class="wrap">
<h1>Pooled CRISPR screening &mdash; cloning and readout</h1>

{info("""A pooled screen has two pieces of chemistry worth drawing. At the front,
<b>how a 20-nt guide gets into a lentiviral vector</b> &mdash; one BsmBI digest and a pair
of annealed oligos. At the back, <b>how it is read out of genomic DNA</b> &mdash; a PCR
whose primers are, unexpectedly, not interchangeable between vectors.""")}

{legend(f"""<b>Everything here is computed from real Addgene maps</b> in
<code>ref/plasmids/</code>, not quoted. The BsmBI overhangs, the
{FILLER_LEN}&nbsp;bp filler, every amplicon size and every claim about which primer works
on which vector is recomputed when this page is built, and asserted in
<code>tools/selftest.py</code>.""")}

{caveat("""<b>Two corrections to the usual account</b>, both confirmed from the maps:
there is <b>no v1-vs-v2 sgRNA scaffold difference</b> &mdash; lentiCRISPR v1, lentiCRISPRv2
and lentiGuide-Puro all carry the identical original 76-nt scaffold; and
<b>lentiCas9-Blast has no BsmBI site at all</b>, being Cas9 delivery only. The real
v1&rarr;v2 difference that breaks primers is the <i>position of the cPPT/CTS element</i>.""")}
"""


def oligos_section() -> str:
    top, bot = guide_oligos(SPACER)
    rows = [
        oligo("guide oligo, top", [seg("", OVERHANG_UPSTREAM, "r3"), seg("", U6_PLUS1, None),
                                   seg("", "N" * 20, "cbc", placeholder=True)]),
        oligo("guide oligo, bottom", [seg("", revcomp(OVERHANG_DOWNSTREAM), "r2"),
                                      seg("", "n" * 20, "cbc", placeholder=True),
                                      seg("", "C", None)]),
        oligo("array oligo (pooled libraries)",
              [seg("", ARRAY_5, "r3"), seg("", "N" * 20, "cbc", placeholder=True),
               seg("", ARRAY_3, "r2")]),
        oligo("GPP P5_ARGON (stagger 0)", [seg("", cs.P5_READ1, "s5"),
                                           seg("", cs.SITE_U6_SHORT, "r3")]),
        oligo("GPP P7_KERMIT", [seg("", cs.P7, "p7"), seg("", "[8-bp i7]", None, placeholder=True),
                                seg("", cs.READ2_IDX, "t7"), seg("", cs.SITE_CPPT, "w1")]),
        oligo("GPP P7_BEAKER", [seg("", cs.P7, "p7"), seg("", "[8-bp i7]", None, placeholder=True),
                                seg("", cs.READ2_IDX, "t7"), seg("", cs.SITE_EFS, "w1")]),
        oligo("Joung NGS-Lib-Rev", [seg("", cs.P7, "p7"), seg("", "[8-bp i7]", None, placeholder=True),
                                    seg("", cs.READ2_IDX, "t7"), seg("", cs.SITE_SCAFFOLD, "r2")]),
    ]
    return f"""<h2>Oligos</h2>
<seq>
{chr(10).join(rows)}
</seq>
{info(f"""The guide oligo pair is the whole of small-scale cloning: anneal them and they
present exactly the overhangs a BsmBI-cut vector leaves. Note the three reverse readout
primers differ <b>only</b> in their 3' tail &mdash; <w1>cPPT</w1>, <w1>EFS</w1> or
<r2>scaffold</r2> &mdash; and that is what decides which vector each one works on.""")}
"""


def step_vector() -> str:
    sc = Scene()
    sc.strand("vector", [
        ph("lead", ELLIPSIS), seg("U6 tail", U6_TAIL, "r3"),
        ph("cut5", "[", "me"),
        ph("filler", f"{OVERHANG_UPSTREAM}--- {FILLER_LEN}-bp filler ---", "me"),
        ph("cut3", "]", "me"),
        seg("scaffold", SCAF_HEAD, "r2"), ph("tail", ELLIPSIS)], label="")
    sc.mark("vector", "cut5", "BsmBI cuts here")
    sc.mark("vector", "cut3", "and here")
    rows = sc.rows()
    return f"""<h2>Part 1 &mdash; getting a guide into the vector</h2>
<h3>(1) The empty vector</h3>
{panel(rows, cls="small",
       caption="Between the U6 promoter and the sgRNA scaffold sits a %d-bp filler "
               "flanked by two BsmBI sites. The filler is what makes an uncut vector "
               "harmless: it cannot express a guide." % FILLER_LEN)}
{info(f"""Computed from the map: <b>exactly two BsmBI sites</b> in each of lentiCRISPR v1,
lentiCRISPRv2 and lentiGuide-Puro, <b>{FILLER_LEN}&nbsp;bp</b> apart in all three.""")}
"""


def step_digest() -> str:
    # The upstream cut: the U6-side fragment keeps a recessed top strand, so its bottom
    # strand runs 4 nt further and presents a 5' overhang -- revcomp(CACC), 5'->3'.
    u6 = seg("U6 tail", U6_TAIL, "r3")
    sl = Scene()
    sl.strand("top", [ph("lead", ELLIPSIS), u6], label="")
    sl.anneal("bottom", [seg("overhang", revcomp(OVERHANG_UPSTREAM), "r3"),
                         *complement_segments([u6]), ph("lead'", ELLIPSIS)],
              to="top", pair=("U6 tail'", "U6 tail"), label="")
    sl.mark("bottom", "overhang", "bottom strand protrudes: pairs with the insert's "
                                  + OVERHANG_UPSTREAM)
    # The downstream cut: the scaffold-side fragment's TOP strand protrudes instead, and
    # the 4 nt it protrudes by are GTTT, the scaffold's own first four bases.
    scaf = seg("scaffold", SCAF_HEAD[len(OVERHANG_DOWNSTREAM):], "r2")
    sr = Scene()
    sr.strand("top", [seg("overhang", OVERHANG_DOWNSTREAM, "r2"), scaf, ph("tail", ELLIPSIS)],
              label="")
    sr.anneal("bottom", [ph("tail'", ELLIPSIS), *complement_segments([scaf])],
              to="top", pair=("scaffold'", "scaffold"), label="")
    sr.mark("top", "overhang", "top strand protrudes: " + OVERHANG_DOWNSTREAM
            + ", the scaffold's own first 4 nt")
    left, right = sl.rows(), sr.rows()
    return f"""<h3>(2) BsmBI digestion</h3>
{panel(left, cls="small",
       caption="The upstream cut. BsmBI is CGTCTC(1/5) -- it cuts outside its own "
               "recognition site, so the overhang it leaves is whatever sequence happens "
               "to sit there.")}
{panel(right, cls="small", caption="The downstream cut.")}
{info(f"""The two overhangs are <b>{OVERHANG_UPSTREAM}</b> and <b>{OVERHANG_DOWNSTREAM}</b>,
and neither is arbitrary: <b>{OVERHANG_UPSTREAM}</b> completes the U6 promoter's 3' end and
<b>{OVERHANG_DOWNSTREAM}</b> is the <i>scaffold's own first four bases</i>. Because the two
differ, the vector cannot re-close on itself and the insert can only go in one way round.
Both BsmBI recognition sites leave with the filler, so the product cannot be re-cut.""")}
"""


def step_ligate() -> str:
    # The annealed pair: top oligo CACC + G + spacer, bottom oligo AAAC + revcomp(spacer)
    # + C (crispr.guide_oligos). The two 5' ends overhang at opposite ends; the bottom
    # oligo's terminal C is what pairs with the +1 G.
    spacer = ph("spacer", "N" * SPACER_LEN, "cbc")
    sc = Scene()
    sc.strand("top", [seg("CACC", OVERHANG_UPSTREAM, "r3"), seg("+1 G", U6_PLUS1),
                      spacer], label="", mod5="p")
    sc.anneal("bottom", [seg("AAAC", revcomp(OVERHANG_DOWNSTREAM), "r2"),
                         *complement_segments([spacer]), seg("C", "C")],
              to="top", pair=("spacer'", "spacer"), label="", mod5="p")
    sc.mark("top", "CACC", "5' overhang " + OVERHANG_UPSTREAM)
    sc.mark("bottom", "C", "pairs with the +1 G")
    sc.mark("bottom", "AAAC", "5' overhang " + revcomp(OVERHANG_DOWNSTREAM)
            + " at the other end, for the vector's " + OVERHANG_DOWNSTREAM)
    rows = sc.rows()

    top = [ph("lead", ELLIPSIS), seg("U6", U6_3PRIME, "r3"), seg("+1 G", U6_PLUS1),
           ph("spacer", "N" * SPACER_LEN, "cbc"), seg("scaffold", SCAF_HEAD, "r2"),
           ph("tail", ELLIPSIS)]
    jn = Scene()
    jn.strand("top", top, label="")
    jn.anneal("bottom", complement_segments(top), to="top", pair=("spacer'", "spacer"),
              label="")
    jn.mark("top", "U6", "U6 promoter")
    jn.mark("top", "+1 G", "+1 G")
    jn.mark("top", "spacer", "%d-nt spacer" % SPACER_LEN)
    jn.mark("top", "scaffold", "scaffold")
    joined = jn.rows()
    return f"""<h3>(3) Anneal the guide oligos and ligate</h3>
{panel(rows, cls="small",
       caption="The annealed oligo pair presents exactly the overhangs the digest left.")}
{panel(joined, cls="small",
       caption="The product. U6, the +1 G, the spacer and the scaffold are now contiguous "
               "-- verified against the real map for both G- and non-G-initiated spacers.")}
{caveat(f"""<b>The +1 G is a standing off-by-one.</b> Pol III starts at a G, so the guide
must begin with one. A spacer that already starts with G supplies it; one that does not gets
a G <i>appended</i>. So <code>{U6_3PRIME}</code> is the promoter and the G after it is the
transcript's first base &mdash; treat the motif <code>{U6_3PRIME}G</code> as "promoter" and
add a G of your own and you have counted it twice. Every amplicon below shifts by one base
depending on the spacer's first letter.""")}
"""


def step_array() -> str:
    sc = Scene()
    sc.strand("array oligo", [seg("left arm", ARRAY_5, "r3"),
                              ph("spacer", "N" * SPACER_LEN, "cbc"),
                              seg("right arm", ARRAY_3, "r2")], label="")
    sc.mark("array oligo", "left arm",
            "%d nt: U6 3' end + the +1 G" % len(ARRAY_5), ch="-")
    sc.mark("array oligo", "spacer", "spacer", ch="-")
    sc.mark("array oligo", "right arm",
            "%d nt: scaffold 1-%d" % (len(ARRAY_3), len(ARRAY_3)), ch="-")
    rows = sc.rows()
    return f"""<h3>(4) For a pooled library, the same join is made by Gibson assembly</h3>
{panel(rows, cls="long",
       caption="The synthesised array oligo, %d nt. Flanks verified against the real "
               "vector: %d nt of U6 on the left, %d nt of scaffold on the right."
               % (len(ARRAY_5) + SPACER_LEN + len(ARRAY_3), len(ARRAY_5), len(ARRAY_3)))}
{info("""Amplified to a 141-bp insert that carries <b>60 bp of homology on each arm</b>, so
a pool of 10&#8309; guides can be assembled in one reaction rather than cloned one at a time.
The junction it produces is identical to the oligo-pair route above.""")}
"""


PAIRS = (("Shalem F1/R1", cs.SHALEM_F1, cs.SHALEM_R1),
         ("GPP ARGON + KERMIT", GPP_F, KERMIT),
         ("GPP ARGON + BEAKER", GPP_F, BEAKER),
         ("Joung one-step", JF, JR))


def pair_size(vec: str, fwd: str, rev: str) -> int | None:
    """Length of the product on the G-spacer clone, or None if the pair gives none."""
    a = amplify(CLONED_G[vec], fwd, rev)
    return a[0].length if a else None


def pair_cell(n: int | None) -> str:
    """One cell of the size table: computed, with the verdict following from the size."""
    if n is None:
        return "no product"
    if n > cs.MAX_USABLE_AMPLICON:
        return f"<b>no product</b> &mdash; {n:,} bp the long way round"
    if n > 500:
        return f"{n} bp &mdash; unusable"
    return f"<b>{n} bp</b>"


def step_readout() -> str:
    sz = {(nm, v): pair_size(v, f, r) for nm, f, r in PAIRS for v in CLONED_G}
    vecs = ("lentiCRISPR v1", "lentiCRISPRv2", "lentiGuide-Puro")
    sc = Scene()
    sc.strand("amplicon", [ph("U6", "[ U6 ]"), ph("spacer", "[spacer]", "cbc"),
                           ph("scaffold", "[-- scaffold --]", "r2"),
                           ph("vector", "[---- vector-specific ----]", "w1")], label="")
    sc.mark("amplicon", "U6", "GPP P5_ARGON and every other forward primer land HERE")
    sc.mark("amplicon", "scaffold", "Joung reverse: lands HERE, in the scaffold")
    sc.mark("amplicon", "vector", "KERMIT / BEAKER: land HERE, outside it")
    rows = sc.rows()
    return f"""<h2>Part 2 &mdash; reading the guide back out of genomic DNA</h2>
<h3>(5) Where the primers land, and why it matters</h3>
{panel(rows, cls="small",
       caption="All four published forward primers anneal in the U6 promoter. The reverse "
               "primers do not agree with each other.")}
{info("""<r2>Joung's reverse primer anneals inside the scaffold</r2>, which is identical in
every guide vector &mdash; so that pair is <b>vector-independent</b>. <w1>KERMIT and BEAKER
anneal outside it</w1>, in the cPPT and the EFS leader respectively, so each works only on
the vectors that carry its site in the right place.""")}
{caveat(f"""<b>Site presence is the wrong test; product length is the right one.</b> Computed
against the real maps: BEAKER's site <i>is</i> present in lentiGuide-Puro, but gives a
{sz[("GPP ARGON + BEAKER", "lentiGuide-Puro")]}&nbsp;bp product rather than the intended
short one. And KERMIT on lentiCRISPRv2 "works" only in the sense of producing a
<b>{sz[("GPP ARGON + KERMIT", "lentiCRISPRv2")] / 1000:.1f}&nbsp;kb wrap-around</b>, because
the cPPT there sits upstream of U6 and the two primers end up pointing the wrong way round
the plasmid.""")}
{table(["primer pair", *vecs],
       [[nm, *(pair_cell(sz[(nm, v)]) for v in vecs)] for nm, _f, _r in PAIRS])}
{info(f"""Every size in the table is computed from the Addgene map, for a G-initiated
spacer; add 1&nbsp;nt for one that is not. GPP's published figure for lentiCRISPRv2 is
285&nbsp;nt &mdash; reproduced here as
{sz[("GPP ARGON + BEAKER", "lentiCRISPRv2")]}&nbsp;nt.""")}
"""


def step_final() -> str:
    out = []
    for title, con in AMPS:
        rows = [strand_row(con, "top"), strand_row(con, "bottom")] + \
               annotation_rows(con, prefix_width=PRE)
        out.append(f"<h3>{title} &mdash; {len(con)} bp</h3>\n"
                   + panel(rows, cls="long"))
    return "<h3>(6) Final library structures</h3>\n" + "\n".join(out) + info(
        """Compare the third with the first two: because Joung's reverse primer sits in the
        <r2>scaffold</r2>, its amplicon contains no <w1>vector-specific</w1> segment at all,
        which is exactly why the same pair works on every backbone.""")


def step_twostep() -> str:
    """Part 3: why a second PCR exists, and what each round produces."""
    # --- why: where the cPPT sits relative to U6, per vector -----------------
    # Both maps in one Scene, so the two strands share a left edge and a label gutter.
    dash = lambda n="-": ph("", n)
    sc = Scene()
    V2, LG = "lentiCRISPRv2", "lentiGuide-Puro / v1"
    sc.strand(V2, [ph("cPPT", "[cPPT]", "w1"), dash("--"), ph("U6", "[ U6 ]", "s5"), dash(),
                   ph("spacer", "[spacer]", "cbc"), dash(),
                   ph("scaffold", "[scaffold]", "r2"), dash(),
                   ph("EFS", "[EFS/EF-1a]", "w1")])
    sc.mark(V2, "cPPT", "KERMIT anneals here, but points AWAY from the guide")
    sc.mark(V2, "U6", "ARGON anneals here, pointing -->")
    sc.strand(LG, [ph("U6", "[ U6 ]", "s5"), dash(), ph("spacer", "[spacer]", "cbc"), dash(),
                   ph("scaffold", "[scaffold]", "r2"), dash(), ph("cPPT", "[cPPT]", "w1")])
    sc.mark(LG, "U6", "ARGON anneals here, pointing -->")
    sc.mark(LG, "cPPT", "KERMIT anneals here, pointing <-- back at the guide")
    sc.blank(before=LG)                                        # blank line between the two
    down = sc.rows()

    # --- the reverse primer's two halves -------------------------------------
    pr = Scene()
    pr.strand("v2-adaptor reverse", [seg("graft", GRAFT, "w1"), seg("anneal", ANNEAL, "s5")],
              label="")
    pr.mark("v2-adaptor reverse", "graft",
            "%d nt of cPPT: NOT in the template" % len(GRAFT), ch="-")
    pr.mark("v2-adaptor reverse", "anneal",
            "%d nt that anneal in the vector" % len(ANNEAL), ch="-")
    prim = pr.rows()

    # every number quoted in the prose below, computed once here
    TWO = ("lentiCRISPRv2", "lentiGuide-Puro")
    r2_size = {v: len(step2_construct(v)) for v in TWO}
    vspec = {v: len([x for x in step1_construct(v).segments
                     if x.name == "vector-specific"][0].top) for v in TWO}
    secs = []
    for vec in TWO:
        for n, con in ((1, step1_construct(vec)), (2, step2_construct(vec))):
            rows = [strand_row(con, "top"), strand_row(con, "bottom")] + \
                   annotation_rows(con, prefix_width=PRE)
            secs.append(f"<h3>{vec} &mdash; round {n}, {len(con)} bp</h3>\n"
                        + panel(rows, cls="long"))

    # --- sizes, cloned vs leftover filler, computed both ways ----------------
    rows_tbl, r1_sizes = [], []
    for vec in TWO:
        a_cl = amplify(CLONED[vec], V2A_F, V2A_R)
        a_st = amplify(VEC[vec], V2A_F, V2A_R)
        a_sh = amplify(CLONED[vec], cs.SHALEM_F1, cs.SHALEM_R1)
        r1_sizes += [a[0].length for a in (a_cl, a_st) if a]
        rows_tbl.append([
            vec,
            f"<b>{a_cl[0].length} bp</b>" if a_cl else "no product",
            f"<b>{a_st[0].length:,} bp</b>" if a_st else "no product",
            f"{r2_size[vec]} bp",
            (f"{a_sh[0].length} bp" if a_sh and a_sh[0].length < cs.MAX_USABLE_AMPLICON
             else "<b>no usable product</b>")])

    return f"""<h2>Part 3 &mdash; the two-step readout</h2>
<h3>(7) Why a second PCR exists at all</h3>
{panel(down, cls="small",
       caption="The same site, in two different places. A one-step readout needs a reverse "
               "primer site downstream of the guide; on lentiCRISPRv2 the cPPT is upstream "
               "of U6 instead.")}
{info("""On <b>lentiGuide-Puro</b> and <b>lentiCRISPR v1</b> the cPPT sits
<w1>downstream</w1> of the scaffold, so <b>ARGON + KERMIT in one step</b> is all that is
needed. On <b>lentiCRISPRv2</b> the cPPT sits <w1>upstream of U6</w1>: both primers still
find their sites, but they point apart, and the only product runs the long way round the
plasmid. That is the gap the two-step readout fills.""")}

<h3>(8) Round 1 grafts the missing site on</h3>
{panel(prim, cls="small",
       caption="The v2-adaptor reverse primer. Only its 3' half anneals; the 5' half is "
               "carried into the product as new sequence.")}
{info(f"""This is the whole trick: <b>round 1 does not amplify a site, it
<i>installs</i> one</b>. The 5' {len(GRAFT)}&nbsp;nt are non-templated, so every round-1
product ends in a cPPT sequence the plasmid never had &mdash; and round 2 can then use the
ordinary <w1>KERMIT</w1> primer that lentiCRISPRv2 would otherwise refuse.""")}
{info("""Its forward primer is <b>identical to Shalem's</b> <code>F1</code>. The two
systems differ in the reverse primer <i>only</i>, which makes them easy to confuse at the
bench &mdash; and the failure is silent, because the wrong one still gives a band.""")}

<h3>(9) What each round actually produces</h3>
{info(f"""Round 1 carries <b>no Illumina sequence at all</b> &mdash; it is a pure vector
PCR whose only job is to install the <t7>cPPT</t7>. Round 2 adds <p5>P5</p5>,
<p7>P7</p7>, the read sites and the index. Note where the two vectors differ: in round 1
the <w1>vector-specific</w1> segment is {vspec["lentiCRISPRv2"]}&nbsp;nt on lentiCRISPRv2
and {vspec["lentiGuide-Puro"]}&nbsp;nt on lentiGuide-Puro, and by round 2 it is gone from
both.""")}
{chr(10).join(secs)}

<h3>(10) Sizing it, and what the gel tells you</h3>
{table(["vector", "round 1, cloned", "round 1, filler still present",
        "round 2 (either)", "Shalem one-step, cloned"], rows_tbl)}
{caveat(f"""<b>Round 2 is {r2_size["lentiCRISPRv2"]}&nbsp;bp on both vectors, so it cannot
tell them apart.</b> Round 1's reverse primer grafts the round-2 reverse site on, so the
round-2 amplicon runs U6 &rarr; grafted tail and contains none of the stretch that differs
between the backbones. <b>Round 1 is the only informative size</b> &mdash; and it separates
a cloned library from leftover filler by
~{(len(VEC["lentiCRISPRv2"]) - len(CLONED["lentiCRISPRv2"])) / 1000:.1f}&nbsp;kb, because
the filler lies inside the amplicon. Size it against a <b>100&nbsp;bp ladder</b>: the
diagnostic range is {min(r1_sizes)}&nbsp;bp to {max(r1_sizes) / 1000:.1f}&nbsp;kb and all
the information is at the bottom of the gel.""")}
{info("""Sizes are computed from the real Addgene maps with a G-initiated spacer; add
1&nbsp;nt for a spacer that is not. The filler column uses the empty vector directly &mdash;
both primer sites lie outside the excised region, so the amplicon shrinks by exactly as
much as the plasmid does.""")}
"""


def step_stagger() -> str:
    sc = Scene()
    longest = max(cs.GPP_STAGGERS, key=len)
    for st in cs.GPP_STAGGERS:
        name = f"stagger {len(st)}"
        sc.strand(name, [seg("stagger", st, "me"), seg("U6", cs.SITE_U6_SHORT, "r3"),
                         ph("spacer", "N" * SPACER_LEN, "cbc")], label=name)
        if st == "":
            sc.note(name, "read 1 cycle 1 starts here, whatever the stagger")
            sc.mark(name, "spacer", "spacer, no stagger")
        elif st == longest:
            sc.mark(name, "spacer", "spacer, stagger %d: %d cycles later"
                    % (len(st), len(st)))
    rows = sc.rows()
    return f"""<h2>Sequencing</h2>
<h3>(11) The stagger ladder</h3>
{panel(rows, cls="small",
       caption="Read 1 begins immediately after the sequencing primer. Without a stagger "
               "every cluster on the flow cell would be reading the identical U6 sequence "
               "at the same cycle -- a monotemplate the instrument cannot base-call. Eight "
               "primers of different length push the spacer to eight different cycles.")}
{info("""GPP's ladder is <b>0, 1, 2, 3, 4, 6, 7, 8</b> &mdash; there is no 5. Zhang/Joung
use ten staggers of 9&ndash;18&nbsp;nt instead. The ladders are <b>not interchangeable</b>:
a counting script tuned for one expects the guide at the wrong cycle in the other. Zhang's
<code>count_spacers.py</code> searches read positions 30&ndash;55, correct for its own
ladder and wrong for GPP's.""")}
{info(f"""Sequencing is an 8-cycle index 1 plus a read 1 long enough to clear the longest
stagger &mdash; Joung specify <b>80 cycles</b>, and &gt;500 reads per sgRNA for a screen.
<b>No custom sequencing primer is needed</b>: every readout primer on this page ends in a
stock TruSeq handle, which is what the next section checks, one library at a time.""")}

{chr(10).join(sp.section(con, cs.SEQ_PRIMERS,
                         heading=f"({12 + i}) Sequencing primers — {title}",
                         intro="Read 1 starts in the U6 promoter, so its first bases are "
                               "the stagger and the promoter's 3' end, not the spacer."
                               if "Joung" not in title else
                               "This pair's read 2 starts inside the <r2>scaffold</r2>, "
                               "which is why it is backbone-independent.")
              for i, (title, con) in enumerate(AMPS))}
</div>
"""


def main() -> None:
    parts = [head("Pooled CRISPR Screening Chemistry"), preamble(), oligos_section(),
             step_vector(), step_digest(), step_ligate(), step_array(),
             step_readout(), step_final(), step_twostep(), step_stagger()]
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
