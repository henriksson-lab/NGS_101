#!/usr/bin/env python3
"""Generate crisprumi.html -- the CRISPR-UMI (Schmierer) chemistry page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))
sys.path.insert(0, str(HERE.parents[1] / "crispr-mip__10.1101+2024.03.28.587082" / "tools"))

import crisprmip as cm
import crisprumi as cu
import seqprimers as sp
from chemdraw import (Construct, Scene, Segment, annotation_rows, complement_segments, panel,
                      strand_row, tm)
from crispr import SCAFFOLD_V1, SCAFFOLD_V1_HEAD, clone_guide
from page import caveat, head, info, legend, table
from plasmid import amplify

OUT = HERE.parent / "crisprumi.html"
PRE = 4

AU_FLIP_HEAD = cu.SCAFFOLD_AU_FLIP[:len(SCAFFOLD_V1_HEAD)]


def vector_sizes():
    """Both vectors cloned with the same stand-in guide, plus their PCR1 sizes.

    Computed on demand rather than at import, because this reads the Addgene parent map,
    which is third-party and not committed (see ref/MANIFEST.md). The selftest imports this
    module to check the drawings; those need no map, and must still run on a fresh clone.
    """
    cl_par = clone_guide(cu.parent_vector(), cu.DEMO_SPACER)
    cl_mod = cu.rebuild_vector(cloned=True)
    return (cl_par, cl_mod,
            amplify(cl_par, cu.OUR_PCR1_FW, cu.OUR_PCR1_REV, 18)[0].length,
            amplify(cl_mod, cu.OUR_PCR1_FW, cu.OUR_PCR1_REV, 18)[0].length)


def seg(n, t, tag=None, **kw):
    return Segment(name=n, top=t, tag=tag, **kw)


def scaffold_segs(s: str, diff) -> list:
    """`s` cut into segments so that every substituted base is a segment of its own.

    That is what lets a Scene mark the substitutions by NAME: nothing in this page
    computes, or types, a column number.
    """
    out, prev = [], 0
    for k, i in enumerate(diff):
        if i > prev:
            out.append(seg(f"same{k}", s[prev:i], "r2"))
        out.append(seg(f"sub{k}", s[i], "me"))
        prev = i + 1
    if prev < len(s):
        out.append(seg("rest", s[prev:], "r2"))
    return out


def cassette() -> Construct:
    return Construct([
        seg("U6 3' end", cu.U6_OVERLAP[-16:], "r3"),
        seg("sgRNA spacer", "N" * cu.SPACER_LEN, "cbc", placeholder=True),
        seg("AU-flip scaffold", cu.SCAFFOLD_AU_FLIP, "r2"),
        seg("term", cu.TERMINATOR, "me"),
        seg("Illumina Read 2 / Index 1", cu.ILLUMINA_ADAPTER, "t7"),
        seg("RSL", "N" * cu.RSL_LEN, "umi", placeholder=True),
        seg("vector", cu.DOWNSTREAM, "w1"),
    ], name="integrated cassette")


CAS = cassette()


def preamble() -> str:
    return f"""<div class="wrap">
<h1>CRISPR-UMI (Schmierer) &mdash; a UMI cloned into the guide library</h1>

{info("""Each virion carries one guide <i>and</i> one random 6-bp label, so every transduced
cell founds a uniquely tagged lineage. Counting distinct labels counts <b>clones</b>, which
is a different and better-behaved quantity than counting reads.""")}

{info("""Schmierer B, Botla SK, Zhang J, Turunen M, Kivioja T, Taipale J. <i>CRISPR/Cas9
screening using unique molecular identifiers.</i> Mol Syst Biol 2017;13(10):945.
<a href="https://doi.org/10.15252/msb.20177834">doi:10.15252/msb.20177834</a>. The authors
call the label a <b>Random Sequence Label (RSL)</b>, which is the less ambiguous name.""")}

{caveat("""<b>Two different methods are called "CRISPR-UMI".</b> This one clones the UMI
<b>into the library plasmid</b>, so it labels a transduced cell lineage and it
<i>changes the vector</i>. <a href="../crispr-mip__10.1101+2024.03.28.587082/crisprmip.html">CRISPR-MIP</a> carries its
UMI <b>on a capture probe</b> applied to genomic DNA afterwards, so it labels a captured
molecule and leaves the plasmid alone. Both improve a screen's statistics; they intervene at
opposite ends of the experiment.""")}

{legend(f"""<b>Why this page exists.</b> The vector is lentiGuide-Puro with three edits, and
each one quietly breaks something designed against the parent. Everything below is recomputed
from the real lentiGuide-Puro map each time this page is built.""")}
"""


def what_it_buys() -> str:
    return f"""<h2>What the label buys</h2>
{table(["", "analysis", "what it does"],
       [[f"<b>{a}</b>", b, c] for a, b, c in cu.ANALYSES])}
{info("""The motivation for <b>LDA</b> in one sentence: a read count cannot distinguish
"many surviving clones" from "one clone that expanded", and those mean opposite things in a
dropout screen. Counting distinct RSLs separates them &mdash; which is why precision improves
without needing more cells.""")}
"""


AU_FLIP_DIFF = [i for i, (x, y) in enumerate(zip(SCAFFOLD_V1, cu.SCAFFOLD_AU_FLIP)) if x != y]


def auflip_scene() -> Scene:
    """Original and AU-flip scaffolds, one above the other, with the two substitutions
    marked by segment name. The two are the same length and the same sense, so this is an
    alignment rather than a duplex -- but the columns still have to line up, and a Scene
    is what makes that true by construction."""
    o, a = SCAFFOLD_V1, cu.SCAFFOLD_AU_FLIP
    sc = Scene()
    sc.strand("original", scaffold_segs(o, AU_FLIP_DIFF))
    sc.strand("AU-flip", scaffold_segs(a, AU_FLIP_DIFF))
    for k, i in enumerate(AU_FLIP_DIFF):
        sc.mark("AU-flip", f"sub{k}", f"{o[i]}->{a[i]} at scaffold position {i + 1}")
    return sc


def vector_section() -> str:
    diff = AU_FLIP_DIFF
    sc = auflip_scene()
    return f"""<h2>The vector: {cu.VECTOR_NAME}</h2>
{info("""Built from <b>lentiGuide-Puro (Addgene #52963)</b> by, in the authors' words,
"replacing the sequence <code>gttttagagctagaaatagcaagttaaaa……TTTTTT</code> with
<code>gtttAagagctagaaatagcaagttTaaa……TTTTTTcgtctct</code> to create an AU-flip
(Chen <i>et al</i>, 2013) and an additional BsmBI site downstream of the tracrRNA".""")}
<h3>The AU-flip is exactly two substitutions</h3>
{panel(sc.rows(), cls="long",
       caption="Computed by diffing the two scaffolds: positions "
               + " and ".join(str(i + 1) for i in diff) + ", nothing else. The two strands "
               "are drawn by a Scene, so the carets sit under the substituted bases by "
               "construction, not by counting.")}
{info("""It removes a <code>TTTT</code> run that Pol III can read as a premature terminator
and restores pairing in the lower stem &mdash; the same motivation as the F+E scaffold, by a
much smaller edit.""")}
{caveat(f"""<b>And it silently breaks anything aimed at the scaffold's 5' half.</b>
<code>{SCAFFOLD_V1_HEAD}</code> &mdash; the diagnostic head of the original scaffold, and
the annealing site of several published primers and probes &mdash; is <b>absent from this
vector</b>. Primers that bind the scaffold's <i>3'</i> end
(<code>&hellip;{SCAFFOLD_V1[-16:]}</code>) are unaffected, which is why the Joung readout primer
survives the change and others do not.""")}
"""


def read1_scene() -> Scene:
    """The custom Read-1 primer on the finished library, drawn by pairing.

    Everything here is templated: PCR2's forward primer grafts the Read-1 site on, so by
    the time the library is sequenced the primer's 5' end has a partner too. The Scene
    checks every drawn column, which is what verifies the claim below -- the primer's 3'
    base is the Pol III +1 G, so cycle 1 reads spacer base 1 with no leader to skip.
    """
    top = [seg("Read 1 site", cu.READ1_SITE[-6:], "t7"),
           seg("U6 3' end", cu.U6_OVERLAP[-23:], "r3"),
           seg("sgRNA spacer", "N" * cu.SPACER_LEN, "cbc", placeholder=True),
           seg("scaffold", cu.SCAFFOLD_AU_FLIP[:12], "r2")]
    sc = Scene()
    sc.strand("library", complement_segments(top), label="library (template)", rev=True)
    sc.anneal("CRIPSRSEQ", [seg("Read 1 site", cu.CUSTOM_SEQ_PRIMER[:6], "t7"),
                            seg("U6 anneal", cu.CUSTOM_SEQ_PRIMER[6:], "r3")],
              to="library", pair=("U6 anneal", "U6 3' end'"))
    sc.arrow("CRIPSRSEQ", f"Read 1, {cu.READ1_CYCLES} cycles")
    sc.mark("library", "sgRNA spacer'", "read from cycle 1: no leader to skip")
    return sc


def insert_section() -> str:
    rows = [strand_row(CAS, "top"), strand_row(CAS, "bottom")] + \
        annotation_rows(CAS, prefix_width=PRE)
    return f"""<h2>The library insert</h2>
{panel(rows, cls="long",
       caption="The integrated cassette, as published. Note what sits after the terminator.")}
{caveat(f"""<b>The Illumina Read-2 / Index-1 adapter is built into the construct</b>, and the
<umi>6-bp RSL</umi> sits immediately 3' of it &mdash; that is, <b>exactly where an i7 index
would be</b>, so it is read by the Index 1 read without needing its own sequencing primer.
<br><br>This is the structural reason the readout differs: on lentiCRISPRv2 and plain
lentiGuide-Puro the Read-2 site has to be <i>grafted on</i> by a PCR primer's non-templated
tail. Here it is templated. A readout designed to add that site is solving a problem this
vector does not have.""")}
<h3>The custom Read-1 primer stops on the +1 G</h3>
{panel(read1_scene().rows(), cls="small",
       caption="The paper's own read primer (its table spells it CRIPSRSEQ) annealed to the "
               "finished library. Placed by pairing, so every drawn column is checked.")}
{info(f"""The run is {cu.READ1_CYCLES} cycles of Read 1 and then <b>two index reads</b> of
{cu.I5_CYCLES} cycles each &mdash; and the i7 is the <umi>RSL</umi>, not a sample index.
Which primer lands where, and what each read reports first, is computed from the finished
library in <a href="#seqprimers">Sequencing primers</a> below.""")}
{info(f"""Assembly: the array oligo is annealed to a single 119-bp oligo carrying the RSL and
the Illumina site, double-stranded with outer primers, and cloned by <b>Gibson assembly</b>.
Library: <b>{cu.LIBRARY['guides']:,} guides</b> over {cu.LIBRARY['genes']:,} genes, with
{cu.LIBRARY['non_targeting']} non-targeting controls.""")}
"""


def consequences() -> str:
    CL_PAR, CL_MOD, A_PAR, A_MOD = vector_sizes()
    return f"""<h2>What changes, relative to the parent vector</h2>
{table(["", "lentiGuide-Puro", cu.VECTOR_NAME],
       [["cloned vector", f"{len(CL_PAR):,} bp",
         f"<b>{len(CL_MOD):,} bp</b> (+{len(CL_MOD) - len(CL_PAR)})"],
        ["PCR1 with CRISPR_PCR1-F/R", f"{A_PAR} bp", f"<b>{A_MOD} bp</b>"],
        ["original scaffold head present", "yes", "<b>no</b>"],
        ["Illumina Read-2 adapter in vector", "no", "<b>yes</b>"],
        ["CRISPR-MIP ligation arm present", "no", "no"]])}
{info(f"""A <b>{A_MOD - A_PAR} bp size difference in PCR1 alone distinguishes the two
vectors</b>, which makes it a cheap diagnostic if you are unsure which one you have. Sanger
across the scaffold settles it outright: <code>{SCAFFOLD_V1_HEAD}</code> for the parent,
<code>{AU_FLIP_HEAD}</code> for this one.""")}
{info("""A padlock probe <i>can</i> be designed for this vector &mdash; see below.""")}
"""


LIB = cu.final_library()


def sequencing() -> str:
    rows = [strand_row(LIB, "top")] + annotation_rows(LIB, prefix_width=PRE)
    return f"""<h2>The finished library</h2>
{info(f"""Three nested PCRs off genomic DNA (14, 19 and 14 cycles) give the paper's
gel-purified <b>{cu.LIBRARY_LEN_PUBLISHED}&nbsp;bp</b> product. PCR3's forward primer adds
P5 and the i5 index; PCR2/PCR3's reverse primers add P7. The drawing is the same
construct the table below is computed from, and the selftest checks it base-for-base
against the product amplified off the rebuilt vector.""")}
{panel(rows, cls="long",
       caption=f"The PCR3 product, {len(LIB)} bp, top strand. The "
               "Read-2 / Index-1 site is the vector's own, not a primer tail.")}
<a id="seqprimers"></a>
{sp.section(LIB, cu.SEQ_PRIMERS,
            intro="Read 1 needs the paper's own spiked-in primer, the i7 read is the "
                  "<umi>RSL</umi>, and all three stock TruSeq primers that are not the "
                  "Index-1 one fail here &mdash; each for a different reason.")}
"""


def padlock_design() -> str:
    return f"""<h2>A padlock probe that works on this vector</h2>
{caveat(f"""<b>The constraint that actually shapes the design is not the AU-flip.</b> The
CRISPR-MIP probe backbone carries <code>{cm.READ2_SITE}</code>; this vector carries
<code>{cu.ILLUMINA_ADAPTER}</code> &mdash; <b>a substring of it</b>. A capture spanning the
vector's copy therefore puts two Read-2 primer sites into one amplicon, which no readout can
disambiguate.""")}
{caveat(f"""<b>But first: the RSL is {cu.RSL_DISTANCE}&nbsp;nt away, and that cannot be
shortened.</b> Between the guide and the RSL sit {len(cu.SCAFFOLD_AU_FLIP) + len(cu.TERMINATOR)
+ len(cu.ILLUMINA_ADAPTER)}&nbsp;nt of fixed construct &mdash; the {len(cu.SCAFFOLD_AU_FLIP)}-nt
scaffold, the {len(cu.TERMINATOR)}-nt terminator and the {len(cu.ILLUMINA_ADAPTER)}-nt
built-in adapter &mdash; and the extension
arm must be upstream of the guide for the guide to be captured at all. So
<b>{cu.RSL_DISTANCE}&nbsp;nt is the floor</b> for any single padlock that captures both.
The published CRISPR-MIP probe fills <b>{cu.VALIDATED_GAP}&nbsp;nt</b>, so this is
<b>{cu.RSL_DISTANCE / cu.VALIDATED_GAP:.2f}&times;</b> the only gap length the chemistry is
known to work at. That is inside the range MIP designs normally quote, but there is
<b>no efficiency data for this probe at that length</b>, and gap-fill efficiency degrades as
the gap grows. Treat option&nbsp;A as a candidate to test, not a drop-in.""")}
<h3>Option A &mdash; guide <i>and</i> RSL (gap {cu.PADLOCK_GAP} nt, untested)</h3>
{table(["", "sequence", "note"],
       [["extension arm", f"<code>{cm.EXT_ARM}</code>", "unchanged; U6 is untouched"],
        ["ligation arm", f"<code>{cu.PADLOCK_LIG_ARM}</code>",
         f"{len(cu.PADLOCK_LIG_ARM)} nt, Tm {tm(cu.PADLOCK_LIG_ARM):.1f} &deg;C, "
         f"{tm(cu.PADLOCK_LIG_ARM) - tm(cm.EXT_ARM):.1f} &deg;C above the extension arm"]])}
{info(f"""Verified: one capture site, gap <b>{cu.PADLOCK_GAP}&nbsp;nt</b>, circle
<b>{cu.PADLOCK_CIRCLE}&nbsp;nt</b>, and the captured region holds the +1&nbsp;G, the
<cbc>{cu.SPACER_LEN}-nt spacer</cbc>, the full scaffold, the terminator, the vector's adapter
<b>and the <umi>RSL</umi></b>. The existing <code>P7_tracrRNA_rev</code> site still lies
inside the capture, so that primer is unchanged &mdash; but the probe backbone must
<b>drop its own Read-2 site</b> and let the vector's serve. The RSL then sits in the index
position and is read as the i7, exactly as Schmierer designed. One amplicon then yields
molecule counts (the probe's UMI) and lineage labels (the RSL).""")}
<h3>Option B &mdash; guide only (gap {cu.PADLOCK_GAP_UNIVERSAL} nt, shorter than the
validated {cu.VALIDATED_GAP})</h3>
{info(f"""<code>{cu.PADLOCK_LIG_ARM_UNIVERSAL}</code> sits in the scaffold's 3' half at
positions {SCAFFOLD_V1.find(cu.PADLOCK_LIG_ARM_UNIVERSAL) + 1}&ndash;{SCAFFOLD_V1.find(cu.PADLOCK_LIG_ARM_UNIVERSAL) + len(cu.PADLOCK_LIG_ARM_UNIVERSAL)}, clear of both AU-flip positions. Verified to give exactly one capture,
gap {cu.PADLOCK_GAP_UNIVERSAL}, circle 192&nbsp;nt, on <b>lentiCRISPR v1, lentiCRISPRv2,
lentiGuide-Puro and this vector alike</b>. It does not capture the RSL, and its
{cu.PADLOCK_GAP_UNIVERSAL}-nt capture is too short to hold both existing PCR primer sites,
so that pair needs redesigning. If option&nbsp;A's gap proves too long this is the fallback
&mdash; but the guide and the RSL are then no longer in one molecule, and pairing them from
two independent assays is not possible.""")}
</div>
"""


def main() -> None:
    # No silent skip here: half this page is numbers measured off the real parent map, and a
    # page with those holes punched in it would be worse than no page. Refuse, in one line.
    try:
        parts = [head("CRISPR-UMI Chemistry"), preamble(), what_it_buys(), vector_section(),
                 insert_section(), consequences(), sequencing(), padlock_design()]
    except cu.MissingMap as e:
        sys.exit(f"cannot build {OUT.name}: {e}")
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
