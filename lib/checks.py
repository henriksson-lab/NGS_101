"""
Shared self-test harness, plus the checks that hold for every protocol.

`run_common()` is retained for older protocol suites. New protocol suites keep only
source-boundary behavior that cannot be enforced by the model; shared primitives are
checked once by the repository-level suite.

The point of these is regression, not discovery: every identity below was verified once
by hand, and is encoded here so that an edit to a shared sequence cannot silently break
a protocol that depends on it.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import illumina as il
import nextera as nx
from plasmid import Plasmid, amplify, bind
import rt
from chemdraw import (Construct, Scene, Segment, bases_pair, bridge, complement,
                      complement_segments, revcomp, tm)


class Check:
    def __init__(self) -> None:
        self.fails: list[str] = []
        self.skips: list[str] = []
        self.n = 0

    def __call__(self, label: str, got, want=True) -> None:
        self.n += 1
        ok = (got == want)
        if not ok:
            self.fails.append(f"{label}\n      got:  {got!r}\n      want: {want!r}")
        print(f"  {'PASS' if ok else 'FAIL'}  {label}")

    def section(self, title: str) -> None:
        print(f"\n{title}\n" + "-" * len(title))

    def skip(self, label: str, reason: str) -> None:
        """Record a check that could not run. Not a failure: a missing third-party source
        file means the clone is incomplete, not that the chemistry is wrong."""
        self.skips.append(f"{label}\n      {reason}")
        print(f"  SKIP  {label}")

    def raises(self, label: str, fn) -> None:
        try:
            fn()
        except (ValueError, KeyError):
            self(label, True)
        else:
            self(label, False)

    def report(self, exit_on_fail: bool = True) -> int:
        print("\n" + "=" * 62)
        if self.fails:
            print(f"{len(self.fails)} of {self.n} checks FAILED:\n")
            for f in self.fails:
                print("  * " + f)
            if exit_on_fail:
                sys.exit(1)
            return 1
        if self.skips:
            print(f"{len(self.skips)} check group(s) SKIPPED -- third-party source material "
                  "is not in the repository:\n")
            for sk in self.skips:
                print("  - " + sk)
            print()
        tail = f" ({len(self.skips)} skipped)" if self.skips else ""
        print(f"all {self.n} checks passed{tail}")
        return 0


@dataclass(frozen=True)
class Source:
    """A third-party file the repo does not ship: a paper's supplementary table, a vendor
    manual, a plasmid map. `origin` must say exactly where to get it again."""
    path: Path
    origin: str

    def __str__(self) -> str:
        return str(self.path)


def have(check: "Check", *sources: Source, label: str = "") -> bool:
    """True if every `source` is on disk. Otherwise record one skip naming each missing
    file and where it comes from, and return False so the caller can skip its checks.

    Third-party material is gitignored (it is not ours to redistribute), so a fresh clone
    legitimately lacks it. The checks that read it then do not run -- and say so loudly --
    rather than failing as if the chemistry were wrong.
    """
    missing = [s for s in sources if not s.path.exists()]
    if not missing:
        return True
    for s in missing:
        try:
            where = s.path.relative_to(Path(__file__).resolve().parents[1])
        except ValueError:
            where = s.path
        check.skip(label or f"checks needing {s.path.name}",
                   f"missing {where} -- obtain it from: {s.origin}")
    return False


def _raises(fn) -> bool:
    try:
        fn()
    except ValueError:
        return True
    return False


def _first_strand_overhang() -> bool:
    """RT's untemplated CCC must render left of the mRNA's 5' end, not under it."""
    m = [Segment("body", "X" * 6, placeholder=True), Segment("polyA", "AAAA")]
    sc = Scene()
    sc.strand("m", m)
    sc.anneal("c", [Segment("dT", "TTTT"), *complement_segments(m[:1]), Segment("CCC", "CCC")],
              to="m", pair=("dT", "polyA"))
    return sc.strands["c"].span("CCC")[1] == sc.strands["m"].col


def run_common(check: Check) -> None:
    """Checks that must hold regardless of which protocol is being documented."""

    check.section("chemdraw primitives")
    check("complement preserves order", complement("AATGATACGG"), "TTACTATGCC")
    check("revcomp round-trips", revcomp(revcomp("GATCGGAAGAGC")), "GATCGGAAGAGC")
    check("complement is not revcomp", complement("GATC") != revcomp("GATC"))
    check("Tm is sane for a 20-mer", 45 < tm("ACACTCTTTCCCTACACGAC") < 65)

    check.section("chemdraw guard rails")
    check.raises("non-DNA without placeholder=True is rejected",
                 lambda: Segment(name="bad", top="AAZZQQ"))
    check.raises("mismatched bottom-strand length is rejected",
                 lambda: Segment(name="bad", top="ACGT", bottom="AC"))
    check.raises("duplicate segment names are rejected",
                 lambda: Construct([Segment("x", "ACGT"), Segment("x", "ACGT")]))
    check.raises("unknown segment name is rejected",
                 lambda: Construct([Segment("x", "ACGT")]).span("nope"))
    check("placeholders are not complemented as if they were bases",
          Segment("bc", "AAAAAAAA", placeholder=True).bottom_text() != "TTTTTTTT")

    check.section("Scene: strands placed by pairing, never by indent (lib/chemdraw.py)")
    check("A pairs T, G pairs C", bases_pair("A", "T") and bases_pair("G", "C"))
    check("T does not pair X (the hand-indent error class)", bases_pair("T", "X"), False)
    check("placeholder X pairs its stand-in x", bases_pair("X", "x"))
    check("IUPAC: W pairs A, N pairs anything", bases_pair("W", "A") and bases_pair("N", "G"))
    _m = [Segment("body", "X" * 10, placeholder=True), Segment("polyA", "A" * 12)]
    _sc = Scene()
    _sc.strand("m", _m)
    _sc.anneal("dT", [Segment("dT", "T" * 8)], to="m", pair=("dT", "polyA"), shift=2)
    check("anneal computes the column itself: dT starts 2 into the polyA",
          _sc.strands["dT"].col, 12)
    check("a shift that puts T over X is refused, not drawn",
          _raises(lambda: Scene.anneal(_sc, "bad", [Segment("dT", "T" * 8)], to="m",
                                       pair=("dT", "polyA"), shift=-3)))
    check("an untemplated overhang hangs OFF the end rather than over the body",
          _first_strand_overhang())
    check("complement_segments reverse-complements real bases",
          complement_segments([Segment("a", "AACG")])[0].top, "CGTT")
    check("...and lowercases placeholders instead of complementing them",
          complement_segments([Segment("b", "XXN", placeholder=True)])[0].top, "nxx")
    check("rendered rows carry no hand-typed indent", all(r.indent == 0 for r in _sc.rows()))
    _d = Scene.duplex([Segment("a", "ACGTAC"), Segment("b", "X" * 4, placeholder=True)])
    check("Scene.duplex pairs the full complement (drawn 3'->5')",
          _d.strands["bottom"].text(), "TGCATGxxxx")
    _d.footer("end", "top", "b")
    check("footer goes last, at the segment's column",
          _d.rows()[-1].plain().index("end") - _d.rows()[0].plain().index("ACGT"), 6)
    _d.stack("bottom", "top")
    check("stack reorders whole strands", _d.rows()[0].plain().lstrip().startswith("3'"))
    _d.blank(before="top")
    check("blank inserts an empty row above a strand", _d.rows()[1].plain(), "")
    check("rows(omit=) drops a strand's row", len(_d.rows(omit=["top"])), len(_d.rows()) - 1)
    check("complement_segments honours an explicit bottom (T:A junction)",
          complement_segments([Segment("j", "T", placeholder=True, bottom="A")])[0].top, "A")

    check.section("sequencing primers: one source, located by computation (lib/seqprimers.py)")
    import seqprimers as sp
    check("TruSeq Index 2 (RC workflow) is revcomp of TruSeq Read 1",
          il.INDEX2_PRIMER_RC, revcomp(il.TRUSEQ_READ1))
    check("the standard sets reference lib constants, not copies",
          (sp.TRUSEQ["R1"].seq is il.TRUSEQ_READ1, sp.NEXTERA["R1"].seq is nx.READ1_PRIMER),
          (True, True))
    _lib = Construct([Segment("P5", il.P5), Segment("R1", il.TRUSEQ_READ1[4:]),
                      Segment("insert", "X" * 20, placeholder=True)])
    check("an exact primer is located, and reads the base after its 3' end",
          sp.locate(_lib, sp.TRUSEQ["R1"]).reads_from, "insert")
    _short = Construct([Segment("R1", il.TRUSEQ_READ1[:-1]), Segment("x", "AGAGACAG")])
    check("a 3'-terminal mismatch is NOT located (extension needs the 3' end)",
          sp.locate(_short, sp.TRUSEQ["R1"]), None)
    check("...but a 5' flap is, if the paired 3' footprint is >= MIN_ANNEAL, and is reported",
          sp.locate(Construct([Segment("R1", il.TRUSEQ_READ1[6:]), Segment("x", "GGGG")]),
                    sp.TRUSEQ["R1"]).free5, 6)
    check("...and a footprint shorter than MIN_ANNEAL is not a site",
          sp.locate(Construct([Segment("R1", il.TRUSEQ_READ1[-10:]), Segment("x", "GGGG")]),
                    sp.TRUSEQ["R1"]), None)
    check("verify() demands all four roles", any("Read 2" in e for e in
          sp.verify(_lib, [sp.TRUSEQ["R1"]])))

    check.section("bridge / loop rendering (lib/chemdraw.py)")
    rows = bridge(4, 40, [[("hello", None, False)]], label="backbone")
    check("a bridge is arc + content + riser", len(rows), 3)
    check("its risers land on the columns asked for",
          rows[0].plain().index(".") == 4 and rows[0].plain().rindex(".") == 40)
    check("every row is the same width", len({len(r.plain()) for r in rows}), 1)
    check("the label is carried in the arc", "backbone" in rows[0].plain())
    check("content sits inside the bars", rows[1].plain().strip().startswith("|")
          and rows[1].plain().rstrip().endswith("|"))
    check.raises("a bridge narrower than 4 columns is rejected",
                 lambda: bridge(0, 2, []))
    check.raises("content too wide for the bridge is rejected",
                 lambda: bridge(0, 8, [[("far too much text", None, False)]]))

    check.section("PCR prediction (lib/plasmid.py)")
    T = "AAGCTTGGATCCGAATTCACGTTCAGGTACCCTGCAGTCTAGAGGCCTA"   # no internal repeats
    lin = Plasmid("lin", T, False, [])
    fwd, rev = T[:16], revcomp(T[-16:])
    check("each test primer binds exactly once",
          len(bind(lin, fwd, 16)) == 1 and len(bind(lin, rev, 16)) == 1)
    a = amplify(lin, fwd, rev, min_anneal=16)
    check("an amplicon running to the last base of a LINEAR template is not lost", bool(a))
    check("...and spans the whole template", a[0].length if a else 0, len(T))
    check("a linear template never reports a wrap", not a[0].wraps_origin if a else False)
    check("the same template as a circle still amplifies",
          bool(amplify(Plasmid("circ", T, True, []), fwd, rev, min_anneal=16)))
    check("primers facing away on a linear template give nothing",
          amplify(lin, T[-16:], revcomp(T[:16]), min_anneal=16), [])
    check("...but on a circle they amplify the long way round",
          bool(amplify(Plasmid("circ", T, True, []), T[-16:], revcomp(T[:16]),
                       min_anneal=16)))

    check.section("canonical Illumina sequences")
    check("P5 is 29 nt", len(il.P5), 29)
    check("P7 is 24 nt", len(il.P7), 24)
    check("TruSeq Read 1 primer is 33 nt", len(il.TRUSEQ_READ1), 33)
    check("TruSeq Read 2 primer is 34 nt", len(il.TRUSEQ_READ2), 34)
    check("the 12-bp stem self-pairs", revcomp(il.STEM_COMPLEMENT), il.STEM)
    check("Read-2 trim string == revcomp(Read 1 primer)",
          revcomp(il.TRUSEQ_READ1), il.TRIM_SEEN_IN_READ2)
    check("Read-1 trim string == revcomp(Read 2 primer), less the final base",
          revcomp(il.TRUSEQ_READ2)[:-1], il.TRIM_SEEN_IN_READ1)
    check("both trim strings start with the dA tail",
          il.TRIM_SEEN_IN_READ1[0] == "A" and il.TRIM_SEEN_IN_READ2[0] == "A")

    check.section("NEBNext hairpin adaptor")
    check("hairpin is 65 nt", len(il.NEBNEXT_HAIRPIN), 65)
    check("the dU sits at position 32 (1-based), exactly mid-loop",
          il.NEBNEXT_HAIRPIN[il.NEBNEXT_HAIRPIN_DU_POS], "U")
    check("USER splits the hairpin into exactly the two published arms",
          il.NEBNEXT_ARM_READ2 + "U" + il.NEBNEXT_ARM_READ1, il.NEBNEXT_HAIRPIN)
    check("the Read-1 arm IS the full TruSeq Read 1 primer",
          il.NEBNEXT_ARM_READ1, il.TRUSEQ_READ1)
    check("the Read-2 arm is truncated by exactly 2 nt",
          len(il.CANONICAL_ARM_READ2) - len(il.NEBNEXT_ARM_READ2),
          il.NEBNEXT_READ2_TRUNCATION)
    check("the canonical Read-2 arm == revcomp(Read 2 primer) less the dA",
          revcomp(il.TRUSEQ_READ2)[1:], il.CANONICAL_ARM_READ2)
    check("the hairpin's two halves pair over the 12-bp stem",
          revcomp(il.NEBNEXT_HAIRPIN[:12]), il.NEBNEXT_HAIRPIN[-13:-1])

    check.section("NEBNext index primers")
    check("4 index sets present",
          all(len(d) > 0 for d in (il.NEBNEXT_I5_SET1, il.NEBNEXT_I7_SET1,
                                   il.NEBNEXT_I5_SET2, il.NEBNEXT_I7_SET2)))
    for name, d in (("i5 set 1", il.NEBNEXT_I5_SET1), ("i7 set 1", il.NEBNEXT_I7_SET1),
                    ("i5 set 2", il.NEBNEXT_I5_SET2), ("i7 set 2", il.NEBNEXT_I7_SET2)):
        check(f"{name}: all indices are 8 nt of real DNA",
              all(len(v) == 8 and set(v) <= set("ACGT") for v in d.values()))
        check(f"{name}: no duplicate index sequences", len(set(d.values())), len(d))
    check("i7 primer assembles to P7 + index + Read 2 primer",
          il.nebnext_i7_primer("AGGAGGAA"), il.P7 + "AGGAGGAA" + il.TRUSEQ_READ2)
    check("i5 primer assembles to P5 + index + Read 1 primer",
          il.nebnext_i5_primer("TTGCTTGC"), il.P5 + "TTGCTTGC" + il.TRUSEQ_READ1)
    check("the universal primer carries no index",
          il.NEBNEXT_UNIVERSAL_PRIMER, il.P5 + il.TRUSEQ_READ1[4:])

    check.section("Nextera / Tn5 (lib/nextera.py)")
    check("mosaic end is 19 bp", len(nx.ME), 19)
    check("s5 entry point is 14 nt", len(nx.S5), 14)
    check("s7 entry point is 15 nt", len(nx.S7), 15)
    check("ME_RC is what you read through into", nx.ME_RC, "CTGTCTCTTATACACATCT")
    check("Read 1 primer == s5 + ME", nx.READ1_PRIMER, nx.S5 + nx.ME)
    check("Read 2 primer == s7 + ME", nx.READ2_PRIMER, nx.S7 + nx.ME)
    check("i7 index primer == revcomp(ME) + revcomp(s7)",
          nx.INDEX1_PRIMER, "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC")
    check("i5 index primer == revcomp(ME) + revcomp(s5)",
          nx.INDEX2_PRIMER, "CTGTCTCTTATACACATCTGACGCTGCCGACGA")
    check("N5xx primer == P5 + i5 + s5",
          nx.n5xx_primer("GGGGGGGG"), il.P5 + "GGGGGGGG" + nx.S5)
    check("N7xx primer == P7 + i7 + s7",
          nx.n7xx_primer("GGGGGGGG"), il.P7 + "GGGGGGGG" + nx.S7)
    check("tagmentation leaves a 9-bp gap", nx.TAGMENTATION_GAP, 9)
    check("only the s5/s7 heteroduplex amplifies",
          [o[2] for o in nx.TAGMENTATION_OUTCOMES], [False, False, True])
    check("amplifiable() agrees with the outcome table",
          all(nx.amplifiable(a, b) == amp for a, b, amp, _ in nx.TAGMENTATION_OUTCOMES))
    check("s5/s7 order does not matter", nx.amplifiable("s7", "s5"))

    check.section("reverse transcription / template switching (lib/rt.py)")
    check("SMART / ISPCR handle is 23 nt", len(rt.SMART_HANDLE), 23)
    check("MMLV adds three untemplated C", rt.UNTEMPLATED_TAIL, "CCC")
    check("the TSO G tail pairs with that CCC overhang",
          rt.TSO_G_TAIL.count("G"), len(rt.UNTEMPLATED_TAIL))
    check("the SMART-seq2 TSO locks its terminal G", rt.TSO_G_TAIL_LNA.endswith("+G"))
    check("oligo_dt builds an anchored primer",
          rt.oligo_dt(5, "VN", "HANDLE"), "HANDLETTTTTVN")
    check("oligo_dt can be unanchored", rt.oligo_dt(4, ""), "TTTT")
    check("tso() places the UMI between handle and G tail",
          rt.tso("HANDLE", umi="NNNNNNNN"), "HANDLENNNNNNNNrGrGrG")

    check.section("NEBNext FS sizing table")
    check("longer 37 C incubation gives smaller fragments",
          [il.NEBNEXT_FS_SIZING[t][1] for t in sorted(il.NEBNEXT_FS_SIZING, reverse=True)],
          sorted((v[1] for v in il.NEBNEXT_FS_SIZING.values())))
