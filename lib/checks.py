"""
Shared self-test harness, plus the checks that hold for every protocol.

Each protocol's selftest.py calls `run_common()` and then adds its own checks through
the same `Check` recorder, so one runner prints one report.

The point of these is regression, not discovery: every identity below was verified once
by hand, and is encoded here so that an edit to a shared sequence cannot silently break
a protocol that depends on it.
"""

from __future__ import annotations

import sys

import illumina as il
from chemdraw import Construct, Segment, complement, revcomp, tm


class Check:
    def __init__(self) -> None:
        self.fails: list[str] = []
        self.n = 0

    def __call__(self, label: str, got, want=True) -> None:
        self.n += 1
        ok = (got == want)
        if not ok:
            self.fails.append(f"{label}\n      got:  {got!r}\n      want: {want!r}")
        print(f"  {'PASS' if ok else 'FAIL'}  {label}")

    def section(self, title: str) -> None:
        print(f"\n{title}\n" + "-" * len(title))

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
        print(f"all {self.n} checks passed")
        return 0


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

    check.section("NEBNext FS sizing table")
    check("longer 37 C incubation gives smaller fragments",
          [il.NEBNEXT_FS_SIZING[t][1] for t in sorted(il.NEBNEXT_FS_SIZING, reverse=True)],
          sorted((v[1] for v in il.NEBNEXT_FS_SIZING.values())))
