#!/usr/bin/env python3
"""
Self-test for the Atrandi SPC + PTA chemistry.

Runs the shared checks (lib/checks.py) and then everything specific to this protocol.
Every assertion here was verified by hand once; encoded so an edit cannot quietly undo it.

Run:  python3 atrandi-wgs__10.1101+2025.06.20.660799/tools/selftest.py
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import atrandi as A
import illumina as il
from checks import Check, run_common
from chemdraw import complement, revcomp, tm

check = Check()
run_common(check)

# ------------------------------------------- Atrandi ligation adapter geometry
check.section("Atrandi ligation adapter (verified base-by-base)")
# The adapter is derived from lib constants; pin that derivation to Atrandi's printed oligos,
# which appear once, as the citation comment in atrandi.py.
import re  # noqa: E402
_printed = dict(re.findall(r"#\s+(Top|Bottom):\s+/\w+/([ACGT*]+)",
                           (HERE / "atrandi.py").read_text()))
check("derived top strand == Atrandi's printed top strand (minus the *)",
      A.LIGADAPT_TOP, _printed["Top"].replace("*", ""))
check("derived bottom strand == Atrandi's printed bottom strand",
      A.LIGADAPT_BOT, _printed["Bottom"])
check("bottom strand's first 12 nt revcomp to the top stem",
      revcomp(A.LIGADAPT_BOT[:12]), A.LIGADAPT_STEM)
check("the stem is the canonical Illumina 12-bp core", A.LIGADAPT_STEM, il.STEM)
check("bottom adapter leaves a 1-nt 3' overhang", len(A.LIGADAPT_BOT) - 12, 1)
check("that overhang is a T (TA ligation to the dA tail)", A.LIGADAPT_BOT[12:], "T")
check("top adapter is 32 nt", len(A.LIGADAPT_TOP), 32)
check("single-stranded arm is 20 nt", len(A.LIGADAPT_ARM), 20)
check("arm revcomps to the first 20 nt of TruSeq Read 1",
      revcomp(A.LIGADAPT_ARM), il.TRUSEQ_READ1[:20])
check("bottom adapter is the last 13 nt of TruSeq Read 1",
      A.LIGADAPT_BOT, il.TRUSEQ_READ1[-13:])
check("post-ligation top strand == Illumina's Read-2 trim sequence",
      "A" + A.LIGADAPT_TOP, il.TRIM_SEEN_IN_READ2)
check("it is NOT a forked adapter: one arm, not two",
      len(A.LIGADAPT_BOT) - 12 + 0, 1)
check("the adapter carries NO Read-2/P7 arm",
      il.TRUSEQ_READ2[:12] not in A.LIGADAPT_TOP + A.LIGADAPT_BOT)

# ------------------------------------------------- Atrandi indexing primers
check.section("Atrandi indexing primers")
check("i7 primer begins with Illumina P7", A.ATRANDI_I7.startswith(il.P7))
check("i7 index is 6 nt, not NEBNext's 8", len(A.ATRANDI_I7_INDEX), 6)
check("i7 3' portion is the first 27 nt of TruSeq Read 2",
      A.ATRANDI_I7_ANNEAL, il.TRUSEQ_READ2[:27])
check("i7 is truncated by exactly 7 nt vs canonical",
      len(il.TRUSEQ_READ2) - len(A.ATRANDI_I7_ANNEAL), 7)
check("P5 primer begins with Illumina P5", A.ATRANDI_P5.startswith(il.P5))
check("P5 3' portion is the first 20 nt of TruSeq Read 1",
      A.ATRANDI_P5_ANNEAL, il.TRUSEQ_READ1[:20])
check("DESIGN RULE: P5 annealing length == adapter arm length, exactly",
      len(A.ATRANDI_P5_ANNEAL), len(A.LIGADAPT_ARM))
check("P5 3' portion revcomps to the adapter arm (primes with zero slack)",
      revcomp(A.ATRANDI_P5_ANNEAL), A.LIGADAPT_ARM)
check("P5 5' tail that PCR copies in is 25 nt", len(A.ATRANDI_P5_TAIL), 25)
check("the 54 C anneal is set by the SHORTER P5 arm, not the i7 arm",
      tm(A.ATRANDI_P5_ANNEAL) < tm(A.ATRANDI_I7_ANNEAL))
check("P5 primer Tm rounds to the protocol's 54 C", round(tm(A.ATRANDI_P5_ANNEAL)), 54)
check("NEB's own primers would NOT work here (too long for the 20-nt arm)",
      len(il.TRUSEQ_READ1) > len(A.LIGADAPT_ARM))

check.section("our indexing primers (ref/our_index_primers.tsv)")
_tsv = [l.split("\t") for l in (HERE.parent / "ref" / "our_index_primers.tsv")
        .read_text().splitlines()[1:]]
_tsv = {r[0]: r for r in _tsv}
for k, idx in A.OUR_I7_INDEX.items():
    row = _tsv[f"{k}_p7"]
    check(f"{k}: model == ordered oligo (minus the * )", A.OUR_I7[k], row[2].replace("*", ""))
    check(f"{k}: index is {A.OUR_I7_INDEX_LEN} nt", len(idx), A.OUR_I7_INDEX_LEN)
    check(f"{k}: sample-sheet i7 is the revcomp of the oligo's index", revcomp(idx), row[4])
    check(f"{k}: only the index differs from Atrandi's i7",
          A.OUR_I7[k].replace(idx, A.ATRANDI_I7_INDEX), A.ATRANDI_I7)
check("the four i7 indices are distinct", len(set(A.OUR_I7_INDEX.values())), 4)
check("our P5 == ordered oligo", A.OUR_P5, _tsv["generic_P5"][2].replace("*", ""))
check("our P5 is Atrandi's P5, unchanged", A.OUR_P5, A.ATRANDI_P5)

check.section("rendered primer lines reconstruct the real oligos")
check("P5 primer line == ATRANDI_P5 (the shared 'ACAC' is not double-counted)",
      A.P5_SEQ + A.ATRANDI_P5[len(A.P5_SEQ):], A.ATRANDI_P5)
check("i7 primer line == ATRANDI_I7",
      A.P7_SEQ + A.ATRANDI_I7_INDEX + A.ATRANDI_I7_ANNEAL, A.ATRANDI_I7)
check("P5 and the Read-1 arm really do overlap by ACAC",
      il.P5.endswith("ACAC") and A.ATRANDI_P5_ANNEAL.startswith("ACAC"))

# -------------------------------------- read layout must match Bascet, independently
check.section("Read 2 layout vs Bascet (crates/bascet-cli/src/barcode/)")
r2 = A.r2_read()
for name, want in (("BC-D", 0), ("BC-C", 12), ("BC-B", 24), ("BC-A", 36)):
    check(f"{name} anchors at R2 offset {want}", r2.offset(name), want)
check("trim_bcread_len == 45  (8+4+8+4+8+4+8+1)", r2.offset("insert"), 45)
check("four barcodes", sum(1 for s in r2 if s.name.startswith("BC-")), 4)
check("three linkers between them",
      sum(1 for s in r2 if s.tag == "r2" and s.placeholder), 3)

# ----------------------------------------------------------- construct integrity
check.section("construct integrity")
lib = A.final_library()
check("top and bottom strands are the same length", len(lib.top()), len(lib.bottom()))
for s in lib:
    if not s.placeholder and s.bottom is None:
        check(f"segment {s.name or '-'!r}: bottom is the complement of top",
              s.bottom_text(), complement(s.top))
check("barcode A placeholder is not complemented to 8 thymines",
      lib.get("BC-A").bottom_text() != "TTTTTTTT")
check("placeholder bottom strand is lowercased, distinct from the dA/dT junction",
      lib.get("BC-A").bottom_text(), "aaaaaaaa")
check("the dA/dT junction base IS complemented (T over A)",
      lib.segments[10].top == "T" and lib.segments[10].bottom_text() == "A")
check("every linker is flagged inferred, for <inf> rendering",
      all(s.inferred for s in lib if s.tag == "r2" and s.placeholder))

# ----------------------------------- PTA amplicon ends (06_pta.md)
check.section("PTA amplicon")
amp = A.pta_amplicon()
check("5' end is the random primer", amp.segments[0].name, "random primer")
check("3' end is the terminator", amp.segments[-1].name, "terminator")
check("primer length is 6-9 nt", 6 <= A.PTA_PRIMER_LEN <= 9)
check("two 3'-terminal phosphorothioate linkages", A.PTA_PRIMER_PS_LINKAGES, 2)
check("both termini carry a note explaining why they matter",
      all(amp.segments[i].note for i in (0, -1)))

# ------------------------------------- the open design question stays switchable
check.section("open design question: 27-nt vs 34-nt round-D arm")
base = len(A.final_library())
A.D_ARM_INCLUDES_CCGATCT = True
grown, alt = len(A.final_library()), A.r2_read()
A.D_ARM_INCLUDES_CCGATCT = False
check("flipping to the 34-nt arm leaves R2 offsets unchanged", alt.offset("BC-A"), 36)
check("...and the construct grows by exactly 7 nt", grown - base, 7)
check("flipping back restores the original length", len(A.final_library()), base)

# ------------------------------------ drawings: every multi-strand panel is a checked Scene
check.section("drawings (chemdraw.Scene verifies every paired column)")
import build_page as B  # noqa: E402
for flag in (False, True):
    A.D_ARM_INCLUDES_CCGATCT = flag
    B.lib = A.final_library()
    try:
        html = "".join(f() for f in (B.step_pta, B.step_endprep, B.step_barcode_a,
                                     B.step_barcode_bcd, B.step_release_frag, B.step_adapter,
                                     B.step_pcr, B.sequencing))
        check(f"all panels build with {'34' if flag else '27'}-nt round-D arm", bool(html))
    except ValueError as e:
        check(f"all panels build with {'34' if flag else '27'}-nt round-D arm: {e}", False)
A.D_ARM_INCLUDES_CCGATCT = False
B.lib = A.final_library()

check.section("P5-side dA junction and sequencing-primer sites")
t = lib.top()
check("the insert's dA is drawn once, between insert and adapter",
      t[t.index("XXX...XXX") + 9:].startswith("A" + A.LIGADAPT_TOP))
check("A + adapter top == Illumina's Read-2 trim sequence, in the library itself",
      A.TRIM_SEEN_IN_R2 in t)
check("TruSeq Read 1 primer is fully templated (its 3'-T pairs the dA)",
      revcomp(il.TRUSEQ_READ1) in t)
_pp = B.pre_pcr_scene().strands
check("pre-PCR bottom strand ends (5'-amino) at the adapter stem",
      _pp["bottom"].text().endswith(complement(A.LIGADAPT_STEM)))
check("...so the 20-nt P5 arm is single-stranded", _pp["top"].end() - _pp["bottom"].end(),
      len(A.LIGADAPT_ARM))
import seqprimers  # noqa: E402
for flag in (False, True):
    A.D_ARM_INCLUDES_CCGATCT = flag
    _l = A.final_library()
    check(f"every declared sequencing primer lands where declared ({'34' if flag else '27'}-nt "
          "arm, lib/seqprimers.py)", seqprimers.verify(_l, A.SEQ_PRIMERS), [])
    check(f"stock Read 2 declared mismatching iff the arm is 27 nt ({flag})",
          bool([p for p in A.SEQ_PRIMERS if p.role == "Read 2"][0].expect_mismatch), not flag)
A.D_ARM_INCLUDES_CCGATCT = False
_i2 = seqprimers.locate(lib, [p for p in A.SEQ_PRIMERS if p.role.startswith("Index 2")][0])
check("no i5 index: the Index 2 read runs straight into P5", _i2.reads_from, "Illumina P5")
check("stock Read 2 primer is fully templated only with the 34-nt arm",
      (il.TRUSEQ_READ2 in t), A.D_ARM_INCLUDES_CCGATCT)

check.report()
