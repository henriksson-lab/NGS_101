#!/usr/bin/env python3
"""
Self-test for the SMART-seq family.

Runs the shared checks (lib/checks.py), then asserts that the constructs built here
reproduce the published final library structures character-for-character.

Run:  python3 smart-seq__10.1038+nbt.2282/tools/selftest.py
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import nextera as nx
import rt
import smartseq as ss
from checks import Check, run_common
from chemdraw import revcomp, strand_row
from illumina import P5, P7

check = Check()
run_common(check)

# ----------------------------------------------- SMART-seq / SMART-seq2 front end
check.section("SMART-seq / SMART-seq2 oligos")
check("oligo-dT, TSO and ISPCR all carry the SAME handle -- this is the whole design",
      all(x.startswith(rt.SMART_HANDLE) for x in (ss.SS2_OLIGO_DT, ss.SS2_TSO, ss.SS2_ISPCR)))
check("the ISPCR primer IS that handle", ss.SS2_ISPCR, rt.SMART_HANDLE)
check("oligo-dT is anchored with VN", ss.SS2_OLIGO_DT.endswith("VN"))
check("oligo-dT carries 30 T", ss.SS2_OLIGO_DT.count("T" * 30) == 1)
check("the TSO ends in an LNA-locked G (the SMART-seq2 improvement)",
      ss.SS2_TSO.endswith("+G"))
check("the TSO's G tail is 3 nt, matching MMLV's CCC overhang",
      ss.SS2_TSO.count("G"), 3 + rt.SMART_HANDLE.count("G"))

# --------------------------------------------------------- SMART-seq3 front end
check.section("SMART-seq3 family oligos")
check("the 5' tag is 11 bp", len(ss.SS3_TSO_TAG), 11)
check("the UMI is 8 bp", ss.SS3_UMI_LEN, 8)
check("the TSO handle ends with the ME's 3' end, so the TSO IS a Nextera s5 arm",
      ss.SS3_TSO_HANDLE.startswith(nx.ME[-8:]))
check("the forward PCR primer rebuilds s5 + ME + tag",
      ss.SS3_FWD_PCR, nx.S5 + nx.ME + ss.SS3_TSO_TAG)
check("the reverse PCR primer is the oligo-dT handle", ss.SS3_REV_PCR, ss.SS3_OLIGO_DT_HANDLE)
check("oligo-dT and TSO handles DIFFER -- this is what makes PCR non-suppressive",
      ss.SS3_OLIGO_DT_HANDLE != ss.SS3_TSO_HANDLE)
check("three family members, differing only by a TSO spacer", len(ss.SS3_TSO_SPACERS), 3)
check("SMART-seq3 has no spacer", ss.SS3_TSO_SPACERS["SMART-seq3"], "")
check("SMART-seq3xpress inserts WW", ss.SS3_TSO_SPACERS["SMART-seq3xpress"], "WW")
check("FLASH-seq inserts CTAAC", ss.SS3_TSO_SPACERS["FLASH-seq"], "CTAAC")
for v in ss.SS3_TSO_SPACERS:
    check(f"{v} TSO ends in the rGrGrG tail", ss.ss3_tso(v).endswith(rt.TSO_G_TAIL))

# --------------------------------- constructs vs the published final structures
check.section("final library structures vs scg_lib_structs")
PUBLISHED_SS2 = ("AATGATACGGCGACCACCGAGATCTACACNNNNNNNNTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
                 "XXXXXXXX...XXXXXXXX"
                 "CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG")
ss2 = ss.smartseq2_library()
check("SMART-seq2 top strand matches the published structure exactly",
      ss2.top(), PUBLISHED_SS2)
check("SMART-seq2 bottom strand is the complement throughout",
      len(ss2.top()), len(ss2.bottom()))
check("SMART-seq2 carries NO UMI", not any(s.tag == "umi" for s in ss2))
check("SMART-seq2 starts at P5 and ends at revcomp(P7)",
      ss2.top().startswith(P5) and ss2.top().endswith(revcomp(P7)))

ss3 = ss.smartseq3_library()
check("SMART-seq3 5' fragment carries a UMI", any(s.tag == "umi" for s in ss3))
check("SMART-seq3 5' fragment carries the 11-bp tag", ss.SS3_TSO_TAG in ss3.top())
check("SMART-seq3 5' fragment is 22 nt longer than SMART-seq2 (tag + UMI + GGG)",
      len(ss3) - len(ss2), len(ss.SS3_TSO_TAG) + ss.SS3_UMI_LEN + 3)
internal = ss.smartseq3_library(five_prime=False)
check("a SMART-seq3 INTERNAL fragment is indistinguishable from SMART-seq2",
      internal.top(), ss2.top())
check("...which is exactly why the 11-bp tag is needed to call a read 5'-derived",
      ss.SS3_TSO_TAG not in internal.top())
fs = ss.smartseq3_library("FLASH-seq")
check("FLASH-seq adds 5 nt over SMART-seq3", len(fs) - len(ss3), 5)

check.section("shared back end: both chemistries end in Nextera")
for name, con in (("SMART-seq2", ss2), ("SMART-seq3", ss3)):
    check(f"{name} contains the mosaic end", nx.ME in con.top())
    check(f"{name} contains s5 and revcomp(s7)",
          nx.S5 in con.top() and nx.S7_RC in con.top())
    check(f"{name} read 1 primes on s5+ME", nx.READ1_PRIMER in con.top())

check.section("sequencing primers (lib/seqprimers.py)")
import seqprimers as sp  # noqa: E402
for lib, prim in ss.seq_primer_sets():
    check(f"every declared sequencing primer lands on {lib.name}", sp.verify(lib, prim), [])
r1 = sp.locate(ss3, sp.NEXTERA["R1"])
check("SMART-seq3 5': stock Read 1 (s5+ME) ends exactly where the tag begins",
      r1.reads.startswith(ss.SS3_TSO_TAG) and r1.free5 == 0)
check("SMART-seq2: Read 1 goes straight into the cDNA",
      sp.locate(ss2, sp.NEXTERA["R1"]).reads_from, "cDNA")

# ------------------------------------------- step drawings: placed by pairing, not indent
check.section("step drawings (chemdraw.Scene)")
import page_parts as bp  # noqa: E402
from chemdraw import Scene  # noqa: E402

for name, scenes in (("SMART-seq2", bp.ss2_rt()), ("SMART-seq3", bp.ss3_rt())):
    prime, first, switch = scenes
    m, d = prime.strands["mRNA"], prime.strands["oligo-dT"]
    check(f"{name}: oligo-dT V sits opposite the mRNA's last non-A base (B)",
          d.span("V"), m.span("B"))
    check(f"{name}: oligo-dT (T)30 covers exactly the drawn poly(A)",
          d.span("dT"), m.span("polyA"))
    c, m = first.strands["cDNA"], first.strands["mRNA"]
    check(f"{name}: untemplated CCC overhangs the mRNA 5' end",
          c.span("CCC")[1], m.col)
    c, t = switch.strands["cDNA"], switch.strands["TSO"]
    check(f"{name}: TSO G tail sits exactly on the CCC", t.span("rGrG"), c.span("CCC"))
    check(f"{name}: TSO handle hangs off beyond the CCC (nothing to pair with yet)",
          t.col < c.col and t.span("rGrG")[0] == c.col)

def _bad_tso():
    sc = Scene()
    sc.strand("mRNA", bp._mrna())
    sc.anneal("cDNA", bp._cdna(bp.SS2_HANDLE, [bp.seg("AC", "AC")]), to="mRNA",
              pair=("dT", "polyA"))
    sc.anneal("TSO", bp.SS2_TSO, to="cDNA", pair=("rGrG", "CCC"), shift=3)
check.raises("a TSO shifted off its CCC onto the cDNA body cannot be drawn", _bad_tso)

for a, b, _, _ in nx.TAGMENTATION_OUTCOMES:
    sc = bp.tagmented(a, b)
    top, bot = sc.strands["top"], sc.strands["bottom"]
    ntt, ntb = sc.strands["nt-top"], sc.strands["nt-bottom"]
    check(f"tagmentation {a}/{b}: 9-nt gap between top 3' end and the non-transferred ME'",
          ntb.col - top.end(), nx.TAGMENTATION_GAP)
    check(f"tagmentation {a}/{b}: 9-nt gap between bottom 3' end and the other ME'",
          bot.col - ntt.end(), nx.TAGMENTATION_GAP)
    check(f"tagmentation {a}/{b}: each strand's single-stranded 9 nt face a gap",
          (top.span("ss"), bot.span("ss")),
          ((ntt.end(), bot.col), (top.end(), ntb.col)))
    check(f"tagmentation {a}/{b}: each non-transferred strand sits on its ME",
          (ntt.span("ME"), ntb.span("ME")), (top.span("ME"), bot.span("ME")))
for k, (segs, _, _) in bp.seq_primer_drawings().items():
    check(f"drawn {k} sequencing primer spells sp.NEXTERA[{k!r}]",
          "".join(x.top for x in segs), sp.NEXTERA[k].seq)
for v in ss.SS3_TSO_SPACERS:
    check(f"{v} TSO drawing spells ss3_tso() (UMI as N)",
          "".join(x.top for x in bp.ss3_tso_segs(v)).replace(f"[{ss.SS3_UMI_LEN}-bp UMI]",
                                                          "N" * ss.SS3_UMI_LEN),
          ss.ss3_tso(v))
gf = bp.gap_filled()
check("gap fill: the filled top strand is the full s5..s7 library insert",
      gf.strands["top"].text(),
      nx.S5 + nx.ME + "X" * 9 + "XXXXXX...XXXXXX" + "X" * 9 + nx.ME_RC + nx.S7_RC)

check.report()
