#!/usr/bin/env python3
"""
Self-test for Small-seq.

Runs the shared checks (lib/checks.py), then asserts every fact that was verified by
hand against Hagemann-Jensen 2018 (Nat Protoc 13:2407). The interesting ones are the
relationships that must hold if the model is right and that nobody fed the model:

  * RTP is the exact reverse complement of RA3 -- the whole digestion step depends on it
  * RP1 is P5 (from lib/illumina.py, a vendor document) plus the first 21 nt of RA5
  * RP1 stops 5 nt short of RA5's 3' end, so RP1 cannot be the read-1 primer
  * the nearest-neighbour Tm of RTP lands on the published PCR-1 annealing temperature
  * the unpublished SRX arm is bracketed, to +/-5 bp, by two independent published sizes
    (the miRNA size-selection window and the 5.8S rRNA peak)

Run:  python3 small-seq__10.1038+nbt.3701/tools/selftest.py
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import nextera as nx
import rt
import smallseq as sm
from checks import Check, run_common
from chemdraw import complement, revcomp, tm
from illumina import P5, P7, P7_RC, TRUSEQ_READ1, TRUSEQ_READ2

check = Check()
run_common(check)

# ------------------------------------------------------------------ the oligo set
check.section("Small-seq oligos (Reagent setup, p. 2414)")
check("RA5 adapter proper is 26 nt", len(sm.RA5_HANDLE), 26)
check("RA5 as ordered is 36 nt (26 + 8 UMI + CA)", len(sm.ra5_oligo()), 36)
check("RA3 is 21 nt", len(sm.RA3), 21)
check("RTP is 21 nt", len(sm.RTP), 21)
check("RP1 is 50 nt", len(sm.RP1), 50)
check("the 5.8S masking oligo is 76 nt", len(sm.MASK_58S), 76)
check("every published oligo is real DNA/RNA, no placeholders smuggled in",
      all(set(s) <= set("ACGT") for s in (sm.RA5_HANDLE, sm.RA3, sm.RTP, sm.RP1,
                                          sm.MASK_58S)))
check("the UMI is 8 nt", sm.RA5_UMI_LEN, 8)
check("UMI bases are H (rA/rU/rC), NOT N -- this is deliberate, to reduce mispriming",
      sm.RA5_UMI_BASE == "H" and "N" not in sm.ra5_oligo())
check("8 H positions give the 6,561 UMIs the paper quotes", sm.RA5_UMI_DIVERSITY, 6561)
check("...which is strictly fewer than 4^8 would give", sm.RA5_UMI_DIVERSITY < 4 ** 8)
check("the UMI is followed by the CA linker, at the adapter's 3' end",
      sm.ra5_oligo().endswith(sm.RA5_UMI_BASE * 8 + "CA"))

# ---------------------------------------------------------- how the oligos interlock
check.section("how the four oligos interlock")
check("RTP is the EXACT reverse complement of RA3 -- the digestion step needs this",
      sm.RTP, revcomp(sm.RA3))
check("...so RTP + RA3 form a blunt 21-bp duplex, lambda exonuclease's substrate",
      len(sm.RTP), len(sm.RA3))
check("RP1 (built as P5 + RA5[:21]) equals the string printed in the paper, kept verbatim "
      "in the comment in smallseq.py",
      f"RP1 ({sm.RP1})" in (HERE / "smallseq.py").read_text())
check("RP1's 5' 29 nt are the canonical P5 (from lib/illumina.py)",
      sm.RP1[:len(P5)], P5)
check("RP1's 3' 21 nt anneal inside RA5", sm.RP1[len(P5):], sm.RA5_HANDLE[:21])
check("RP1 anneals over 21 nt", sm.RP1_ANNEAL_NT, len(sm.RP1) - len(P5))
check("RP1 stops 5 nt short of RA5's 3' end -- those 5 are CGATC",
      sm.RA5_HANDLE[sm.RP1_ANNEAL_NT:], "CGATC")
check("so RP1 is NOT a read-1 primer: priming on it would read CGATC before the UMI",
      not sm.RP1.endswith(sm.RA5_HANDLE))
check("RA3 begins TGG, which is why the pipeline script is remove_reads_with_genomic_TGG",
      sm.RA3[:3], "TGG")
check("the RT primer's 3' end is where the cDNA starts, i.e. on RA3",
      revcomp(sm.RTP) in sm.RA3)

# ----------------------------------------------------------------- thermodynamics
check.section("thermodynamics (lib/chemdraw.tm, SantaLucia 1998)")
rtp_tm = tm(sm.RTP)
check("RTP / RA3 melt at the same temperature -- they are the same duplex",
      round(tm(sm.RA3), 3), round(rtp_tm, 3))
check("RTP melts at 60 C to within a degree", 59.0 < rtp_tm < 61.0)
check("...which is exactly the published PCR-1 annealing temperature",
      abs(sm.PCR1_ANNEAL_C - rtp_tm) < 1.0)
check("full-length RP1 melts far higher, so it is never limiting",
      tm(sm.RP1) > rtp_tm + 10)
check("the published PCR-2 annealing temperature sits below the full-length Tm",
      sm.PCR2_ANNEAL_C < tm(sm.RP1))
check("PCR 2 anneals hotter than PCR 1, because by then both primers bind full-length",
      sm.PCR2_ANNEAL_C > sm.PCR1_ANNEAL_C)
check("RP1's 21-nt annealing foot alone melts BELOW 60 C -- cycle 1 is the hard one",
      tm(sm.RP1[len(P5):]) < sm.PCR1_ANNEAL_C)
check("the masking oligo stays annealed through the 72 C lysis",
      tm(sm.MASK_58S) > sm.LYSIS_TEMP_C)
check("...with better than 10 C of margin", tm(sm.MASK_58S) > sm.LYSIS_TEMP_C + 10)

# ------------------------------------------------------------------ the mask target
check.section("5.8S rRNA masking oligonucleotide")
check("its target is its reverse complement, and round-trips",
      revcomp(sm.MASK_58S_TARGET), sm.MASK_58S)
check("mask and target are the same length", len(sm.MASK_58S_TARGET), len(sm.MASK_58S))
check("it covers about half of a 156-nt 5.8S rRNA",
      0.45 < len(sm.MASK_58S) / sm.RRNA_58S_NT < 0.55)
check("it is not an Illumina sequence in disguise",
      not any(x in sm.MASK_58S for x in (P5, P7, TRUSEQ_READ1, TRUSEQ_READ2, nx.ME)))
check("it does not anneal to RA3, so it cannot sequester the 3' adapter",
      revcomp(sm.RA3) not in sm.MASK_58S and sm.RA3 not in sm.MASK_58S)

# --------------------------------------------------------------- read-1 arithmetic
check.section("read layout")
lib = sm.library()
check("the UMI starts at offset 55, immediately after P5 + RA5",
      lib.offset("UMI"), len(P5) + len(sm.RA5_HANDLE))
check("sequencing starts from the UMI, so 10 nt are trimmed before mapping",
      sm.TRIM_5P, 10)
check("a 51-bp read leaves 41 bp, the paper's own figure",
      sm.READ_LEN - sm.TRIM_5P, 41)
check("...and 41 is the precursor cut-off", sm.READ_LEN - sm.TRIM_5P, sm.PRECURSOR_MIN_NT)
check("a small RNA is anything shorter than that", sm.SMALLRNA_MAX_NT + 1,
      sm.PRECURSOR_MIN_NT)
check("the paper's worked 76-bp example gives 66 bp, as stated", 76 - sm.TRIM_5P, 66)
check("the analysis window is 18-40 nt", (sm.SMALLRNA_MIN_NT, sm.SMALLRNA_MAX_NT),
      (18, 40))
check("read 1 must be primed from a site ending at the last base of RA5",
      (lib.top()[:lib.offset("UMI")]).endswith(sm.RA5_HANDLE))
check("that site is a suffix of P5 + RA5 and of nothing else",
      lib.top()[:lib.offset("UMI")], P5 + sm.RA5_HANDLE)
check("the TruSeq DNA read-1 primer site is ABSENT from this library",
      TRUSEQ_READ1 not in lib.top() and revcomp(TRUSEQ_READ1) not in lib.top())
check("so is the TruSeq DNA read-2 primer site",
      TRUSEQ_READ2 not in lib.top() and revcomp(TRUSEQ_READ2) not in lib.top())

# -------------------------------------------------------------- the ligated molecule
import seqprimers  # noqa: E402
check("every declared sequencing primer lands where declared (lib/seqprimers.py)",
      seqprimers.verify(lib, sm.SEQ_PRIMERS), [])
r1 = seqprimers.locate(lib, sm.SEQ_PRIMERS[0])
check("the read-1 site ends exactly where the UMI begins", r1.end, lib.offset("UMI"))
rp1 = seqprimers.locate(lib, sm.SEQ_PRIMERS[1])
check("RP1 as a read-1 primer would read CGATC first", rp1.reads[:5], "CGATC")

check.section("sequential ligation")
lig = sm.ligated_rna()
check("after both ligations the molecule is RA5 + UMI + CA + insert + RA3",
      [s.name for s in lig], ["RA5", "UMI", "CA", "small RNA", "RA3"])
check("no P5, no P7, no index yet -- adapters only",
      not any(x in lig.top() for x in (P5, P7, P7_RC)))
check("the 3' adapter is ligated FIRST (T4 RNA ligase 2 truncated, no ATP)",
      sm.ENZYMES[0][1].startswith("T4 RNA ligase 2"))
check("the 5' adapter is ligated LAST (T4 RNA ligase 1, with ATP)",
      [e[1] for e in sm.ENZYMES].index("T4 RNA ligase 1 + ATP") > 0)
check("deadenylase comes before lambda exonuclease -- it makes the 5'-phosphate",
      [e[1] for e in sm.ENZYMES].index("5' deadenylase")
      < [e[1] for e in sm.ENZYMES].index("lambda exonuclease"))
check("the ligated molecule is 69 nt for a 22-nt miRNA",
      len(lig), 26 + 8 + 2 + 22 + 21)

# ------------------------------------------------------------------- final library
check.section("the final library")
check("left arm (P5 + RA5 + UMI + CA + RA3) is 86 bp", sm.LEFT_ARM_NT, 86)
check("segment table reproduces that arithmetic",
      sum(len(lib.get(n)) for n in ("P5", "RA5", "UMI", "CA", "RA3")), 86)
check("the library starts with P5", lib.top().startswith(P5))
check("the library ends with revcomp(P7)", lib.top().endswith(P7_RC))
check("bottom strand is the same length as the top", len(lib.bottom()), len(lib.top()))
check("the UMI placeholder is never complemented as if it were bases",
      "h" * 8 in lib.bottom() and "T" * 8 not in lib.bottom())
check("the published left arm is drawn as known, the SRX arm as inferred",
      [s.inferred for s in lib],
      [False, False, False, False, False, False, True, True, True])
check("there is no cell barcode -- one cell per well, identified by its i7 index",
      not any(s.name == "cell barcode" for s in lib))
check("exactly one UMI segment", sum(s.tag == "umi" for s in lib), 1)

check.section("this is neither a tagmentation nor a template-switching library")
check("no Nextera mosaic end", nx.ME not in lib.top() and nx.ME_RC not in lib.top())
check("no s5 / s7 entry point", nx.S5 not in lib.top() and nx.S7 not in lib.top())
check("no SMART / ISPCR handle", rt.SMART_HANDLE not in lib.top())
check("no oligo-dT, because small RNAs are not polyadenylated", "TTTTTTTT" not in lib.top())
check("no rGrGrG tail: the 5' handle is ligated, not template-switched",
      "GGG" not in sm.ra5_oligo())

# ---------------------------------------------------- sizes, and the SRX open question
check.section("library sizes, and bracketing the unpublished SRX arm")
check("SRX's own sequence is deliberately absent", sm.SRX_SEQUENCE, None)
check("its index is 8 bp", sm.SRX_INDEX_NT, 8)
check("192 of them, i.e. two 96-well plates", sm.SRX_N_INDICES, 2 * 96)
check("the arm it adds is P7 + index + linker", sm.SRX_RIGHT_ARM_NT,
      len(P7) + sm.SRX_INDEX_NT + sm.SRX_LINKER_NT)

lo_mi = sm.SIZE_SELECT_TARGET_BP[0] - sm.LEFT_ARM_NT - sm.TYPICAL_MIRNA_NT
hi_mi = sm.SIZE_SELECT_TARGET_BP[1] - sm.LEFT_ARM_NT - sm.TYPICAL_MIRNA_NT
lo_rr = sm.RRNA_PEAK_BP[0] - sm.LEFT_ARM_NT - sm.RRNA_58S_NT
hi_rr = sm.RRNA_PEAK_BP[1] - sm.LEFT_ARM_NT - sm.RRNA_58S_NT
check("the miRNA size-selection window alone brackets the arm to 37-47 nt",
      (lo_mi, hi_mi), (37, 47))
check("the 5.8S rRNA peak independently brackets it to 28-48 nt", (lo_rr, hi_rr), (28, 48))
check("the two windows overlap at all -- they did not have to",
      max(lo_mi, lo_rr) <= min(hi_mi, hi_rr))
check("the modelled arm sits inside BOTH published windows",
      lo_mi <= sm.SRX_RIGHT_ARM_NT <= hi_mi and lo_rr <= sm.SRX_RIGHT_ARM_NT <= hi_rr)

check("a 22-nt miRNA library is 152 bp", sm.library_bp(sm.TYPICAL_MIRNA_NT), 152)
check("...which lands in the published 145-155 bp cut window",
      sm.SIZE_SELECT_TARGET_BP[0] <= sm.library_bp(22) <= sm.SIZE_SELECT_TARGET_BP[1])
check("...and inside the Pippin 130-160 bp window too",
      sm.PIPPIN_WINDOW_BP[0] <= sm.library_bp(22) <= sm.PIPPIN_WINDOW_BP[1])
check("a full-length 5.8S rRNA library would be 286 bp",
      sm.library_bp(sm.RRNA_58S_NT), 286)
check("...landing in the 270-290 bp peak seen when the mask is left out",
      sm.RRNA_PEAK_BP[0] <= sm.library_bp(sm.RRNA_58S_NT) <= sm.RRNA_PEAK_BP[1])
check("the whole insert range 18-40 nt fits the 100-300 bp Bioanalyzer profile",
      all(sm.LIBRARY_PROFILE_BP[0] <= sm.library_bp(n) <= sm.LIBRARY_PROFILE_BP[1]
          for n in range(sm.SMALLRNA_MIN_NT, sm.SMALLRNA_MAX_NT + 1)))
check("the adapter dimer is 130 bp", sm.ADAPTER_DIMER_BP, 130)
check("...which is INSIDE the Pippin window, so size selection cannot remove it",
      sm.PIPPIN_WINDOW_BP[0] <= sm.ADAPTER_DIMER_BP)
check("...and inside the 120-200 bp gel cut as well",
      sm.GEL_CUT_BP[0] <= sm.ADAPTER_DIMER_BP <= sm.GEL_CUT_BP[1])
check("the shortest real library is only 18 bp bigger than the dimer",
      sm.library_bp(sm.SMALLRNA_MIN_NT) - sm.ADAPTER_DIMER_BP, 18)
check("library length is linear in insert length, by construction",
      sm.library_bp(40) - sm.library_bp(18), 22)
check("the construct's own length agrees with library_bp()",
      len(sm.library(22)), sm.library_bp(22))

check.section("protocol constants")
check("lysis is 20 min at 72 C", (sm.LYSIS_TEMP_C, sm.LYSIS_MIN), (72, 20))
check("3' ligation is 6 h at 30 C, then overnight at 4 C",
      (sm.LIG3_TEMP_C, sm.LIG3_HOURS), (30, 6))
check("5' ligation is 1 h at 37 C", (sm.LIG5_TEMP_C, sm.LIG5_MIN), (37, 60))
check("RT is 1 h at 42 C with SuperScript II", (sm.RT_TEMP_C, sm.RT_MIN), (42, 60))
check("both PCRs run 13 cycles", (sm.PCR1_CYCLES, sm.PCR2_CYCLES), (13, 13))
check("PCR 1 uses RP1 plus the leftover RT primer", sm.PCR1_PRIMERS[0], "RP1")
check("PCR 2 uses the same forward primer and an SRX index primer",
      sm.PCR2_PRIMERS[0], "RP1")
check("SRX is supplied at 50 uM", sm.SRX_CONC_UM, 50)

# ------------------------------------------------------------ the page's drawings
check.section("drawings: every multi-strand panel is a chemdraw.Scene")
import build_page as bp  # noqa: E402

for fn in (bp.scene_rtp_ra3, bp.scene_rp1, bp.scene_mask, bp.scene_lig3, bp.scene_digest,
           bp.scene_lig5, bp.scene_rt, bp.scene_pcr1, bp.scene_pcr2, bp.scene_read1):
    try:
        fn()
        ok = True
    except ValueError:
        ok = False
    check(f"{fn.__name__}: every drawn column base-pairs", ok)

src = (HERE / "build_page.py").read_text()
check("build_page.py positions no strand by hand (no indent=, no markers([...]))",
      "indent=" not in src and "markers([" not in src)

srx = bp.srx_segs()
check("SRX is drawn with an RA3-annealing 3' end (fixed: it used to be the top strand "
      "reversed, with nothing that could anneal)",
      "".join(s.top for s in srx).endswith(sm.RTP))
check("...and its 5' end is P7, as the clustering requires",
      "".join(s.top for s in srx).startswith(P7))
a, _ = bp.scene_pcr2()
st, top = a.strands["SRX"], a.strands["top"]
check("in PCR 2, SRX's 3' end lies over RA3 and its link/i7/P7 tail overhangs the template",
      st.span("RA3'") == top.span("RA3") and st.span("SRX link'")[0] >= top.end())
c, _ = bp.scene_pcr1()
check("in PCR 1, RP1 anneals to the cDNA (antiparallel), not alongside the RNA",
      c.strands["RP1"].partner == "cDNA" and c.strands["RP1"].rev != c.strands["cDNA"].rev)
r = bp.scene_read1()
check("read 1 primer anneals to the bottom strand and ends where the UMI begins",
      r.strands["read 1 primer"].partner == "bottom"
      and r.strands["read 1 primer"].end() == r.strands["top"].span("UMI")[0])

check.report()
