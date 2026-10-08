#!/usr/bin/env python3
"""Self-test for CRISPR-UMI (Schmierer). Run: python3 crispr-umi-schmierer__10.15252+msb.20177834/tools/selftest.py"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))
sys.path.insert(0, str(HERE.parents[1] / "crispr-mip__10.1101+2024.03.28.587082" / "tools"))

import crisprmip as cm
import crisprumi as cu
import michlits as mi
from checks import Check, Source, have, run_common
from crispr import SCAFFOLD_V1, clone_guide
import seqprimers as sp
import illumina as il
from illumina import TRUSEQ_READ2
from chemdraw import Scene, Segment, complement_segments, revcomp
from plasmid import Plasmid, amplify, find_both, read_genbank

# ------------------------------------------------- third-party material this suite reads
# None of it is committed (it is Addgene's / the depositors' annotation, not ours to
# redistribute) -- see ref/MANIFEST.md. Each group of checks that reads one of these files
# is wrapped in `have(...)`, so on a fresh clone it SKIPS and says where to get the file
# instead of failing as though the chemistry were wrong. Everything else still runs.
_SCREEN = (HERE.parents[1] / "lenticrispr-v1-screening__10.1126+science.1247005"
           / "ref" / "plasmids")

# cu.PARENT_MAP is already whichever candidate map is on disk, or the preferred one if none.
PARENT = Source(cu.PARENT_MAP, cu.PARENT_MAP_ORIGIN)
SCREEN_MAPS = [
    ("lentiCRISPR v1", Source(_SCREEN / "addgene-49535_lentiCRISPRv1.gb",
                              "Addgene #49535 (lentiCRISPR v1), full sequence")),
    ("lentiCRISPRv2", Source(_SCREEN / "addgene-52961_lentiCRISPRv2.gb",
                             "Addgene #52961 (lentiCRISPRv2), full sequence")),
    ("lentiGuide-Puro", Source(_SCREEN / "addgene-52963_lentiGuide-Puro.gb",
                               "Addgene #52963 (lentiGuide-Puro), full sequence")),
]
UMI_FA = {name: Source(HERE.parents[0] / "ref" / f"addgene-{meta['addgene']}_{name}.fa",
                       f"Addgene #{meta['addgene']} ({name}), full sequence")
          for name, meta in mi.VECTORS.items()}

SPACER = cu.DEMO_SPACER
check = Check()
run_common(check)

check.section("the AU-flip is two substitutions")
diff = [(i, a, b) for i, (a, b) in enumerate(zip(SCAFFOLD_V1, cu.SCAFFOLD_AU_FLIP)) if a != b]
check("exactly two bases differ from the original scaffold", len(diff), 2)
check("at positions 5 and 26 (1-based)", [d[0] + 1 for d in diff], [5, 26])
check("T->A then A->T", [(d[1], d[2]) for d in diff], [("T", "A"), ("A", "T")])
check("same length as the original", len(cu.SCAFFOLD_AU_FLIP), len(SCAFFOLD_V1))
check("the original scaffold head does NOT match the AU-flip vector",
      not cu.SCAFFOLD_AU_FLIP.startswith(SCAFFOLD_V1[:24]))
check("the scaffold 3' end is untouched, so primers there still work",
      cu.SCAFFOLD_AU_FLIP[-30:], SCAFFOLD_V1[-30:])

check.section("the library insert")
ins = cu.insert()
check("insert contains the U6 overlap ending on the +1 G", cu.U6_OVERLAP.endswith("ACCG"))
check("insert carries the AU-flip scaffold, not the original",
      cu.SCAFFOLD_AU_FLIP in ins and SCAFFOLD_V1 not in ins)
check("RSL is 6 bp", cu.RSL_LEN, 6)
check("the Illumina adapter is built INTO the construct", cu.ILLUMINA_ADAPTER in ins)
check("that adapter reverse-complements to the TruSeq Read 2 primer",
      revcomp(cu.ILLUMINA_ADAPTER), TRUSEQ_READ2[:len(cu.ILLUMINA_ADAPTER)])
check("the RSL sits immediately 3' of it -- i.e. in the i7 index position",
      ins.index("N" * cu.RSL_LEN, ins.index(cu.ILLUMINA_ADAPTER)),
      ins.index(cu.ILLUMINA_ADAPTER) + len(cu.ILLUMINA_ADAPTER))
check("three named analyses: TCA, LDA, IRA", [a[0] for a in cu.ANALYSES],
      ["TCA", "LDA", "IRA"])

check.section("vs the parent vector, computed from the real map")
if have(check, PARENT, label="the parent vector and its PCR1 sizes"):
    p = cu.parent_vector()
    check("the parent map carries the ORIGINAL scaffold",
          p.seq.find(SCAFFOLD_V1 + cu.TERMINATOR) >= 0)
    mod = cu.rebuild_vector()
    check("the edit adds 39 bp (33-nt Illumina adapter + 6-nt RSL)", len(mod) - len(p), 39)
    check("after the edit the original scaffold head is gone",
          not find_both(mod.seq, SCAFFOLD_V1[:24], True))
    check("...and the AU-flip head is there",
          bool(find_both(mod.seq, cu.SCAFFOLD_AU_FLIP[:24], True)))
    check("the Illumina Read-2 adapter is now in the vector",
          bool(find_both(mod.seq, cu.ILLUMINA_ADAPTER, True)))
    check("a scaffold-3' primer still binds (Joung's readout site)",
          bool(find_both(mod.seq, revcomp(SCAFFOLD_V1[42:65]), True)))
    check("the CRISPR-MIP ligation arm is absent here too",
          not find_both(mod.seq, cm.LIG_ARM, True))

    cl_parent, cl_mod = clone_guide(p, SPACER), cu.rebuild_vector(cloned=True)
    a_par = amplify(cl_parent, cu.OUR_PCR1_FW, cu.OUR_PCR1_REV, min_anneal=18)
    a_mod = amplify(cl_mod, cu.OUR_PCR1_FW, cu.OUR_PCR1_REV, min_anneal=18)
    check("the parent gives a 557 bp PCR1", a_par[0].length, 557)
    check("the AU-flip vector gives 596 bp -- 39 bp longer", a_mod[0].length, 596)
    check("so a size difference alone distinguishes the two vectors",
          a_mod[0].length - a_par[0].length, 39)

check.section("padlock design for this vector")
from padlock import Padlock, capture
from design_padlock import schmierer_vector
from chemdraw import tm
check("the CRISPR-MIP backbone's Read-2 site CONTAINS the vector's adapter -- "
      "so a capture spanning it would duplicate a primer site",
      cu.ILLUMINA_ADAPTER in cm.READ2_SITE)
check("the guide-to-RSL span is irreducible: scaffold + terminator + adapter",
      len(cu.SCAFFOLD_AU_FLIP) + len(cu.TERMINATOR) + len(cu.ILLUMINA_ADAPTER), 115)
check("so capturing both needs at least 141 nt",
      20 + 115 + cu.RSL_LEN, cu.RSL_DISTANCE)
check("...which is 1.26x the gap the published probe is known to fill",
      round(cu.RSL_DISTANCE / cu.VALIDATED_GAP, 2), 1.26)
check("the guide-only option is SHORTER than the validated gap",
      cu.PADLOCK_GAP_UNIVERSAL < cu.VALIDATED_GAP)
check("...the ligation arm is hotter than the extension arm, by 3-6 C",
      3 <= tm(cu.PADLOCK_LIG_ARM) - tm(cm.EXT_ARM) <= 6)
check("universal arm sits in the scaffold, clear of both AU-flip positions",
      SCAFFOLD_V1.find(cu.PADLOCK_LIG_ARM_UNIVERSAL) > 25)
check("...and is present in BOTH scaffold variants",
      cu.PADLOCK_LIG_ARM_UNIVERSAL in SCAFFOLD_V1
      and cu.PADLOCK_LIG_ARM_UNIVERSAL in cu.SCAFFOLD_AU_FLIP)

rec = Padlock("rec", ext_arm=cm.EXT_ARM, lig_arm=cu.PADLOCK_LIG_ARM, backbone=cm.backbone())
uni = Padlock("uni", ext_arm=cm.EXT_ARM, lig_arm=cu.PADLOCK_LIG_ARM_UNIVERSAL,
              backbone=cm.backbone())
if have(check, PARENT, label="capture of both probes on the rebuilt vector"):
    sv = schmierer_vector()
    c = capture(sv, rec)
    check("option A: exactly one capture site", len(c), 1)
    check(f"...gap is {cu.PADLOCK_GAP} nt", c[0].gap, cu.PADLOCK_GAP)
    check(f"...circle is {cu.PADLOCK_CIRCLE} nt", len(c[0].circle), cu.PADLOCK_CIRCLE)
    check("...the capture contains the spacer", SPACER in c[0].fill)
    check("...and the RSL", "N" * cu.RSL_LEN in c[0].fill)
    check("...the existing P7_tracrRNA site still lies inside the capture",
          revcomp(cm.P7_TRACR_REV[24:]) in c[0].fill)

    cs = capture(sv, uni)
    check("universal arm: same gap on the Schmierer vector too",
          len(cs) == 1 and cs[0].gap == cu.PADLOCK_GAP_UNIVERSAL)
    check("...but it does NOT capture the RSL", "N" * cu.RSL_LEN not in cs[0].fill)
    check("...and its capture is too short for the existing PCR pair",
          revcomp(cm.P7_TRACR_REV[24:]) not in cs[0].fill)

# The point of the universal arm is that it works on the OTHER screening vectors unchanged,
# so this one needs their maps -- kept with the lentiCRISPR screen notes, one directory over.
for nm, src in SCREEN_MAPS:
    if not have(check, src, label=f"universal arm on {nm}"):
        continue
    v = clone_guide(read_genbank(str(src.path)), SPACER)
    cc = capture(v, uni)
    check(f"universal arm: one capture on {nm}, gap {cu.PADLOCK_GAP_UNIVERSAL}",
          len(cc) == 1 and cc[0].gap == cu.PADLOCK_GAP_UNIVERSAL)

check.section("drafted protocol (02_padlock_protocol.md)")
bbA = cu.probe_A_backbone()
check("option A backbone carries NO Read-2 site", cm.READ2_SITE not in bbA)
check("...so the amplicon cannot contain two copies of it",
      (cu.ILLUMINA_ADAPTER in bbA) is False)
prA = Padlock("A", ext_arm=cm.EXT_ARM, lig_arm=cu.PADLOCK_LIG_ARM, backbone=bbA)
check(f"probe A is {cu.PROBE_LEN_A} nt -- ext arm + backbone + lig arm, which moves with "
      "the CRISPR-MIP backbone's own lengths", len(prA), cu.PROBE_LEN_A)
check("the custom seq primer is the 3' 20 nt of the vector's built-in adapter",
      cu.ILLUMINA_ADAPTER.endswith(cu.CUSTOM_SEQ_PRIMER_A))
check("the reverse primer's 3' end is the cooler of the pair -- flagged in the protocol",
      tm(cu.PCR_REV_A[24:]) < tm(cu.PCR_FWD_A[29:]))
if have(check, PARENT, label="probe A's circle, library and cycle positions"):
    cA = capture(schmierer_vector(), prA)
    check(f"probe A: one capture, {cu.CIRCLE_LEN_A}-nt circle",
          len(cA) == 1 and len(cA[0].circle) == cu.CIRCLE_LEN_A)
    ampA = amplify(cA[0].circle, cu.PCR_FWD_A, cu.PCR_REV_A, min_anneal=15)
    check(f"probe A library is {cu.LIBRARY_LEN_A} bp", ampA[0].length, cu.LIBRARY_LEN_A)
    check(f"...which is the circle less the {cu.OUTSIDE_PRIMERS_A}-nt arc outside the primer "
          "pair, plus the non-templated P5 and P7 tails",
          len(cA[0].circle) - (ampA[0].length - len(il.P5) - len(il.P7)),
          cu.OUTSIDE_PRIMERS_A)
    sA = ampA[0].seq
    check("the guide sits at read-1 cycles 21-40",
          sA.find(SPACER) - (sA.find(cm.READ1_SITE) + len(cm.READ1_SITE)), cu.SPACER_LEN)
    check("reading from it, the RSL is at cycles 1-6",
          sA.find("N" * cu.RSL_LEN)
          - (sA.find(cu.CUSTOM_SEQ_PRIMER_A) + len(cu.CUSTOM_SEQ_PRIMER_A)), 0)
    check("...and the UMI within 50 cycles",
          sA.find(cm.UMI) + len(cm.UMI)
          - (sA.find(cu.CUSTOM_SEQ_PRIMER_A) + len(cu.CUSTOM_SEQ_PRIMER_A)) <= 50)




# ---------------------------------------------------------------- Michlits
from plasmid import Plasmid as _P

check.section("CRISPR-UMI (Michlits)")
check("barcode is 10 nt, experimental index 6", (mi.BARCODE_LEN, mi.EXP_INDEX_LEN), (10, 6))
check("the two labels differ in length, so a read of each identifies which is which",
      mi.BARCODE_LEN != mi.EXP_INDEX_LEN)
check("the vector's P7 differs from canonical at exactly one base",
      sum(a != b for a, b in zip(mi.P7_RC_CANONICAL, mi.P7_RC_IN_VECTOR)), 1)
check("...and that change removes a BbsI site",
      "GAAGAC" in mi.P7_RC_CANONICAL or "GTCTTC" in mi.P7_RC_CANONICAL)
check("...which the readout primer restores, carrying canonical P7",
      mi.PCR_REV.startswith(il.P7))
check("the forward readout primer is P5 + 6-nt index + a U6 anneal",
      mi.PCR_FWD.startswith(il.P5)
      and "N" * mi.EXP_INDEX_LEN in mi.PCR_FWD)
check("guide exclusions cover both Golden Gate enzymes",
      "GAAGAC" in mi.FORBIDDEN_IN_GUIDE and "CGTCTC" in mi.FORBIDDEN_IN_GUIDE)
check("a guide containing a BbsI site is rejected", not mi.guide_allowed("ACGAAGACGTACGTACGTAC"))
check("a clean guide is accepted", mi.guide_allowed(SPACER))
check("a guide ending CTCGA is rejected", not mi.guide_allowed("ATCGATCGATCGATCGCTCGA"))

check("the enrichment is 3-4 orders of magnitude", mi.ENRICHMENT_FOLD, (1000, 10000))


def _fasta(src: Source) -> str:
    return "".join(l.strip() for l in src.path.read_text().splitlines()
                   if not l.startswith(">"))


for name, meta in mi.VECTORS.items():
    src = UMI_FA[name]
    if not have(check, src, label=f"{name} as deposited"):
        continue
    p = _P(name, _fasta(src), True, [])
    check(f"{name}: {meta['bp']} bp as deposited", len(p), meta["bp"])
    check(f"{name}: carries the built-in Illumina i7 site",
          bool(find_both(p.seq, mi.ILLUMINA_I7_SITE, True)))
    check(f"{name}: carries the BbsI-modified P7, not the canonical one",
          bool(find_both(p.seq, mi.P7_RC_IN_VECTOR, True))
          and not find_both(p.seq, mi.P7_RC_CANONICAL, True))
    check(f"{name}: the CRISPR-MIP extension arm (U6) is present",
          bool(find_both(p.seq, cm.EXT_ARM, True)))
    check(f"{name}: the CRISPR-MIP ligation arm is NOT",
          not find_both(p.seq, cm.LIG_ARM, True))

if have(check, UMI_FA["pLenti-UMI"], label="the PacI enrichment fragment"):
    pl = _P("pLenti-UMI", _fasta(UMI_FA["pLenti-UMI"]), True, [])
    pac = sorted(x[0] for x in find_both(pl.seq, mi.PACI, True))
    check("pLenti-UMI has exactly 2 PacI sites", len(pac), 2)
    check("...giving a fragment within 10 bp of the paper's 589",
          abs((pac[1] - pac[0]) - mi.PACI_FRAGMENT_BP) <= 10)

check.section("the PUBLISHED readout, reproduced on the rebuilt vector")
check("20 + 6 + 6 cycles, as published", (cu.READ1_CYCLES, cu.I5_CYCLES, cu.I7_CYCLES),
      (20, 6, 6))
check("Read 1 covers the whole spacer", cu.READ1_CYCLES, 20)
lib = None
if have(check, PARENT, label="the three nested PCRs and the 288-bp product"):
    sv2 = schmierer_vector()
    r1 = amplify(sv2, cu.PCR1_FW, cu.PCR1_REV, min_anneal=18)
    check("PCR1 gives a single product", len(r1), 1)
    r2 = amplify(Plasmid("r1", r1[0].seq, False, []), cu.PCR2_FW, cu.PCR2_REV, min_anneal=18)
    check("PCR2 nests inside PCR1", len(r2) == 1 and r2[0].length < r1[0].length)
    r3 = amplify(Plasmid("r2", r2[0].seq, False, []), cu.pcr3_fw("ACGTAC"), cu.PCR3_REV,
                 min_anneal=18)
    check(f"PCR3 gives the published {cu.LIBRARY_LEN_PUBLISHED} bp",
          r3[0].length, cu.LIBRARY_LEN_PUBLISHED)
    lib = r3[0].seq
    check("...and it contains the guide", SPACER in lib)
    check("...and the 6-nt RSL", "N" * cu.RSL_LEN in lib)
    check("the custom read primer's 3' end abuts the spacer exactly",
          lib.find(cu.CUSTOM_SEQ_PRIMER) + len(cu.CUSTOM_SEQ_PRIMER),
          lib.find(SPACER))

check.section("the final library, and the sequencing primers (lib/seqprimers.py)")
I5 = "ACGTAC"
drawn = cu.final_library(spacer=SPACER, i5=I5)
check("...so the page's segment names describe the real amplicon",
      [s.name for s in drawn][:3], ["Illumina P5", "i5 sample index", "TruSeq Read 1 site"])
# The one claim the drawing cannot make on its own: that it IS the simulated amplicon.
if lib is not None:
    check("the drawn library is the PCR3 product, base for base", drawn.top(), lib)
L = cu.final_library()
check("every declared sequencing primer lands where declared", sp.verify(L, cu.SEQ_PRIMERS), [])
r1p, i1p, i2p = cu.SEQ_PRIMERS[0], sp.TRUSEQ["I1"], cu.SEQ_PRIMERS[3]
check("Read 1 (CRIPSRSEQ) reads the spacer from cycle 1",
      sp.locate(L, r1p).reads_from, "sgRNA spacer")
check("...with no unpaired 5' base: PCR2-F templated the whole site",
      sp.locate(L, r1p).free5, 0)
check("the i7 read is the RSL -- the stock TruSeq Index-1 primer, on the vector's own site",
      sp.locate(L, i1p).reads_from, "RSL")
check("...and that site is the templated adapter, not a primer tail",
      sp.locate(L, i1p).covers, ("Illumina Read 2 / Index 1 site",))
check("the i5 read runs off the flow-cell P5 oligo and starts on the index",
      sp.locate(L, i2p).reads_from, "i5 sample index")
check("the RC-workflow i5 primer has NO site: PCR3-F omits the ACAC its 3' end needs",
      sp.locate(L, sp.TRUSEQ["I2"]), None)
check("...and the stock Read-2 primer has none either: the templated site is the 33-nt "
      "Index-1 oligo, one A short of Read 2's 34-nt complement",
      sp.locate(L, sp.TRUSEQ["R2"]) is None
      and revcomp(TRUSEQ_READ2) == "A" + cu.ILLUMINA_ADAPTER)
# The stock Read-1 primer is a near-miss, and the reason matters: its 3' end pairs, so it
# WOULD prime -- it just starts read 1 in U6, 23 nt short of the spacer.
# Measured on `drawn`, which the check above pins base-for-base to the simulated amplicon
# whenever the parent map is present -- so these three need no map of their own.
LIBSEQ = drawn.top()
check("the stock Read-1 primer's 3' 29 nt do match the library",
      il.TRUSEQ_READ1[4:] in LIBSEQ)
check("...and seqprimers locates it, reporting the 5' ACAC as a 4-nt flap",
      sp.locate(L, sp.TRUSEQ["R1"]).free5, 4)
check("...and priming there would start Read 1 in U6, not on the spacer",
      LIBSEQ.find(SPACER) - (LIBSEQ.find(il.TRUSEQ_READ1[4:]) + len(il.TRUSEQ_READ1[4:])),
      len(cu.U6_OVERLAP[-23:]))
check("all four roles are declared", sorted({q.role for q in cu.SEQ_PRIMERS}),
      sorted(sp.ROLES))

check.section("the published readout vs the primers in to_debug/")
# Martin's CRISPR_PCR1-F/-R, defined once in crisprumi.py and used by the page too.
USR_F, USR_R = cu.OUR_PCR1_FW, cu.OUR_PCR1_REV
check("CRISPR_PCR1-F is Schmierer's PCR1_FW with a 3-nt 5' extension",
      USR_F, "AAT" + cu.PCR1_FW)
check("CRISPR_PCR1-R's 5' tail starts with Schmierer's PCR1_REV 3' 15 nt",
      USR_R[:15], cu.PCR1_REV[-15:])
# On plain lentiGuide-Puro the same three PCRs run, but 39 bp shorter -- exactly the
# Illumina adapter + RSL that Schmierer inserted. That difference is the diagnostic.
if have(check, PARENT, label="the same readout run on plain lentiGuide-Puro"):
    gp = clone_guide(cu.parent_vector(), SPACER)
    g1 = amplify(gp, cu.PCR1_FW, cu.PCR1_REV, min_anneal=18)
    g2 = amplify(Plasmid("g1", g1[0].seq, False, []), cu.PCR2_FW, cu.PCR2_REV, min_anneal=18)
    g3 = amplify(Plasmid("g2", g2[0].seq, False, []), cu.pcr3_fw("ACGTAC"), cu.PCR3_REV,
                 min_anneal=18)
    check("the readout also runs on plain lentiGuide-Puro", len(g3), 1)
    check("...but 39 bp shorter -- the adapter + RSL Schmierer inserted",
          cu.LIBRARY_LEN_PUBLISHED - g3[0].length,
          len(cu.ILLUMINA_ADAPTER) + cu.RSL_LEN)
    check("...which is why the product size alone identifies the vector",
          g3[0].length, 249)

check.section("the drawings are placed by pairing, not by hand-typed columns")
import importlib.util                                               # noqa: E402

# by path: crispr-mip__10.1101+2024.03.28.587082/tools is also on sys.path and has a build_page of its own
_spec = importlib.util.spec_from_file_location("crisprumi_build_page", HERE / "build_page.py")
bp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bp)

r1 = bp.read1_scene()
pri, lib = r1.strands["CRIPSRSEQ"], r1.strands["library"]
check("the read primer is drawn as the paper's own sequence", pri.text(),
      cu.CUSTOM_SEQ_PRIMER)
check("its 3' base abuts the spacer: cycle 1 reads spacer base 1",
      pri.end(), lib.span("sgRNA spacer'")[0])
check("the primer's annealing block sits on the U6 3' end, column for column",
      pri.span("U6 anneal"), lib.span("U6 3' end'"))
# And a Scene refuses the off-by-one that a hand-typed indent would have produced.
try:
    bad = Scene()
    bad.strand("library", complement_segments(
        [Segment("U6 3' end", cu.U6_OVERLAP[-23:]),
         Segment("sgRNA spacer", "N" * cu.SPACER_LEN, placeholder=True)]), rev=True)
    bad.anneal("primer", [Segment("U6 anneal", cu.CUSTOM_SEQ_PRIMER[6:])], to="library",
               pair=("U6 anneal", "U6 3' end'"), shift=1)
    shifted = "drawn without complaint"
except ValueError:
    shifted = "refused"
check("a one-column shift of that primer is refused by the Scene, not drawn", shifted,
      "refused")

check.report()
