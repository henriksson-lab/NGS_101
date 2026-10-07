#!/usr/bin/env python3
"""
Self-test for pooled CRISPR screening.

This one earns its keep: it checks published primer and amplicon claims against the real
Addgene maps in ../ref/plasmids/, rather than restating them.

Those maps are third-party material and may not be in the clone (see
ref/plasmids/MANIFEST.md). Everything that reads them is therefore behind `have(...)`:
with the files present every check runs; without them the map-dependent groups SKIP, the
map-independent ones still run, and the suite still passes.

Run:  python3 lenticrispr-gecko-screen__10.1126+science.1247005/tools/selftest.py
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import crisprscreen as cs
import illumina as il
import seqprimers
from checks import Check, Source, have, run_common
from chemdraw import revcomp
from crispr import (FILLER_LEN, OVERHANG_DOWNSTREAM, OVERHANG_UPSTREAM, SCAFFOLD_FE_HEAD,
                    SCAFFOLD_V1, SCAFFOLD_V1_HEAD, SPACER_LEN, U6_3PRIME, U6_PLUS1,
                    clone_guide, guide_context, guide_oligos, scaffold_of)
from plasmid import BsmBI, amplify, bind, digest, find_both, read_genbank

PL = HERE.parents[0] / "ref" / "plasmids"
G_SPACER = cs.DEMO_SPACER_G             # starts with G: no +1 appended
N_SPACER = cs.DEMO_SPACER               # does not: +1 G appended

MANIFEST = "see ref/plasmids/MANIFEST.md for the exact URL"

# The five Addgene maps every computed claim on the page rests on. `origin` names the
# Addgene ID so a skip line says what to download, not just what is missing.
MAP_FILES = (
    ("addgene-49535_lentiCRISPRv1.gb", "lentiCRISPR v1",
     "Addgene #49535 (lentiCRISPR v1), depositor full sequence"),
    ("addgene-52961_lentiCRISPRv2.gb", "lentiCRISPRv2",
     "Addgene #52961 (lentiCRISPRv2), full sequence"),
    ("addgene-52963_lentiGuide-Puro.gb", "lentiGuide-Puro",
     "Addgene #52963 (lentiGuide-Puro), full sequence"),
    ("addgene-52962_lentiCas9-Blast.gb", "lentiCas9-Blast",
     "Addgene #52962 (lentiCas9-Blast), full sequence"),
    ("addgene-59702_pXPR_011.gb", "pXPR_011",
     "Addgene #59702 (pXPR_011), full sequence"),
)
MAPS = [Source(PL / f, f"{origin} -- {MANIFEST}") for f, _n, origin in MAP_FILES]

# Optional independent SnapGene maps, kept only as cross-checks: each pair is guarded on
# its own, so one missing map costs one comparison and nothing else.
ALT = PL / "alt-snapgene"
CROSSCHECKS = (
    # (Addgene map name in V, SnapGene file, published rotation offset, Addgene origin)
    ("lentiCRISPRv2", "snapgene-52961_lentiCRISPRv2.gb", 14256),
    ("lentiCas9-Blast", "snapgene-52962_lentiCas9-Blast.gb", 4315),
    ("lentiGuide-Puro", "snapgene-52963_lentiGuide-Puro.gb", 2593),
)
SG_ORIGIN = ("SnapGene's public CRISPR plasmid set, converted with Biopython -- "
             + MANIFEST)
WEISSMAN_AG = Source(
    PL / "addgene-66217_CRISPRi-library-backbone_pU6-sgRNA-EF1Alpha-puro-T2A-BFP.gb",
    "Addgene #66217 (CRISPRi v2 library backbone), depositor GenBank -- " + MANIFEST)
WEISSMAN_SG = Source(ALT / "snapgene-62217_pCRISPRia-v2.gb",
                     "SnapGene's CRISPRi_library_backbone (Addgene #62217) -- " + MANIFEST)

check = Check()
run_common(check)

# ===========================================================================
# Checks that need no plasmid map: everything that follows from lib/crispr.py
# and lib/crisprscreen.py alone. These run on a fresh clone.
# ===========================================================================
check.section("overhangs and guide oligos (no plasmid map needed)")
check("the downstream overhang is the scaffold's own first 4 nt",
      OVERHANG_DOWNSTREAM, SCAFFOLD_V1[:4])
check("revcomp(downstream overhang) is the published AAAC oligo overhang",
      revcomp(OVERHANG_DOWNSTREAM), "AAAC")
top, bot = guide_oligos(G_SPACER)
check("G-initiated spacer takes no extra G", top, OVERHANG_UPSTREAM + G_SPACER)
check("its bottom oligo starts AAAC", bot.startswith(revcomp(OVERHANG_DOWNSTREAM)))
check("non-G spacer gets the +1 G appended",
      guide_oligos(N_SPACER)[0], OVERHANG_UPSTREAM + U6_PLUS1 + N_SPACER)
check("...and its bottom oligo gets the matching C", guide_oligos(N_SPACER)[1].endswith("C"))
# geometry of the annealed pair, as drawn in step (3): the two 5' overhangs are at OPPOSITE
# ends, and the bottom oligo's 3'-terminal C is what pairs with the +1 G -- not, as the
# drawing used to imply, part of the AAAC overhang.
_nt, _nb = guide_oligos(N_SPACER)
check("bottom oligo is AAAC + revcomp(spacer) + C, in that order",
      _nb, revcomp(OVERHANG_DOWNSTREAM) + revcomp(N_SPACER) + revcomp(U6_PLUS1))
check("the paired core is the top oligo minus CACC against the bottom minus AAAC",
      revcomp(_nt[len(OVERHANG_UPSTREAM):]), _nb[len(OVERHANG_DOWNSTREAM):])
check("the bottom oligo's 3'-terminal base is the +1 G's partner, NOT part of the AAAC",
      revcomp(_nb[-1]), U6_PLUS1)
check("so the two 5' overhangs sit at opposite ends of the duplex",
      _nt.startswith(OVERHANG_UPSTREAM) and _nb.startswith(revcomp(OVERHANG_DOWNSTREAM))
      and len(_nt) - len(OVERHANG_UPSTREAM) == len(_nb) - len(OVERHANG_DOWNSTREAM))
check("the +1 G is NOT double-counted (U6_3PRIME excludes it)",
      U6_3PRIME.endswith("CACC") and not U6_3PRIME.endswith("CACCG"))
check("the two scaffolds are genuinely distinguishable",
      SCAFFOLD_V1_HEAD != SCAFFOLD_FE_HEAD and not SCAFFOLD_V1.startswith(SCAFFOLD_FE_HEAD))

check.section("no sequence is retyped: every block is composed from lib/")
# The Illumina blocks are built from lib/illumina.py rather than pasted. illumina.py holds
# the same 58-mer only under its NEBNext name, so pin the two together: if either moves,
# this fails rather than the page quietly disagreeing with lib.
check("P5_READ1 is P5 + TruSeq read 1, overlapping by ACAC (58 nt)",
      (cs.P5_READ1, len(cs.P5_READ1)), (il.TRUSEQ_P5_FULL, 58))
check("P7 and the read-2 block are lib/illumina.py's own constants",
      (cs.P7, cs.READ2_IDX), (il.P7, il.TRUSEQ_READ2))
check("the U6 region ends with lib's U6_3PRIME, so the +1 G is not baked in",
      cs.U6_REGION.endswith(U6_3PRIME) and not cs.U6_REGION.endswith(U6_3PRIME + U6_PLUS1))
check("the three published U6 sites are nested 3' ends of it",
      [len(s) for s in (cs.SITE_U6_SHORT, cs.SITE_U6_24, cs.SITE_U6_LONG)]
      == [22, 24, 33]
      and (cs.U6_REGION + U6_PLUS1).endswith(cs.SITE_U6_SHORT)
      and (cs.U6_REGION + U6_PLUS1).endswith(cs.SITE_U6_24))
check("Joung's reverse site is cut from lib's SCAFFOLD_V1, not retyped",
      revcomp(cs.SITE_SCAFFOLD) in SCAFFOLD_V1)
check("Shalem's reverse primer ends in the same cPPT site, written once",
      cs.SHALEM_R1.endswith(cs.SITE_CPPT[:15]) and len(cs.SHALEM_R1), 44)
check("the v2 adaptor's reverse primer is graft + annealing half, nothing else",
      cs.V2ADAPTOR_GRAFT + cs.V2ADAPTOR_ANNEAL, cs.V2ADAPTOR_R)

check.section("stagger ladders")
check("GPP uses 8 staggers and skips length 5",
      sorted(len(s) for s in cs.GPP_STAGGERS), [0, 1, 2, 3, 4, 6, 7, 8])
check("Joung uses 10 staggers of 9-18 nt",
      sorted(len(s) for s in cs.JOUNG_STAGGERS), list(range(9, 19)))
check("every stagger keeps the same 3' annealing site",
      all(cs.gpp_p5(s).endswith(cs.SITE_U6_SHORT) for s in cs.GPP_STAGGERS))
check("spacer is 20 nt", SPACER_LEN, 20)

# ===========================================================================
# Everything below reads the Addgene maps. They are third-party files the repo
# does not ship, so the whole block is behind one `have()`.
# ===========================================================================
if not all(s.path.exists() for s in MAPS):          # a heading for the skips to sit under
    check.section("checks that read the Addgene plasmid maps in ref/plasmids/")
if have(check, *MAPS):
    V = {n: read_genbank(str(PL / f)) for f, n, _o in MAP_FILES}
    CLONED = {n: clone_guide(p, G_SPACER) for n, p in V.items() if n in cs.VECTOR_READOUT}

    check.section("plasmid files are real and self-consistent")
    for n, p in V.items():
        check(f"{n}: parses, circular, plausible size", p.circular and 5000 < len(p) < 20000)
        check(f"{n}: sequence is pure ACGT(N)", set(p.seq) <= set("ACGTN"))

    check.section("BsmBI cloning site (computed from the maps, not quoted)")
    for n in ("lentiCRISPR v1", "lentiCRISPRv2", "lentiGuide-Puro"):
        cuts = digest(V[n], BsmBI)
        check(f"{n}: exactly 2 BsmBI sites", len(cuts), 2)
        check(f"{n}: upstream overhang {OVERHANG_UPSTREAM}", cuts[0].overhang, OVERHANG_UPSTREAM)
        check(f"{n}: downstream overhang {OVERHANG_DOWNSTREAM}", cuts[1].overhang, OVERHANG_DOWNSTREAM)
        check(f"{n}: filler is {FILLER_LEN} bp", cuts[1].top_cut - cuts[0].top_cut, FILLER_LEN)
    check("lentiCas9-Blast has NO BsmBI site -- Cas9 delivery only",
          len(digest(V["lentiCas9-Blast"], BsmBI)), 0)
    check("lentiCas9-Blast carries no sgRNA scaffold", scaffold_of(V["lentiCas9-Blast"]), None)

    check.section("cloning the guide into the real maps")
    for n in CLONED:
        for sp, lbl in ((G_SPACER, "G"), (N_SPACER, "non-G")):
            check(f"{n}: {lbl} spacer clone reconstitutes U6 -> spacer -> scaffold",
                  guide_context(clone_guide(V[n], sp), sp))

    check.section("scaffold identity")
    for n in ("lentiCRISPR v1", "lentiCRISPRv2", "lentiGuide-Puro"):
        check(f"{n} carries the original v1 scaffold", scaffold_of(V[n]), "v1")
    check("pXPR_011 carries the F+E scaffold instead", scaffold_of(V["pXPR_011"]), "F+E")
    check("there is NO v1-vs-v2 scaffold difference -- both are the original",
          scaffold_of(V["lentiCRISPR v1"]) == scaffold_of(V["lentiCRISPRv2"]))

    check.section("the U6 primer sites really are in the maps")
    # The U6 sites are nested 3' ends of one stretch, and that stretch is really in the maps.
    for n in ("lentiCRISPR v1", "lentiCRISPRv2", "lentiGuide-Puro"):
        check(f"{n}: the U6 stretch the primer sites are cut from is in the map, once",
              len(find_both(V[n].seq, cs.U6_REGION, True)), 1)

    check.section("readout primers: computed amplicons vs published sizes")
    gpp_f = cs.gpp_p5("")
    kermit = cs.gpp_p7(cs.GPP_P7_INDEX_A01, cs.SITE_CPPT)
    beaker = cs.gpp_p7(cs.GPP_P7_INDEX_A01, cs.SITE_EFS)
    jf, jr = cs.joung_fwd(cs.JOUNG_STAGGERS[0]), cs.joung_rev(cs.JOUNG_REV_INDEX_1)

    def size(vec, f, r):
        a = amplify(CLONED[vec], f, r)
        return a[0].length if a else None

    check("GPP BEAKER on lentiCRISPRv2 == the published 285 nt",
          size("lentiCRISPRv2", gpp_f, beaker), 285)
    check("Joung's pair is vector-independent (one size across all three backbones)",
          len({size(v, jf, jr) for v in CLONED}), 1)
    check("Joung's amplicon matches the published ~260 bp", size("lentiCRISPRv2", jf, jr), 259)
    check("GPP KERMIT gives a short product on lentiGuide-Puro",
          size("lentiGuide-Puro", gpp_f, kermit) < cs.MAX_USABLE_AMPLICON)
    check("GPP KERMIT gives a short product on lentiCRISPR v1",
          size("lentiCRISPR v1", gpp_f, kermit) < cs.MAX_USABLE_AMPLICON)
    check("...but NOT on lentiCRISPRv2: the cPPT sits upstream of U6 there, so the product "
          "would run the long way round the plasmid",
          size("lentiCRISPRv2", gpp_f, kermit) > cs.MAX_USABLE_AMPLICON)
    check("Shalem's step-1 pair fails on lentiCRISPRv2 for the same reason",
          size("lentiCRISPRv2", cs.SHALEM_F1, cs.SHALEM_R1) > cs.MAX_USABLE_AMPLICON)
    check("...while working on lentiCRISPR v1",
          size("lentiCRISPR v1", cs.SHALEM_F1, cs.SHALEM_R1) < cs.MAX_USABLE_AMPLICON)
    check("BEAKER gives no usable product on lentiGuide-Puro either",
          size("lentiGuide-Puro", gpp_f, beaker) > 500)
    check("the +1 G shifts every amplicon by exactly 1 nt",
          amplify(clone_guide(V["lentiCRISPRv2"], N_SPACER), gpp_f, beaker)[0].length
          - size("lentiCRISPRv2", gpp_f, beaker), 1)

    check.section("where each primer actually anneals")
    check("GPP P5_ARGON binds the U6 promoter of every guide vector",
          all(bind(CLONED[v], gpp_f) for v in CLONED))
    check("Joung's reverse binds the SCAFFOLD, which is why it is vector-independent",
          all(find_both(CLONED[v].seq, cs.SITE_SCAFFOLD, True) for v in CLONED))
    check("BEAKER's site IS present in lentiGuide-Puro -- it has an EF-1a cassette too, so "
          "the vector/primer rule is about product LENGTH, not about the site being missing",
          bool(find_both(V["lentiGuide-Puro"].seq, cs.SITE_EFS, True)))
    check("the v2 adaptor grafts the cPPT site onto its 5' tail",
          cs.V2ADAPTOR_R.startswith(cs.SITE_CPPT))

    check.section("sizing step-1 PCR: cloned library vs leftover stuffer (ira1.md S10)")
    # The single most informative gel measurement: step 1 of the two-step readout separates a
    # correctly cloned library from uncut/re-ligated backbone by ~1.9 kb, because the filler
    # sits inside the amplicon. Both templates amplify, so only the SIZE tells them apart.
    def stuffer_size(vec, f, r):
        a = amplify(V[vec], f, r)
        return a[0].length if a else None

    check("v2 adaptor on cloned lentiCRISPRv2 is 287 bp",
          size("lentiCRISPRv2", cs.V2ADAPTOR_F, cs.V2ADAPTOR_R), 287)
    check("v2 adaptor on stuffer-bearing lentiCRISPRv2 is 2148 bp",
          stuffer_size("lentiCRISPRv2", cs.V2ADAPTOR_F, cs.V2ADAPTOR_R), 2148)
    # Both primer sites lie outside the excised filler, so the amplicon shrinks by exactly
    # as much as the plasmid does -- that is the invariant, not any FILLER_LEN arithmetic
    # (FILLER_LEN counts between the BsmBI nicks; cloning retains the 4-nt overhangs).
    for _v in ("lentiCRISPRv2", "lentiGuide-Puro"):
        check(f"{_v}: the step-1 amplicon shrinks by exactly as much as the plasmid does",
              stuffer_size(_v, cs.V2ADAPTOR_F, cs.V2ADAPTOR_R)
              - size(_v, cs.V2ADAPTOR_F, cs.V2ADAPTOR_R),
              len(V[_v]) - len(CLONED[_v]))
    check("that shrinkage is ~1.9 kb, so a gel separates the two cases trivially",
          1800 < len(V["lentiCRISPRv2"]) - len(CLONED["lentiCRISPRv2"]) < 1900)
    check("v2 adaptor on cloned lentiGuide-Puro is 556 bp",
          size("lentiGuide-Puro", cs.V2ADAPTOR_F, cs.V2ADAPTOR_R), 556)
    check("v2 adaptor on stuffer-bearing lentiGuide-Puro is 2417 bp",
          stuffer_size("lentiGuide-Puro", cs.V2ADAPTOR_F, cs.V2ADAPTOR_R), 2417)
    check("Shalem F1/R1 on cloned lentiGuide-Puro is 316 bp",
          size("lentiGuide-Puro", cs.SHALEM_F1, cs.SHALEM_R1), 316)
    check("Shalem F1/R1 on stuffer-bearing lentiGuide-Puro is 2177 bp",
          stuffer_size("lentiGuide-Puro", cs.SHALEM_F1, cs.SHALEM_R1), 2177)
    # The observed ~1 kb band matches none of the four sizes the two vectors can give.
    _SIZES = {287, 556, 2148, 2417}
    check("no cloned-or-stuffer combination of these two vectors gives a ~1 kb step-1 product",
          not any(850 <= x <= 1200 for x in _SIZES))
    check("the two step-1 pairs share their forward primer, so only the REVERSE "
          "distinguishes them -- the easy mix-up to check first",
          cs.V2ADAPTOR_F == cs.SHALEM_F1 and cs.V2ADAPTOR_R != cs.SHALEM_R1)

    check.section("the two-step (v2 adaptor) readout, as drawn on the page")
    import build_page as bp  # noqa: E402  (reads the same maps; guarded by have() above)

    check("the reverse primer splits into a 25-nt graft and a 20-nt annealing half",
          (len(bp.GRAFT), len(bp.ANNEAL)), (25, 20))
    check("only the 3' half is present in the template",
          bool(find_both(CLONED["lentiCRISPRv2"].seq, bp.ANNEAL, True)))
    check("the grafted 5' half is NOT downstream of the guide, which is why it is grafted",
          size("lentiCRISPRv2", cs.gpp_p5(""), cs.gpp_p7(cs.GPP_P7_INDEX_A01, cs.SITE_CPPT))
          > cs.MAX_USABLE_AMPLICON)
    for _v, _r1 in (("lentiCRISPRv2", 288), ("lentiGuide-Puro", 557)):
        check(f"round 1 on {_v} is {_r1} bp", len(bp.step1_construct(_v)), _r1)
        check(f"round 2 on {_v} is 351 bp", len(bp.step2_construct(_v)), 351)
    check("round 2 is the SAME size on both vectors, so it cannot identify the backbone",
          len(bp.step2_construct("lentiCRISPRv2")),
          len(bp.step2_construct("lentiGuide-Puro")))
    check("...while round 1 differs by 269 bp, so it can",
          len(bp.step1_construct("lentiGuide-Puro")) - len(bp.step1_construct("lentiCRISPRv2")),
          269)
    check("the round-1 product ends in the grafted cPPT, reverse-complemented",
          bp.step1_construct("lentiCRISPRv2").top().endswith(revcomp(bp.GRAFT)))
    check("that grafted tail is absent from the plasmid at that position -- it is new sequence",
          revcomp(bp.GRAFT) not in
          amplify(CLONED["lentiCRISPRv2"], cs.SHALEM_F1, bp.ANNEAL)[0].seq)
    check("round 1 carries no Illumina P5", cs.P5_READ1 not in bp.step1_construct("lentiCRISPRv2").top())
    check("round 2 does carry P5", bp.step2_construct("lentiCRISPRv2").top().startswith(cs.P5_READ1[:29]))
    check("round 2 carries the sgRNA spacer placeholder",
          "N" * cs.SPACER_LEN if False else
          "N" * 20 in bp.step2_construct("lentiCRISPRv2").top())
    # the vector-specific stretch is what differs, and round 2 removes it
    _v1 = [x for x in bp.step1_construct("lentiCRISPRv2").segments if x.name == "vector-specific"]
    _v2 = [x for x in bp.step1_construct("lentiGuide-Puro").segments if x.name == "vector-specific"]
    check("the vector-specific stretch is 84 nt on lentiCRISPRv2", len(_v1[0].top), 84)
    check("...and 353 nt on lentiGuide-Puro", len(_v2[0].top), 353)
    check("their difference is exactly the round-1 size difference",
          len(_v2[0].top) - len(_v1[0].top), 269)

    check.section("sequencing primers (lib/seqprimers.py locates each one)")
    for _title, _con in bp.AMPS:
        check(f"every declared sequencing primer lands where declared -- {_title}",
              seqprimers.verify(_con, cs.SEQ_PRIMERS), [])
    for _v in ("lentiCRISPRv2", "lentiGuide-Puro"):
        check(f"...and on the two-step round-2 library ({_v})",
              seqprimers.verify(bp.step2_construct(_v), cs.SEQ_PRIMERS), [])
    check("read 1 reports the stagger and U6, never the spacer first",
          all(seqprimers.locate(c, cs.SEQ_PRIMERS[0]).reads[0] != "N" for _t, c in bp.AMPS))
    check("the i7 read reports the reverse complement of the index in the primer",
          seqprimers.locate(bp.AMPS[0][1], cs.SEQ_PRIMERS[1]).reads[:8],
          revcomp(cs.GPP_P7_INDEX_A01))
    check("the page's array-oligo flanks are the vector's own U6 and scaffold",
          bp.ARRAY_5.endswith(cs.SITE_U6_LONG + U6_PLUS1)
          and SCAFFOLD_V1.startswith(bp.ARRAY_3))

    # ----------------------------------------------------------------- cross-checks
    # SnapGene builds its maps from the depositors' data rather than copying Addgene, so
    # agreement is meaningful -- and these are the numbers MANIFEST.md quotes. Optional:
    # each pair is guarded separately, so a missing alt map skips one comparison only.
    check.section("independent SnapGene maps agree with the Addgene ones (optional)")
    for _name, _file, _offset in CROSSCHECKS:
        _src = Source(ALT / _file, SG_ORIGIN)
        if have(check, _src, label=f"{_name}: SnapGene cross-check"):
            _sg = read_genbank(str(_src.path))
            _ag = V[_name].seq
            check(f"{_name}: SnapGene map is the same length", len(_sg), len(_ag))
            check(f"{_name}: ...and an exact circular rotation at the published offset "
                  f"{_offset}", (_ag + _ag).find(_sg.seq) + 1, _offset)

check.section("the Weissman CRISPRi backbone pair (optional cross-check)")
if have(check, WEISSMAN_AG, WEISSMAN_SG,
        label="Weissman CRISPRi backbone cross-check"):
    # MANIFEST.md's claim about the Weissman pair: base-identical apart from the N20 run
    # the depositor left at the guide position.
    _ag = read_genbank(str(WEISSMAN_AG.path)).seq
    _sg = read_genbank(str(WEISSMAN_SG.path)).seq
    check("#66217 and SnapGene's #62217 map are the same length", len(_ag), len(_sg))
    check("#66217 leaves the 20-nt guide position as N", _ag.count("N"), 20)
    check("...and that N run is the ONLY difference between the two files",
          [i for i, (a, b) in enumerate(zip(_ag, _sg)) if a != b],
          [i for i, a in enumerate(_ag) if a == "N"])

check.report()
