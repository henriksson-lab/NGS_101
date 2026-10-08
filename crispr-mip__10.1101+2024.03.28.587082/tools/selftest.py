#!/usr/bin/env python3
"""
Self-test for CRISPR-MIP.

Checks the published probe design against the real lentiCRISPRv2 map that the screens
actually used (Brunello kinome library, Addgene #75314, backbone #52961).

Run:  python3 crispr-mip__10.1101+2024.03.28.587082/tools/selftest.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import crisprmip as cm
from checks import Check, Source, have, run_common
from chemdraw import revcomp
from crispr import SCAFFOLD_V1_HEAD, U6_3PRIME, clone_guide
import illumina as il
import seqprimers
from plasmid import amplify, find_both
from padlock import Padlock, capture, survives_exonuclease
from plasmid import Plasmid, read_genbank

VEC = Source(HERE.parents[1] / "lenticrispr-v1-screening__10.1126+science.1247005" /
             "ref" / "plasmids" / "addgene-52961_lentiCRISPRv2.gb",
             "Addgene #52961 (lentiCRISPRv2, Zhang lab); that directory's "
             "ref/plasmids/MANIFEST.md records exactly how the map was obtained")
XLSX = Source(HERE.parent / "ref" / "TableS2_primers_and_probes.xlsx",
              "Supplementary Table S2 (media-2.xlsx) of doi:10.1101/2024.03.28.587082")
N_SPACER = cm.EXAMPLE_SPACER             # Brunello spacers are not G-initiated
G_SPACER = "GTCGCTGAGTACTTCGAAAT"

check = Check()
run_common(check)

# Third-party source material is not committed (see ref/MANIFEST.md): the checks that read
# it are guarded, so a fresh clone skips them loudly instead of failing as if the chemistry
# were wrong.
HAVE_VEC = have(check, VEC, label="everything measured on the real lentiCRISPRv2 map: "
                                  "capture site, circle, inverse PCR, restriction sites, "
                                  "and the page's drawings")

probe = cm.probe()
model_cap = cm.modeled_capture(N_SPACER)

check.section("probe geometry, from Table S2")
check("probe is 126 nt as actually listed in Table S2", len(probe), 126)
check("...which DISAGREES with the Methods text's '134 bp' by 8 nt",
      cm.PROBE_LEN_IN_TEXT - cm.PROBE_LEN, 8)
check("extension arm is 19 nt", len(cm.EXT_ARM), 19)
check("ligation arm is 23 nt", len(cm.LIG_ARM), 23)
check("backbone is 84 nt", cm.BACKBONE_LEN, 84)
check("backbone decomposes exactly: UMI + read2 site + i7 + read1 site",
      cm.UMI_LEN + len(cm.READ2_SITE) + cm.INDEX_LEN + len(cm.READ1_SITE), cm.BACKBONE_LEN)
check("nine probes, differing ONLY in an 8-nt i7 index", len(cm.PROBE_INDICES), 9)
check("all indices are 8 nt and distinct",
      all(len(i) == 8 for i in cm.PROBE_INDICES) and len(set(cm.PROBE_INDICES)) == 9)

# The sequences above are assembled from lib/ constants (TruSeq sites, P5/P7, scaffold);
# this is what keeps that assembly honest: compare with Table S2 itself.
_s2 = None
S2 = "Table S2 cross-check: every probe and primer, verbatim"
if have(check, XLSX, label=S2):
    try:
        import openpyxl
    except ImportError:
        check.skip(S2, "openpyxl is not installed, so the .xlsx cannot be read (the file "
                       "itself is present) -- install it with: pip3 install openpyxl")
    else:
        _ws = openpyxl.load_workbook(XLSX.path).active
        _s2 = {r[0].strip(): r[1].strip()
               for r in _ws.iter_rows(min_row=2, values_only=True) if r[0]}
if _s2:
    check("all nine probes, as assembled, equal Table S2 verbatim (after the 5' phosphate)",
          [cm.probe(i).sequence() for i in cm.PROBE_INDICES],
          [_s2[f"MIP_probe_{k + 1}"].split("/")[-1] for k in range(9)])
    check("P5_tracrRNA_fwd, as assembled, equals Table S2", cm.P5_TRACR_FWD, _s2["P5_tracrRNA_fwd"])
    check("P7_tracrRNA_rev, as assembled, equals Table S2", cm.P7_TRACR_REV, _s2["P7_tracrRNA_rev"])
    check("Table S2's per-probe variable stretch is 8 nt, not 6 (bases 7-8 differ too)",
          len({_s2[f"MIP_probe_{k + 1}"][-len(cm.EXT_ARM) - len(cm.READ1_SITE) - 2:
                  -len(cm.EXT_ARM) - len(cm.READ1_SITE)] for k in range(9)}) > 1)
check("probes differ only at the index", len({cm.backbone(i).replace(i, "", 1)
                                              for i in cm.PROBE_INDICES}), 1)
check("P5/P7 are NOT in the probe", il.P5[:15] not in probe.sequence()
      and il.P7[:15] not in probe.sequence())
e_tm, l_tm = probe.arm_tms()
check("extension arm Tm is ~54 C (published 53.8)", round(e_tm), 54)
check("ligation arm Tm is ~58 C (published 58.3)", round(l_tm), 58)
check("the LIGATION arm is the hotter one -- otherwise the polymerase displaces it "
      "and the circle never closes", probe.arm_tm_gap() > 0)

check.section("source-independent capture model")
check("modeled gap-fill is 112 nt for a non-G spacer", model_cap.gap, cm.CAPTURED_SPAN)
check("modeled fill contains the complete spacer", N_SPACER in model_cap.fill)
check("modeled fill is +1 G + spacer + canonical scaffold + verified flank",
      model_cap.fill,
      "G" + N_SPACER + cm.SCAFFOLD_V1 + cm.CAPTURE_3PRIME_FLANK)
check("a G-initiated spacer gives one base less",
      cm.modeled_capture(G_SPACER).gap, cm.CAPTURED_SPAN - 1)
model_amp = amplify(model_cap.circle, cm.P5_TRACR_FWD, cm.P7_TRACR_REV, min_anneal=15)
check("modeled circle gives exactly one inverse-PCR product", len(model_amp), 1)
check("modeled final library is 269 bp", model_amp[0].length, 269)

check.section("where the arms land on the real vector")
check("the extension arm ends exactly at the U6 +1 position",
      cm.EXT_ARM.endswith(U6_3PRIME))
if HAVE_VEC:
    p = read_genbank(str(VEC))
    v_n, v_g = clone_guide(p, N_SPACER), clone_guide(p, G_SPACER)
    caps_n = capture(v_n, probe)
    caps_g = capture(v_g, probe)
    check("exactly one capture site in the whole 13 kb vector", len(caps_n), 1)
    check("the probe is not promiscuous on the empty backbone either",
          len(capture(p, probe)) <= 1)
    c = caps_n[0]
    check("source map and source-independent model give the same gap-fill",
          c.fill, model_cap.fill)
    check("gap-fill is the published 112 nt (non-G spacer)", c.gap, cm.CAPTURED_SPAN)
    check("a G-initiated spacer gives one base less", caps_g[0].gap, cm.CAPTURED_SPAN - 1)
    check("the captured sequence contains the complete sgRNA spacer", N_SPACER in c.fill)
    check("the fill starts at the Pol III +1 base", c.fill[0], "G")
    check("the fill runs through the scaffold", SCAFFOLD_V1_HEAD in c.fill)
    check("circle = probe + gap", len(c.circle), len(probe) + c.gap)
    check("the circle is closed", c.circle.circular)

check.section("the circle is the selection")
if HAVE_VEC:
    check("a closed circle survives exonuclease I/III", survives_exonuclease(c.circle))
check("an un-ligated linear probe does not",
      not survives_exonuclease(Plasmid("probe", probe.sequence(), False, [])))
if HAVE_VEC:
    check("linear genomic DNA does not",
          not survives_exonuclease(Plasmid("gDNA", v_n.seq, False, [])))
check("four independent selection steps, none of them a size selection",
      len(cm.SELECTION_STEPS), 4)

if HAVE_VEC:
    check.section("probe arms must be ordered and co-strand")
    backwards = Padlock("swapped", ext_arm=cm.LIG_ARM, lig_arm=cm.EXT_ARM,
                        backbone=cm.backbone())
    check("swapping the arms captures NOTHING at a plausible gap -- the probe would have to "
          "span the long way round the plasmid", capture(v_n, backwards), [])
    check("...and that span is indeed most of the vector",
          capture(v_n, backwards, max_gap=10**9)[0].gap > 10000)
    nonsense = Padlock("absent", ext_arm="ACGTACGTACGTACGTACG",
                       lig_arm=cm.LIG_ARM, backbone=cm.backbone())
    check("a probe whose arm is not present captures nothing", capture(v_n, nonsense), [])
    check("arms are searched on both strands",
          len(capture(Plasmid("rc", revcomp(v_n.seq), True, []), probe)), 1)

if HAVE_VEC:
    check.section("amplification off the circle")
    check("the P5 primer's 3' tail lies INSIDE the captured sequence",
          bool(find_both(c.fill, cm.P5_TRACR_TAIL, False)))
    check("the P7 primer's 3' tail lies INSIDE the captured sequence",
          bool(find_both(c.fill, cm.P7_TRACR_TAIL, False)))
    check("...so neither primer can act on an unextended probe",
          cm.P5_TRACR_TAIL not in probe.sequence()
          and cm.P7_TRACR_TAIL not in probe.sequence())
    amp = amplify(c.circle, cm.P5_TRACR_FWD, cm.P7_TRACR_REV, min_anneal=15)
    check("the circle gives exactly one product", len(amp), 1)
    check("it is inverse PCR -- the product wraps the circle's origin", amp[0].wraps_origin)
    check("final library is 269 bp", amp[0].length, 269)
    check("the library carries the UMI", cm.UMI in amp[0].seq)
    check("the library carries the sgRNA spacer", N_SPACER in amp[0].seq)
    check("the library starts at P5 and ends at revcomp(P7)",
          amp[0].seq.startswith(il.P5) and amp[0].seq.endswith(revcomp(il.P7)))

check.section("sequencing")
check("read 2 is long enough for the 13-nt UMI", cm.READ2_CYCLES >= cm.UMI_LEN)
check("read 1 (60 cycles) covers the 20-nt spacer within the 112-nt capture",
      cm.READ1_CYCLES >= 20)

if HAVE_VEC:
    check.section("plasmid substrate: why MIP failed on supercoiled plasmid (ira1.md S9)")
    # The probe footprint = ext arm + captured region + lig arm. In the cloned #73179 map
    # it spans the GenBank origin, so it must be measured on a rotated sequence.
    _e = find_both(v_n.seq, cm.EXT_ARM, True)[0][0]
    _l = find_both(v_n.seq, cm.LIG_ARM, True)[0][0]
    _rot = max(_e, _l) if abs(_e - _l) > 1000 else min(_e, _l)
    _seq = v_n.seq[_rot:] + v_n.seq[:_rot]
    _e2 = find_both(_seq, cm.EXT_ARM, True)[0][0]
    _l2 = find_both(_seq, cm.LIG_ARM, True)[0][0]
    _lo, _hi = min(_e2, _l2), max(_e2, _l2) + len(cm.LIG_ARM)
    check("the probe footprint wraps the map origin, so it needs rotating to measure",
          abs(_e - _l) > 1000)
    check("rotated, the footprint is contiguous: ext arm + 112 nt + lig arm = 154 bp",
          _hi - _lo, len(cm.EXT_ARM) + 112 + len(cm.LIG_ARM))

    def _cuts(site, seq=_seq):
        """Cut positions on a circular sequence, both strands."""
        hits = set()
        for s in {site, revcomp(site)}:
            hits |= {m.start() for m in re.finditer(f"(?={s})", seq + seq[:len(site)])}
        return sorted(x for x in hits if x < len(seq))

    def _in_footprint(site):
        return [x for x in _cuts(site) if _lo - len(site) + 1 <= x <= _hi - 1]

    # Linearising the plasmid is the fix; the enzyme must not cut inside the footprint.
    for _nm, _site in (("NotI", "GCGGCCGC"), ("PacI", "TTAATTAA"), ("BamHI", "GGATCC")):
        check(f"{_nm} linearises the cloned vector with a single cut", len(_cuts(_site)), 1)
        check(f"...and {_nm} does not cut inside the probe footprint", _in_footprint(_site), [])
    # Both of the preprint's gDNA enzymes also spare the footprint.
    check("HindIII (the preprint's gDNA enzyme) spares the footprint", _in_footprint("AAGCTT"), [])
    # These two do cut the footprint in the vector itself -- they must stay rejected.
    check("EcoRI cuts inside the footprint, so it is unusable", bool(_in_footprint("GAATTC")))
    check("NheI cuts inside the footprint, so it is unusable", bool(_in_footprint("GCTAGC")))

if HAVE_VEC:
    check.section("drawings: every duplex is a chemdraw.Scene (columns computed, pairing checked)")
    import importlib.util
    _spec = importlib.util.spec_from_file_location("crisprmip_build_page", HERE / "build_page.py")
    bp = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(bp)
    from chemdraw import Scene, Segment

    _flat = "".join(x.top for x in bp.pcr_circle())
    _circ = c.circle.seq.replace(N_SPACER, "N" * 20)
    check("the PCR panel draws the whole circle exactly once (old row counted GTGC twice and "
          "dropped the 4 nt CGCT before the ligation arm)",
          len(_flat) == len(c.circle) and _flat in _circ + _circ)
    check("the P5 primer's tail ends 1 nt before the ligation arm -- wholly inside the fill "
          "(old drawing put its 3' CGC over the ligation arm's AGC)",
          c.fill.find(cm.P5_TRACR_TAIL) + len(cm.P5_TRACR_TAIL), len(c.fill) - 1)
    _c1, _c2 = bp.scene_pcr_cycle1(), bp.scene_pcr_cycle2()
    check("cycle 1: the P7 primer's tail sits exactly on the P7 site",
          _c1.strands["P7 primer"].span("P7 tail"), _c1.strands["circle"].span("P7 site"))
    check("cycle 2: the P5 primer's tail sits exactly on the copy's P5 site",
          _c2.strands["P5 primer"].span("P5 tail"), _c2.strands["first strand"].span("P5 site'"))
    check.raises("the P5 primer cannot anneal to the circle itself (same sense)",
                 lambda: (lambda sc: (sc.strand("c", bp.pcr_circle()),
                              sc.anneal("p5", [Segment("t", cm.P5_TRACR_TAIL)], to="c",
                                        pair=("t", "P5 site"))))(Scene()))

    _u, _f = bp.scene_capture(False)[0], bp.scene_capture(True)[0]
    _t = _u.strands["template"]
    check("unfilled: the extension arm's 3' end abuts the first gap base",
          _u.strands["ext arm"].end(), _t.span("+1 site'")[0])
    check("unfilled: the ligation arm's 5' end abuts the last gap base",
          _u.strands["lig arm"].col, _t.span("gap end site'")[1])
    check("filled: the new strand runs continuously from ext arm to lig arm",
          _f.strands["probe"].text(),
          cm.EXT_ARM + "".join(x.top for x in bp.gap_segments()) + cm.LIG_ARM)

    _r1, _i1, _r2 = bp.scene_read1(), bp.scene_index1(), bp.scene_read2()
    check("Read 1 primer's 3' end is immediately before the extension arm",
          _r1.strands["Read 1 primer"].end(), _r1.strands["top"].span("ext arm")[0])
    check("Index 1 primer's 3' end is immediately before the i7",
          _i1.strands["Index 1 primer"].end(), _i1.strands["top"].span("i7")[0])
    check("Read 2 primer's 3' end is immediately after the UMI (reads leftward into it)",
          _r2.strands["Read 2 primer"].col, _r2.strands["top"].span("UMI")[1])
    for _role, _sc in (("Read 1", _r1), ("Index 1 (i7)", _i1), ("Read 2", _r2)):
        _st = _sc.strands[_role.split(" (")[0] + " primer"]
        check(f"the {_role} drawing is the declared primer (same lib constant), 5'->3'",
              "".join(x.top for x in _st.segs), bp.primer(_role).seq)
    check("every declared sequencing primer lands where declared (lib/seqprimers.py)",
          seqprimers.verify(bp.LIB, cm.SEQ_PRIMERS), [])
    check.raises("the Read 2 site's own sense cannot be drawn annealed to the top strand "
                 "(old drawing: Read 2 primer written backwards beside the wrong strand)",
                 lambda: (lambda sc: (sc.strand("top", list(bp.LIB)),
                                      sc.anneal("r2", [Segment("r2", cm.READ2_SITE)], to="top",
                                                pair=("r2", "Read 2 site"))))(Scene()))

check.report()
