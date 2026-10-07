#!/usr/bin/env python3
"""
Design a padlock ligation arm for the Schmierer vector (pLenti-Puro-AU-flip-3xBsmBI).

The extension arm is kept -- it sits in U6, which is unchanged. Only the ligation arm has
to move, and it has to dodge three things at once:

  1. the AU-flip substitutions, at scaffold positions 5 and 26;
  2. the vector's BUILT-IN Illumina Read-2 adapter -- the CRISPR-MIP backbone already
     carries that same sequence, so a capture spanning it puts two copies of one primer
     site into a single amplicon;
  3. the usual padlock requirement that the ligation arm be the hotter of the two.

Run: python3 crispr-umi-schmierer__10.15252+msb.20177834/tools/design_padlock.py
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))
sys.path.insert(0, str(HERE.parents[1] / "crispr-mip__10.1101+2024.03.28.587082" / "tools"))

import crisprmip as cm
import crisprumi as cu
from chemdraw import tm
from padlock import Padlock, capture
from plasmid import Plasmid, find_both

ARM_LEN = 23
SPACER = cu.DEMO_SPACER


def schmierer_vector(cloned: bool = True) -> Plasmid:
    """pLenti-Puro-AU-flip-3xBsmBI, rebuilt from the real lentiGuide-Puro map.

    One line, because the rebuild itself lives in crisprumi.py and is shared with the page
    and the selftest.
    """
    return cu.rebuild_vector(cloned=cloned, spacer=SPACER)


def candidates(v: Plasmid, lo: int, hi: int, tm_lo=56.0, tm_hi=60.5):
    """Ligation arms whose start lies `lo`..`hi` nt past the extension arm."""
    ext_end = (find_both(v.seq, cm.EXT_ARM, True)[0][0] + len(cm.EXT_ARM)) % len(v)
    ext_tm = tm(cm.EXT_ARM)
    out = []
    for gap in range(lo, hi):
        start = (ext_end + gap) % len(v)
        s = v.sub(start, start + ARM_LEN)
        if len(s) != ARM_LEN or set(s) - set("ACGT"):
            continue                              # skip the N's of the RSL
        t = tm(s)
        if not (tm_lo <= t <= tm_hi and t > ext_tm):
            continue
        if len(find_both(v.seq, s, True)) != 1:
            continue                              # must be unique in the vector
        out.append((gap, t, s))
    return out


def report(title, cands, v, note=""):
    print(f"\n{title}\n" + "=" * len(title))
    if note:
        print(f"  {note}")
    if not cands:
        print("  none found")
        return
    for gap, t, s in cands[:4]:
        pr = Padlock("x", ext_arm=cm.EXT_ARM, lig_arm=s, backbone=cm.backbone())
        caps = capture(v, pr)
        ok = len(caps) == 1
        print(f"  {s}  Tm {t:4.1f}  gap {gap:>3d} nt  "
              f"captures: {'1 site, circle %d nt' % len(caps[0].circle) if ok else 'FAILED'}")
        if ok:
            f = caps[0].fill
            print(f"      spacer in capture: {SPACER in f} | "
                  f"RSL region in capture: {'N' * cu.RSL_LEN in f} | "
                  f"vector adapter in capture: {cu.ILLUMINA_ADAPTER in f}")


def main() -> None:
    # Every arm below is scored against the real vector sequence. Without the parent map
    # there is nothing to score, so refuse in one line rather than print a shorter report.
    try:
        v = schmierer_vector()
    except cu.MissingMap as e:
        sys.exit(f"cannot design a ligation arm: {e}")
    print(f"target: {v.name}, cloned, {len(v)} bp")
    print(f"extension arm (kept): {cm.EXT_ARM}  Tm {tm(cm.EXT_ARM):.1f}")
    print(f"\nregion map, measured from the end of the extension arm:")
    for lbl, sub in (("+1 G + spacer", SPACER), ("AU-flip scaffold", cu.SCAFFOLD_AU_FLIP[:24]),
                     ("terminator", cu.TERMINATOR), ("vector Illumina adapter", cu.ILLUMINA_ADAPTER),
                     ("RSL", "N" * cu.RSL_LEN), ("downstream vector", cu.DOWNSTREAM[:20])):
        h = find_both(v.seq, sub, True)
        ext_end = (find_both(v.seq, cm.EXT_ARM, True)[0][0] + len(cm.EXT_ARM)) % len(v)
        print(f"   {lbl:24s} starts at gap {(h[0][0]-ext_end) % len(v):>4d}")

    report("OPTION B  -- ligation arm BEFORE the vector's adapter (recommended)",
           candidates(v, 55, 104 - ARM_LEN + 1), v,
           "captures the guide only; no sequence is duplicated")
    report("OPTION A  -- ligation arm AFTER the RSL",
           candidates(v, 143, 210), v,
           "also captures the RSL -- but see the warning below")

    print("\nWARNING for option A")
    print("  The CRISPR-MIP backbone already contains", cm.READ2_SITE)
    print("  The vector contains                     ", cu.ILLUMINA_ADAPTER)
    print("  The second is a substring of the first, so an amplicon that spans both would")
    print("  carry two copies of the Read-2 primer site. To use option A the probe backbone")
    print("  must DROP its own Read-2 site and let the vector's serve instead -- in which")
    print("  case the RSL lands in the index position and is read as the i7.")


if __name__ == "__main__":
    main()
