#!/usr/bin/env python3
"""Extract guide read counts from the conventional PCR readout (no UMI).

The guide is found by ANCHOR, not by offset: it is the 20 bp immediately before the
scaffold start `GTTTTAGAGC`. A fixed offset would be wrong, because the Broad GPP P5
primers deliberately carry staggers of 0-8 nt to maintain flowcell base diversity, so
the guide does not sit at a constant position. (The padlock readout has no stagger, which
is why extract_mip.py can use a fixed offset there.)

Output matches extract_mip.py's guide.tsv so the same downstream tools apply, but with
distinct_umi left as NA -- a PCR library has no molecule counter, which is the whole point
of the comparison.

Usage
    python3 extract_pcr.py R1.fq.gz --whitelist lib.csv --out-prefix out
"""
from __future__ import annotations

import argparse
import gzip
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "lib"))
sys.path.insert(0, str(HERE.parents[1] / "lib"))
from umimodel import gc_fraction  # noqa: E402
from extract_mip import load_whitelist, seqs  # noqa: E402

ANCHOR = "GTTTTAGAGC"
GUIDE_LEN = 20


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("r1", type=Path)
    ap.add_argument("--whitelist", type=Path, required=True)
    ap.add_argument("--anchor", default=ANCHOR)
    ap.add_argument("--out-prefix", type=Path, default=None)
    a = ap.parse_args()

    white = load_whitelist(a.whitelist)
    print(f"whitelist {len(white):,} guides")
    cnt: Counter[str] = Counter()
    tot = no_anchor = short = off_lib = 0
    pos = Counter()
    for s in seqs(a.r1):
        tot += 1
        i = s.find(a.anchor)
        if i < 0:
            no_anchor += 1
            continue
        if i < GUIDE_LEN:
            short += 1
            continue
        pos[i] += 1
        g = s[i - GUIDE_LEN:i]
        if g not in white:
            off_lib += 1
            continue
        cnt[g] += 1
    print(f"reads                 {tot:,}")
    print(f"  no scaffold anchor  {no_anchor:,} ({100 * no_anchor / max(1, tot):.1f}%)")
    print(f"  anchor too close to read start {short:,}")
    print(f"  guide not in library {off_lib:,} ({100 * off_lib / max(1, tot):.1f}%)")
    print(f"  kept                {sum(cnt.values()):,}")
    top = ", ".join(f"{p}:{n:,}" for p, n in pos.most_common(6))
    print(f"  anchor positions (confirms the stagger): {top}")
    print(f"guides observed       {len(cnt):,}")

    pref = a.out_prefix or a.r1.with_suffix("").with_suffix("")
    pg = Path(f"{pref}.guide.tsv")
    with pg.open("w") as fh:
        fh.write("guide\tgc\treads\tdistinct_umi\tm1\tm2\n")
        for g, n in sorted(cnt.items()):
            fh.write(f"{g}\t{gc_fraction(g):.4f}\t{n}\tNA\tNA\tNA\n")
    print(f"wrote {pg}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
