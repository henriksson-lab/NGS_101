#!/usr/bin/env python3
"""Paired-FASTQ extractor for CRISPR-MIP: (guide, UMI) read counts.

Different read layout from the Schmierer data that extract.py handles. There the UMI is an
Illumina index sitting in the FASTQ header; here both reads are real:

    R1[20:40] = the 20-nt sgRNA spacer
    R2[0:13]  = the 13-nt UMI

Those offsets are not arbitrary and are not taken on trust from the authors' script. They
follow from the construct in crispr-mip__10.1101+2024.03.28.587082/tools/crisprmip.py:

  * Read 1 primes on READ1_SITE, and what follows it in the library is
    EXT_ARM (19 nt, in U6) + the Pol III +1 G + the spacer. 19 + 1 = **20**.
  * Read 2 primes on READ2_SITE, and the probe backbone places the 13-nt UMI immediately
    adjacent, so the UMI starts at **0**. READ2_CYCLES is 15, i.e. 13 + 2 spare.

`--verify` re-derives the R1 offset empirically by scanning every offset for whitelist
matches, and refuses to run if the best offset is not the expected one. Use it on each new
dataset: a layout assumption that is wrong by a few bases fails silently, producing a
plausible-looking but near-empty table.
"""
from __future__ import annotations

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "lib"))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

from umimodel import gc_fraction  # noqa: E402

GUIDE_LEN = 20
GUIDE_OFFSET = 20          # ext arm (19) + Pol III +1 G
UMI_LEN = 13


def seqs(path: Path):
    op = gzip.open if path.suffix == ".gz" else open
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline().rstrip("\n")
            fh.readline()
            fh.readline()
            yield s


def load_whitelist(path: Path) -> set[str]:
    import re
    rows = [ln.rstrip("\n").split("\t") for ln in path.open() if ln.strip()]
    if len(rows[0]) < 2:
        rows = [ln.rstrip("\n").split(",") for ln in path.open() if ln.strip()]
    best, bestn = 0, -1
    ncol = max(len(r) for r in rows)
    for c in range(ncol):
        n = sum(1 for r in rows[1:] if len(r) > c
                and re.fullmatch(r"[ACGT]{%d}" % GUIDE_LEN, r[c].strip().upper()))
        if n > bestn:
            best, bestn = c, n
    return {r[best].strip().upper() for r in rows[1:]
            if len(r) > best and re.fullmatch(r"[ACGT]{%d}" % GUIDE_LEN, r[best].strip().upper())}


def verify_offset(r1: Path, white: set[str], n_probe: int = 200_000) -> int:
    """Scan offsets 0..40 and return the one with most whitelist hits."""
    hits = defaultdict(int)
    for i, s in enumerate(seqs(r1)):
        if i >= n_probe:
            break
        for off in range(0, min(41, max(0, len(s) - GUIDE_LEN) + 1)):
            if s[off:off + GUIDE_LEN] in white:
                hits[off] += 1
    return hits, max(hits, key=hits.get) if hits else -1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("r1", type=Path)
    ap.add_argument("r2", type=Path)
    ap.add_argument("--whitelist", type=Path, required=True)
    ap.add_argument("--guide-offset", type=int, default=GUIDE_OFFSET)
    ap.add_argument("--umi-len", type=int, default=UMI_LEN)
    ap.add_argument("--min-guide-reads", type=int, default=0)
    ap.add_argument("--verify", action="store_true",
                    help="re-derive the R1 offset from the data and abort on mismatch")
    ap.add_argument("--out-prefix", type=Path, default=None)
    a = ap.parse_args()

    white = load_whitelist(a.whitelist)
    print(f"whitelist {len(white):,} guides from {a.whitelist.name}")

    if a.verify:
        hits, best = verify_offset(a.r1, white)
        top = sorted(hits.items(), key=lambda kv: -kv[1])[:4]
        print("  offset scan (whitelist hits per offset, 200k reads): "
              + ", ".join(f"{o}:{n:,}" for o, n in top))
        if best != a.guide_offset:
            raise SystemExit(f"  ABORT: best offset is {best}, expected {a.guide_offset}. "
                             f"The read layout is not what this script assumes.")
        print(f"  verified: offset {best} is best, as the construct predicts")

    pair: dict[tuple[str, str], int] = defaultdict(int)
    tot = kept = bad_guide = bad_umi = 0
    for s1, s2 in zip(seqs(a.r1), seqs(a.r2)):
        tot += 1
        g = s1[a.guide_offset:a.guide_offset + GUIDE_LEN]
        u = s2[:a.umi_len]
        if len(u) < a.umi_len or "N" in u:
            bad_umi += 1
            continue
        if g not in white:
            bad_guide += 1
            continue
        pair[(g, u)] += 1
        kept += 1
    print(f"read pairs            {tot:,}")
    print(f"  dropped, UMI bad    {bad_umi:,}")
    print(f"  dropped, guide not in library {bad_guide:,} "
          f"({100 * bad_guide / max(1, tot):.1f}%)")
    print(f"  kept                {kept:,}")

    # Keep the real UMI sequence, not a placeholder: without it a later error-collapse
    # (the authors use UMIClusterer at edit distance 2) is impossible to reproduce.
    per: dict[str, dict[str, int]] = defaultdict(dict)
    for (g, u), n in pair.items():
        per[g][u] = n
    guides = {g: v for g, v in per.items() if sum(v.values()) >= a.min_guide_reads}
    print(f"guides observed       {len(per):,}  (>= {a.min_guide_reads} reads: {len(guides):,})")
    nu = sum(len(v) for v in guides.values())
    print(f"distinct (guide,UMI)  {nu:,}   mean reads/UMI "
          f"{sum(sum(v.values()) for v in guides.values()) / max(1, nu):.3f}")

    pref = a.out_prefix or a.r1.with_suffix("").with_suffix("")
    pu = Path(f"{pref}.guide_umi.tsv")
    with pu.open("w") as fh:
        fh.write("guide\tumi\treads\n")
        for g, v in sorted(guides.items()):
            for u, n in sorted(v.items(), key=lambda kv: -kv[1]):
                fh.write(f"{g}\t{u}\t{n}\n")
    pg = Path(f"{pref}.guide.tsv")
    with pg.open("w") as fh:
        fh.write("guide\tgc\treads\tdistinct_umi\tm1\tm2\n")
        for g, vd in sorted(guides.items()):
            v = list(vd.values())
            d, r = len(v), sum(v)
            fh.write(f"{g}\t{gc_fraction(g):.4f}\t{r}\t{d}\t{r / d:.6f}\t"
                     f"{sum(x * x for x in v) / d:.6f}\n")
    print(f"wrote {pu}\nwrote {pg}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
